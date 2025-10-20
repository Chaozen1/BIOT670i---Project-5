#!/usr/bin/env python3
"""
Spoilage ML Trainer (classification + CFU regression)

- Label: early/unsafe (1) if 'Total mesophilic aerobic flora (log10 CFU.g-1)' >= 7, else late/safe (0)
- Classifier: L1-LogisticRegressionCV (SAGA), adaptive StratifiedKFold
- CFU Regressor: ElasticNetCV on either:
    * ALL numeric features (default), or
    * ONLY 'late-associated' taxa (negative logit coefficients) + Day_numeric (+ interactions) with --late_only_reg
- Exports artifacts your app uses:
    - logit_model.pkl                 (pickle pipeline: scaler + logregCV)
    - clf_pipeline.joblib             (same classifier pipeline, joblib)
    - reg_pipeline.joblib             (CFU regressor pipeline)
    - model_meta.json                 (feature names, reg feature names, thresholds, etc.)
    - metrics.txt, confusion_matrix.txt, roc_curve.png
    - coefficients_logit.csv, topN plots: top{N}_coeffs_bar.png, top{N}_odds_ratios.png
    - threshold_sweep.csv
"""

import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegressionCV, ElasticNetCV, LassoCV, RidgeCV
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report, confusion_matrix, RocCurveDisplay, r2_score
import joblib, pickle

TARGET_SUBSTR = "total mesophilic aerobic flora"  # matches your header text
DROP_CANDS = {"samn_id","srr_id","sample_name","sample_name_y","envtype","samplingtime"}

# ---------- helpers ----------
def find_target_column(cols):
    for c in cols:
        if TARGET_SUBSTR in c.lower():
            return c
    raise ValueError(f"No column containing '{TARGET_SUBSTR}' was found.")

def load_and_prepare(csv_path, microbe_only=False):
    df = pd.read_csv(csv_path)
    target_col = find_target_column(df.columns)

    # label: unsafe/early = 1 if CFU >= 7
    y = (df[target_col] >= 7).astype(int)
    y_cont = df[target_col].copy()

    # drop the target from features
    X = df.drop(columns=[target_col], errors="ignore")

    # keep Day_numeric as numeric if present
    if "Day_numeric" in df.columns:
        X["Day_numeric"] = pd.to_numeric(df["Day_numeric"], errors="coerce")

    # optionally remove Day/Time etc if microbe-only requested
    if microbe_only:
        tod = [c for c in X.columns if "day" in c.lower() or "time" in c.lower()]
        X = X.drop(columns=tod, errors="ignore")

    # drop obvious IDs / strings
    to_drop = [c for c in X.columns if c.lower() in DROP_CANDS]
    if to_drop:
        X = X.drop(columns=to_drop, errors="ignore")

    # numeric only
    X = X.select_dtypes(include=np.number)

    # force-restore Day_numeric if it existed and was numeric
    if "Day_numeric" in df.columns and "Day_numeric" not in X.columns:
        maybe = pd.to_numeric(df["Day_numeric"], errors="coerce")
        if np.issubdtype(maybe.dtype, np.number):
            X = pd.concat([X, maybe.rename("Day_numeric")], axis=1)

    # drop rows with NaNs
    keep = (~X.isna().any(axis=1)) & (~y.isna())
    X = X.loc[keep].reset_index(drop=True)
    y = y.loc[keep].reset_index(drop=True)
    y_cont = y_cont.loc[keep].reset_index(drop=True)

    return df, X, y, y_cont, target_col

def adaptive_cv(y, seed=100):
    counts = np.bincount(np.asarray(y))
    min_class = int(counts.min()) if counts.size > 1 else len(y)
    n_splits = max(3, min(10, min_class))
    return StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

def fit_logit(X, y, class_weight=None, seed=100):
    cv = adaptive_cv(y, seed)
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegressionCV(
            Cs=10, cv=cv, penalty="l1", solver="saga", scoring="roc_auc",
            max_iter=5000, refit=True, random_state=seed, class_weight=class_weight
        ))
    ])
    pipe.fit(X, y)
    clf = pipe.named_steps["clf"]
    return pipe, clf.coef_.ravel(), clf.intercept_[0], float(clf.C_[0])

def save_roc(y_true, y_prob, auc, out_png):
    fig, ax = plt.subplots()
    RocCurveDisplay.from_predictions(y_true, y_prob, ax=ax, name=f"AUC={auc:.3f}")
    ax.set_title("ROC — Early (unsafe=1) vs Late (safe=0)")
    fig.tight_layout(); fig.savefig(out_png, dpi=180); plt.close(fig)

def save_coef_artifacts(outdir: Path, coef_df: pd.DataFrame, topn: int):
    (outdir / "coefficients_logit.csv").write_text(coef_df.to_csv(index=False))
    if topn and topn > 0:
        top = coef_df.assign(abs_coef=lambda d: d["Coefficient"].abs())\
                     .sort_values("abs_coef", ascending=False)\
                     .head(topn).drop(columns="abs_coef")
        top["odds_ratio_per_SD"] = np.exp(top["Coefficient"])
        top.to_csv(outdir / f"top{topn}_coeffs_with_OR.csv", index=False)

        # coef plot
        order = top.sort_values("Coefficient")
        plt.figure(figsize=(8, 0.35*len(order)+3))
        plt.barh(order["Feature"], order["Coefficient"])
        plt.xlabel("Coefficient (per 1 SD)"); plt.title(f"Top {topn} Coefficients")
        plt.tight_layout(); plt.savefig(outdir / f"top{topn}_coeffs_bar.png", dpi=180); plt.close()

        # OR plot
        order = top.sort_values("odds_ratio_per_SD")
        colors = ["#1f77b4" if c>0 else "#d62728" for c in order["Coefficient"]]
        plt.figure(figsize=(8, 0.35*len(order)+3))
        plt.barh(order["Feature"], order["odds_ratio_per_SD"], color=colors)
        plt.axvline(1.0, color="black", linestyle="--")
        plt.xlabel("Odds Ratio (per 1 SD)"); plt.title(f"Top {topn} Odds Ratios")
        plt.tight_layout(); plt.savefig(outdir / f"top{topn}_odds_ratios.png", dpi=180); plt.close()

def threshold_sweep_df(y_true, y_prob, thresholds):
    rows = []
    from sklearn.metrics import precision_recall_fscore_support, accuracy_score
    for t in thresholds:
        yp = (y_prob >= t).astype(int)
        acc = accuracy_score(y_true, yp)
        p,r,f1,_ = precision_recall_fscore_support(y_true, yp, average=None, labels=[0,1], zero_division=0)
        rows.append({"threshold": float(t), "accuracy": float(acc),
                     "precision_0": float(p[0]), "recall_0": float(r[0]), "f1_0": float(f1[0]),
                     "precision_1": float(p[1]), "recall_1": float(r[1]), "f1_1": float(f1[1])})
    return pd.DataFrame(rows)

# ---- Late-only CFU regressor helpers ----
from sklearn.preprocessing import StandardScaler

def get_late_features_from_logit(coef_df: pd.DataFrame, topk: int, include_day: bool, X_cols):
    late = coef_df[coef_df["Coefficient"] < 0].copy()  # negative → safer/late-associated
    late["abs"] = late["Coefficient"].abs()
    feats = late.sort_values("abs", ascending=False).head(topk)["Feature"].tolist()
    if include_day and "Day_numeric" in X_cols and "Day_numeric" not in feats:
        feats.append("Day_numeric")
    return [f for f in feats if f in X_cols]

def add_day_interactions(X: pd.DataFrame):
    """Add (feature * Day_numeric) interaction columns."""
    if "Day_numeric" not in X.columns:
        return X, []
    others = X.drop(columns=["Day_numeric"])
    if others.shape[1] == 0:
        return X, []
    inter = others.mul(X["Day_numeric"], axis=0)
    inter.columns = [f"{c}__x__Day" for c in others.columns]
    X2 = pd.concat([X, inter], axis=1)
    return X2, inter.columns.tolist()

def fit_cfu_regressor_smart(X, y_cont, use_features=None, l1_ratio=0.2):
    # choose subset for regressor
    Xr = X[use_features].copy() if use_features is not None else X.copy()
    # add interactions with Day_numeric
    Xr, inter_cols = add_day_interactions(Xr)
    pipe = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("reg", ElasticNetCV(l1_ratio=l1_ratio, n_alphas=100, cv=5, random_state=42, max_iter=5000))
    ])
    pipe.fit(Xr, y_cont)
    en = pipe.named_steps["reg"]
    return pipe, Xr.columns.tolist(), float(en.alpha_), float(en.l1_ratio)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--test_size", type=float, default=0.25)
    ap.add_argument("--random_state", type=int, default=100)
    ap.add_argument("--microbe_only", action="store_true")
    ap.add_argument("--balanced", action="store_true")
    ap.add_argument("--threshold", type=float, default=0.40)
    ap.add_argument("--topn", type=int, default=20)
    # late-only CFU regressor options
    ap.add_argument("--late_only_reg", action="store_true",
                    help="Use only late-associated taxa (negative logit coefs) + Day_numeric (+ interactions) for CFU regressor.")
    ap.add_argument("--late_topk", type=int, default=30,
                    help="Top-K (by |neg logit coef|) to keep for late-only regressor.")
    ap.add_argument("--reg_alpha_l1", type=float, default=0.2,
                    help="ElasticNet l1_ratio for CFU regressor (0..1).")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    df, X, y, y_cont, target_col = load_and_prepare(args.csv, microbe_only=args.microbe_only)
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=args.test_size,
                                          random_state=args.random_state, stratify=y)

    cw = "balanced" if args.balanced else None
    # --- classifier ---
    clf_pipe, coefs, intercept, bestC = fit_logit(Xtr, ytr, class_weight=cw, seed=args.random_state)
    yprob = clf_pipe.predict_proba(Xte)[:, 1]
    ypred = (yprob >= args.threshold).astype(int)

    acc = accuracy_score(yte, ypred)
    auc = roc_auc_score(yte, yprob)
    cm = confusion_matrix(yte, ypred, labels=[0,1])

    # save classifier artifacts
    coef_df = pd.DataFrame({"Feature": X.columns, "Coefficient": coefs})
    save_coef_artifacts(outdir, coef_df, args.topn)
    save_roc(yte, yprob, auc, outdir / "roc_curve.png")
    threshold_sweep_df(yte, yprob, np.round(np.linspace(0.10, 0.90, 17), 2))\
        .to_csv(outdir / "threshold_sweep.csv", index=False)
    np.savetxt(outdir / "confusion_matrix.txt", cm, fmt="%d")

    metrics_text = (
        f"Model: L1-Logit\n"
        f"Target: {target_col} (>=7 early/unsafe=1)\n"
        f"Best C:{bestC}\n"
        f"Class weight:{cw}\n"
        f"Threshold:{args.threshold}\n"
        f"Acc:{acc:.3f}\nAUC:{auc:.3f}\n"
        f"CM:\n{cm}\n\n" + classification_report(yte, ypred, digits=3)
    )
    (outdir / "metrics.txt").write_text(metrics_text)
    print(metrics_text)

    # also save as pkl + joblib for the app
    joblib.dump(clf_pipe, outdir / "clf_pipeline.joblib")
    with open(outdir / "logit_model.pkl", "wb") as f:
        pickle.dump(clf_pipe, f)

    # --- CFU regressor ---
    # choose features for regressor
    if args.late_only_reg:
        chosen_feats = get_late_features_from_logit(coef_df, args.late_topk, include_day=True, X_cols=X.columns)
    else:
        chosen_feats = list(X.columns)

    reg_pipe, used_cols, en_alpha, en_l1 = fit_cfu_regressor_smart(
        Xtr, y_cont.loc[Xtr.index], use_features=chosen_feats, l1_ratio=args.reg_alpha_l1
    )
    joblib.dump(reg_pipe, outdir / "reg_pipeline.joblib")

    # fit and print simple R2 on test split for info
    # NOTE: must transform test columns same as training set (select used_cols + interactions added inside pipeline)
    # We simply compute predictions on test with reg_pipe:
    r2 = r2_score(y_cont.loc[Xte.index], reg_pipe.predict(
        # Rebuild Xte view with columns used during training (missing ones will be filled by pipeline scaler)
        pd.concat([Xte, pd.DataFrame(index=Xte.index)], axis=1)[used_cols] if all(c in Xte.columns for c in used_cols)
        else Xte.reindex(columns=used_cols, fill_value=0.0)
    ))
    (outdir / "reg_r2.txt").write_text(f"R2 on test (naive): {r2:.3f}\n")

    # --- meta for app ---
    meta = {
        "target_col": target_col,
        "threshold_cfu": 7.0,
        "best_C": bestC,
        "class_weight": cw,
        "threshold": float(args.threshold),
        "feature_names": list(X.columns),           # classifier features (numeric-only)
        "reg_feature_names": used_cols,             # features actually used by CFU regressor (may include interactions)
        "microbe_only": bool(args.microbe_only),
        "late_only_reg": bool(args.late_only_reg),
        "late_topk": int(args.late_topk),
        "reg_en_alpha": en_alpha,
        "reg_en_l1_ratio": en_l1
    }
    (outdir / "model_meta.json").write_text(json.dumps(meta, indent=2))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[ERROR] {e}\n")
        sys.exit(1)
