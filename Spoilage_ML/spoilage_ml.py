#!/usr/bin/env python3
"""
spoilage_ml.py — Master CLI for early-vs-late spoilage modeling

Features:
- Label: early=1 if 'Total mesophilic aerobic flora (log10 CFU.g-1)' >= 7
- No leakage (drop label source column from X)
- Numeric-only predictors; optional microbe-only mode
- Models: L1-LogisticRegressionCV (default) or LassoCV (--model lasso)
- Handles imbalance (--balanced), custom threshold (--threshold)
- Artifacts: metrics.txt, confusion_matrix.txt, coefficients_*.csv, features_used.csv,
             threshold_sweep.csv, roc_curve.png (logit only),
             topN coeffs CSV/plot and odds-ratio CSV/plot (--topn N)
"""

import argparse
import pathlib
import sys
import numpy as np
import pandas as pd

# headless plotting
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegressionCV, LassoCV
from sklearn.metrics import (
    accuracy_score, roc_auc_score, classification_report, r2_score,
    confusion_matrix, precision_recall_fscore_support, RocCurveDisplay
)

TARGET_SUBSTR = "total mesophilic aerobic flora"
DEFAULT_DROP_CANDIDATES = {
    "samn_id", "srr_id", "sample_name", "sample_name_y", "envtype", "samplingtime"
}

# ---------- helpers ----------
def find_target_column(cols):
    lc = {c: c.lower() for c in cols}
    matches = [orig for orig, low in lc.items() if TARGET_SUBSTR in low]
    if not matches:
        raise ValueError(
            f"Could not find a column containing '{TARGET_SUBSTR}'. "
            f"Columns seen: {list(cols)[:12]}"
        )
    return matches[0]

def load_and_prepare(csv_path, drop_ids=True, microbe_only=False):
    df = pd.read_csv(csv_path)
    target_col = find_target_column(df.columns)

    # build label: >=7 -> 1 (early), else 0 (late)
    y = (df[target_col] >= 7).astype(int)

    # drop label-source col to prevent leakage
    X = df.drop(columns=[target_col], errors="ignore")

    # optionally drop obvious non-feature columns
    if drop_ids:
        to_drop = [c for c in X.columns if c.lower() in DEFAULT_DROP_CANDIDATES]
        if to_drop:
            X = X.drop(columns=to_drop, errors="ignore")

    # optionally drop time/day cols
    if microbe_only:
        for pattern in ["day_numeric", "day", "time", "samplingtime"]:
            tod = [c for c in X.columns if c.lower() == pattern]
            if tod:
                X = X.drop(columns=tod, errors="ignore")

    # numeric-only features
    X = X.select_dtypes(include=np.number)

    # drop rows with NaNs
    keep = (~X.isna().any(axis=1)) & (~y.isna())
    X = X.loc[keep]
    y = y.loc[keep]

    return X, y, target_col

def fit_logit(X, y, class_weight=None, random_state=100):
    pipe = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("clf", LogisticRegressionCV(
                Cs=10, cv=10, penalty="l1", solver="saga",
                scoring="roc_auc", max_iter=5000, refit=True,
                random_state=random_state, class_weight=class_weight
            ))
        ]
    )
    pipe.fit(X, y)
    clf = pipe.named_steps["clf"]
    return pipe, clf.coef_.ravel(), clf.intercept_[0], clf.C_[0]

def fit_lasso(X, y, random_state=100):
    pipe = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("reg", LassoCV(cv=10, random_state=random_state, max_iter=5000))
        ]
    )
    pipe.fit(X, y)
    reg = pipe.named_steps["reg"]
    return pipe, reg.coef_, reg.intercept_, reg.alpha_

def threshold_sweep_df(y_true, y_prob, thresholds):
    rows = []
    for t in thresholds:
        y_pred = (y_prob >= t).astype(int)
        acc = accuracy_score(y_true, y_pred)
        p, r, f1, _ = precision_recall_fscore_support(
            y_true, y_pred, average=None, labels=[0, 1], zero_division=0
        )
        rows.append({
            "threshold": float(t), "accuracy": float(acc),
            "precision_0": float(p[0]), "recall_0": float(r[0]), "f1_0": float(f1[0]),
            "precision_1": float(p[1]), "recall_1": float(r[1]), "f1_1": float(f1[1]),
        })
    return pd.DataFrame(rows)

def save_coef_artifacts(outdir, coef_df, topn):
    # full coefficients
    coef_df.to_csv(outdir / "coefficients_logit.csv", index=False)
    if topn and topn > 0:
        # make topN with odds ratios & plots
        top = coef_df.copy()
        top["abs_coef"] = top["Coefficient"].abs()
        top = top.sort_values("abs_coef", ascending=False).head(topn).drop(columns="abs_coef")
        top["odds_ratio_per_SD"] = np.exp(top["Coefficient"])
        top.to_csv(outdir / f"top{topn}_coeffs_with_OR.csv", index=False)

        # coef plot (signed)
        order = top.sort_values("Coefficient")
        plt.figure(figsize=(8, 0.35*len(order)+3))
        plt.barh(order["Feature"], order["Coefficient"])
        plt.xlabel("L1 Logistic Regression Coefficient (per 1 SD)")
        plt.ylabel("Feature")
        plt.title(f"Top {topn} Features by |Coefficient| — Early vs Late Spoilage")
        plt.tight_layout()
        plt.savefig(outdir / f"top{topn}_coeffs_bar.png", dpi=180)
        plt.close()

        # odds ratio plot (color by sign)
        order = top.sort_values("odds_ratio_per_SD")
        colors = ["#1f77b4" if c > 0 else "#d62728" for c in order["Coefficient"]]
        plt.figure(figsize=(8, 0.35*len(order)+3))
        plt.barh(order["Feature"], order["odds_ratio_per_SD"], color=colors)
        plt.axvline(1.0, color="black", linestyle="--")
        plt.xlabel("Odds Ratio (per 1 SD increase)")
        plt.ylabel("Feature")
        plt.title(f"Top {topn} Odds Ratios — Early vs Late Spoilage")
        plt.tight_layout()
        plt.savefig(outdir / f"top{topn}_odds_ratios.png", dpi=180)
        plt.close()

def save_roc(y_true, y_prob, auc, outpath_png):
    fig, ax = plt.subplots()
    RocCurveDisplay.from_predictions(y_true, y_prob, ax=ax, name=f"Logit (AUC={auc:.3f})")
    ax.set_title("ROC — Early vs Late Spoilage")
    fig.tight_layout()
    fig.savefig(outpath_png, dpi=180)
    plt.close(fig)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Early vs late spoilage modeling.")
    ap.add_argument("--csv", required=True, help="Path to fully_combined_project_data.csv")
    ap.add_argument("--outdir", default="results", help="Output directory (created if missing)")
    ap.add_argument("--model", choices=["logit", "lasso"], default="logit",
                    help="Logistic (default) or Lasso (regression on 0/1).")
    ap.add_argument("--test_size", type=float, default=0.25, help="Test fraction")
    ap.add_argument("--random_state", type=int, default=100, help="Random seed")
    ap.add_argument("--microbe_only", action="store_true", help="Drop time/day columns.")
    ap.add_argument("--balanced", action="store_true", help='Use class_weight="balanced".')
    ap.add_argument("--threshold", type=float, default=0.50, help="Decision cutoff for class 1.")
    ap.add_argument("--topn", type=int, default=20,
                    help="If >0 (default 20), save top-N coeff/OR CSVs and plots (logit only).")
    args = ap.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    X, y, target_col = load_and_prepare(args.csv, drop_ids=True, microbe_only=args.microbe_only)
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=args.test_size, random_state=args.random_state,
        shuffle=True, stratify=y
    )

    if args.model == "logit":
        class_weight = "balanced" if args.balanced else None
        pipe, coefs, intercept, best_C = fit_logit(
            X_train, y_train, class_weight=class_weight, random_state=args.random_state
        )
        y_prob = pipe.predict_proba(X_test)[:, 1]
        cutoff = float(args.threshold)
        y_pred = (y_prob >= cutoff).astype(int)

        acc = accuracy_score(y_test, y_pred)
        auc = roc_auc_score(y_test, y_prob)
        cm = confusion_matrix(y_test, y_pred, labels=[0, 1])

        # sweep
        ts = np.round(np.linspace(0.10, 0.90, 17), 2)
        sweep_df = threshold_sweep_df(y_test, y_prob, ts)
        sweep_df.to_csv(outdir / "threshold_sweep.csv", index=False)

        # roc
        save_roc(y_test, y_prob, auc, outdir / "roc_curve.png")

        # coefficients & topN
        coef_df = pd.DataFrame({"Feature": X.columns, "Coefficient": coefs})
        save_coef_artifacts(outdir, coef_df, args.topn)

        # confusion matrix + metrics
        np.savetxt(outdir / "confusion_matrix.txt", cm, fmt="%d")
        metrics_text = (
            f"Model: L1-LogisticRegressionCV\n"
            f"Target-from: '{target_col}' (>=7 => early)\n"
            f"Best C: {best_C}\n"
            f"Class weight: {'balanced' if class_weight else 'None'}\n"
            f"Threshold: {cutoff:.2f}\n"
            f"Accuracy (test): {acc:.3f}\n"
            f"ROC AUC (test): {auc:.3f}\n\n"
            f"Confusion matrix (rows=true, cols=pred) [0,1]:\n{cm}\n\n"
            + classification_report(y_test, y_pred, digits=3)
        )
        (outdir / "metrics.txt").write_text(metrics_text)
        print(metrics_text)

    else:  # lasso
        pipe, coefs, intercept, alpha = fit_lasso(X_train, y_train, random_state=args.random_state)
        y_pred_cont = pipe.predict(X_test)
        r2 = r2_score(y_test, y_pred_cont)

        coef_df = pd.DataFrame({"Feature": X.columns, "Coefficient": coefs})
        coef_df.to_csv(outdir / "coefficients_lasso.csv", index=False)

        metrics_text = (
            f"Model: LassoCV (regression on 0/1 label)\n"
            f"Target-from: '{target_col}' (>=7 => early)\n"
            f"Optimal alpha: {alpha}\n"
            f"R^2 (test): {r2:.3f}\n"
        )
        (outdir / "metrics.txt").write_text(metrics_text)
        print(metrics_text)

    # feature list
    pd.Series(X.columns, name="feature").to_csv(outdir / "features_used.csv", index=False)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        sys.stderr.write(f"[ERROR] {e}\n")
        sys.exit(1)
