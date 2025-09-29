#!/usr/bin/env python3
"""
Spoilage 16S pipeline (matplotlib only; no seaborn / scikit-bio)

Inputs (CSV):
  --features : Cleaned_TableS09_16S_Wide.csv  (rows=Cluster_ID, cols=samples)
  --taxonomy : Cleaned_TableS09_Taxonomy.csv  (Cluster_ID + Domain..Genus)
  --metadata : Cleaned_TableS09_Sample_Metadata_withDay.csv
               (Sample_Name, EnvType, Packaging, Day_numeric, ...)
  --sample-id-col : Sample_Name
  --day-col       : Day_numeric
  --group-col     : EnvType  (Pork/Poultry, …)

Outputs (in --outdir):
  alpha_richness.png, pca.png
  microbes_by_time/stacked_top12_{Env}.png + TSVs
  timecourse_all_days/timecourse_top.png + TSV
  timecourse_all_days/diffbar_top.png + TSV   <-- NEW (Δ vs Day 2)
  packaging_vs_spoilage/{Genus}_vs_Packaging.png + packaging_summary.tsv (with p)
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------- utils ----------------

def _ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def load_tables(features_fp, taxonomy_fp, metadata_fp, sample_id_col):
    ft = pd.read_csv(features_fp)
    tax = pd.read_csv(taxonomy_fp)
    md  = pd.read_csv(metadata_fp)

    if ft.columns[0] != "Cluster_ID":
        raise ValueError("First column in features must be 'Cluster_ID'")
    ft = ft.set_index("Cluster_ID")

    if tax.columns[0] != "Cluster_ID":
        raise ValueError("First column in taxonomy must be 'Cluster_ID'")
    tax = tax.set_index("Cluster_ID")

    if sample_id_col not in md.columns:
        raise ValueError(f"{sample_id_col} not found in metadata")
    md = md.set_index(sample_id_col).sort_index()

    keep = [c for c in ft.columns if c in md.index]
    ft = ft.loc[:, keep]
    md = md.loc[keep]
    return ft, tax, md

def to_genus_rel(ft_counts: pd.DataFrame, tax: pd.DataFrame, md_index):
    """Collapse Cluster_ID -> Genus and return relative abundance (rows=Genus, cols=Sample)."""
    ft_num = ft_counts.select_dtypes(include=[np.number])
    common = ft_num.index.intersection(tax.index)
    if len(common) == 0:
        raise ValueError("No overlap between features and taxonomy Cluster_IDs")
    ft_num = ft_num.loc[common]
    if "Genus" not in tax.columns:
        raise ValueError("'Genus' column missing in taxonomy")

    gmat = ft_num.groupby(tax.loc[common, "Genus"]).sum()
    gmat = gmat.loc[:, gmat.columns.intersection(md_index)]
    rel  = gmat.div(gmat.sum(axis=0).replace(0, np.nan), axis=1).fillna(0.0)
    return rel

# -------------- basic QC --------------

def plot_alpha_richness(ft_counts, md, group_col, outdir: Path):
    rich = (ft_counts > 0).sum(axis=0)
    md = md.copy(); md["richness"] = md.index.map(rich.to_dict())

    plt.figure(figsize=(7,5))
    for g, sub in md.groupby(group_col):
        vals = sub["richness"].dropna().values
        if len(vals) == 0: continue
        plt.hist(vals, bins=20, alpha=0.5, label=str(g))
    plt.xlabel("Observed clusters (richness)"); plt.ylabel("Sample count")
    plt.title("Richness by group"); plt.legend(); plt.tight_layout()
    plt.savefig(outdir / "alpha_richness.png", dpi=160); plt.close()

def plot_pca(ft_counts, md, group_col, outdir: Path):
    if ft_counts.shape[1] < 2: return
    X = np.log1p(ft_counts.T.values)
    X = (X - X.mean(0)) / (X.std(0) + 1e-8)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    pc = U[:, :2] * S[:2]
    md = md.copy(); md["PC1"], md["PC2"] = pc[:,0], pc[:,1]

    plt.figure(figsize=(6,6))
    for g, sub in md.groupby(group_col):
        plt.scatter(sub["PC1"], sub["PC2"], s=20, label=str(g))
    plt.axhline(0, color="k", lw=0.5); plt.axvline(0, color="k", lw=0.5)
    plt.xlabel("PC1"); plt.ylabel("PC2"); plt.title("PCA of log counts")
    plt.legend(); plt.tight_layout(); plt.savefig(outdir / "pca.png", dpi=160); plt.close()

# ----------- microbe-centric -----------

def plot_microbes_by_time(ft_counts, tax, md, day_col, group_col, outdir: Path, top_n=12):
    out = outdir / "microbes_by_time"; _ensure_dir(out)
    genus_rel = to_genus_rel(ft_counts, tax, md.index)

    for env in md[group_col].dropna().unique():
        keep_md = md[md[group_col] == env]
        if keep_md.empty: continue
        rel = genus_rel.loc[:, genus_rel.columns.intersection(keep_md.index)].T
        if rel.empty: continue
        rel["Day"] = pd.to_numeric(keep_md[day_col], errors="coerce").values
        rel = rel[rel["Day"].isin([2,7,15,22])]
        if rel.empty: continue

        mean_by_day = rel.groupby("Day").mean().sort_index()   # day x genus
        top = mean_by_day.mean(0).sort_values(ascending=False).head(top_n).index
        ax = mean_by_day.loc[:, top].plot(kind="bar", stacked=True, figsize=(11,6))
        ax.set_ylabel("Relative abundance"); ax.set_xlabel("Day")
        ax.set_title(f"Top {top_n} genera in {env} by day")
        ax.legend(title="Genus", bbox_to_anchor=(1.02,1), loc="upper left")
        plt.tight_layout(); plt.savefig(out / f"stacked_top{top_n}_{env}.png", dpi=160); plt.close()
        mean_by_day.to_csv(out / f"mean_genus_rel_by_day_{env}.tsv", sep="\t")

def plot_timecourse_all(ft_counts, tax, md, day_col, outdir: Path,
                        spoilage_list=None, top_n=12):
    out = outdir / "timecourse_all_days"; _ensure_dir(out)
    genus_rel = to_genus_rel(ft_counts, tax, md.index)
    days = pd.to_numeric(md[day_col], errors="coerce")
    df = genus_rel.T.copy(); df["Day"] = days.values
    df = df[df["Day"].isin([2,7,15,22])]
    if df.empty:
        (out / "NOTE_no_requested_days.txt").write_text("No samples at days 2/7/15/22.")
        return
    mean_by_day = df.groupby("Day").mean().sort_index()  # day x genus

    # choose genera to show
    if spoilage_list:
        show = [g for g in spoilage_list if g in mean_by_day.columns]
        if not show:
            show = mean_by_day.mean().sort_values(ascending=False).head(top_n).index.tolist()
    else:
        show = mean_by_day.mean().sort_values(ascending=False).head(top_n).index.tolist()

    # line plot (absolute mean rel)
    plt.figure(figsize=(11,6))
    for g in show:
        plt.plot(mean_by_day.index.values, mean_by_day[g].values, marker="o", label=g)
    plt.xticks([2,7,15,22]); plt.xlabel("Day"); plt.ylabel("Mean relative abundance")
    plt.title("Genus time-course (Days 2, 7, 15, 22)")
    plt.legend(title="Genus", bbox_to_anchor=(1.02,1), loc="upper left")
    plt.tight_layout(); plt.savefig(out / "timecourse_top.png", dpi=160); plt.close()

    # --- NEW: differential bar (Δ vs Day 2 baseline) ---
    if 2 in mean_by_day.index:
        base = mean_by_day.loc[2]
        diff = mean_by_day[show].subtract(base[show], axis=1)  # day x genus (Δ vs day2)
        diff = diff.loc[[2,7,15,22].__iter__()]  # keep order if all present
        diff = diff.dropna(axis=0, how="all")     # just in case

        fig, ax = plt.subplots(figsize=(12,6))
        idx = np.arange(len(diff.index))
        width = max(0.75 / len(show), 0.02)
        for i, g in enumerate(show):
            ax.bar(idx + i*width, diff[g].values, width=width, label=g)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(idx + (len(show)-1)*width/2)
        ax.set_xticklabels([str(d) for d in diff.index])
        ax.set_xlabel("Day"); ax.set_ylabel("Δ mean rel. abundance vs Day 2")
        ax.set_title("Differential abundance vs Day 2 (top spoilage genera)")
        ax.legend(title="Genus", bbox_to_anchor=(1.02,1), loc="upper left")
        plt.tight_layout(); plt.savefig(out / "diffbar_top.png", dpi=160); plt.close()

    mean_by_day.to_csv(out / "timecourse_mean_rel_by_day.tsv", sep="\t")

def plot_packaging(ft_counts, tax, md, outdir: Path):
    out = outdir / "packaging_vs_spoilage"; _ensure_dir(out)
    if "Packaging" not in md.columns:
        (out / "NOTE_missing_packaging.txt").write_text("No 'Packaging' column in metadata.")
        return

    genus_rel = to_genus_rel(ft_counts, tax, md.index)
    targets = [
        "Brochothrix","Lactobacillus","Leuconostoc",
        "Pseudomonas","Carnobacterium","Shewanella","Photobacterium"
    ]
    pkg = md["Packaging"].astype(str)
    rows = []

    # optional stats
    try:
        from scipy.stats import kruskal
    except Exception:
        kruskal = None

    for g in targets:
        if g not in genus_rel.index: continue
        sidx = genus_rel.columns.intersection(pkg.index)
        vals = pd.DataFrame({"rel": genus_rel.loc[g, sidx].values,
                             "Packaging": pkg.loc[sidx].values}).dropna()
        if vals.empty: continue

        labels = sorted(vals["Packaging"].unique())
        groups = [vals.loc[vals["Packaging"]==lab, "rel"].values for lab in labels]
        medians = [np.median(x) if len(x)>0 else np.nan for x in groups]

        # stats
        if kruskal and all(len(x)>0 for x in groups) and len(groups) >= 2:
            try:
                stat, p = kruskal(*groups)
            except Exception:
                stat, p = np.nan, np.nan
        else:
            stat, p = np.nan, np.nan

        # plot
        plt.figure(figsize=(7,5))
        bp = plt.boxplot(groups, labels=labels)
        plt.ylabel("Relative abundance"); plt.xlabel("Packaging")
        title = f"{g} relative abundance vs Packaging"
        if not np.isnan(p):
            title += f"  (Kruskal p={p:.3g})"
        plt.title(title)
        for i, m in enumerate(medians, start=1):
            plt.text(i, m, f"median={m:.3f}", ha="center", va="bottom", fontsize=9)
        plt.tight_layout(); plt.savefig(out / f"{g}_vs_Packaging.png", dpi=160); plt.close()

        for lab, m in zip(labels, medians):
            rows.append({"Genus": g, "Packaging": lab, "median_rel": m,
                         "Kruskal_H": stat, "p_value": p})

    pd.DataFrame(rows).to_csv(out / "packaging_summary.tsv", sep="\t", index=False)

# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser(description="Spoilage 16S microbe-centric EDA (matplotlib only)")
    ap.add_argument("--features", required=True)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--sample-id-col", required=True)
    ap.add_argument("--day-col", required=True)
    ap.add_argument("--group-col", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir); _ensure_dir(outdir)

    ft_counts, tax, md = load_tables(args.features, args.taxonomy, args.metadata, args.sample_id_col)
    print(f"[align] samples={ft_counts.shape[1]} features={ft_counts.shape[0]}")

    plot_alpha_richness(ft_counts, md, args.group_col, outdir)
    plot_pca(ft_counts, md, args.group_col, outdir)

    plot_microbes_by_time(ft_counts, tax, md, args.day_col, args.group_col, outdir, top_n=12)
    plot_timecourse_all(
        ft_counts, tax, md, args.day_col, outdir,
        spoilage_list=["Brochothrix","Lactobacillus","Leuconostoc",
                       "Pseudomonas","Carnobacterium","Shewanella","Photobacterium"],
        top_n=12
    )
    plot_packaging(ft_counts, tax, md, outdir)

    print("DONE:", outdir.resolve())

if __name__ == "__main__":
    main()
