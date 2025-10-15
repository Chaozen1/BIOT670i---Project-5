#!/usr/bin/env python3
import json, pickle
from pathlib import Path
from io import BytesIO

import joblib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st
import plotly.express as px  # Plotly for Microbiome + Interpretations

# ---------------- UI / THEME ----------------
st.set_page_config(page_title="Spoilage Predictor", page_icon="🥓", layout="wide")
st.markdown("""
<style>
  .stApp { background: linear-gradient(180deg, #ffeef8 0%, #fff9fc 100%); color:#4d004d; }
  h1,h2,h3 { color:#e75480; }
  .stButton button{ background:#ffb6c1; color:#4d004d; border:none; border-radius:12px;
    padding:.55rem 1.2rem; font-weight:700; box-shadow:0 1px 8px rgba(231,84,128,.25);}
  .stButton button:hover{ background:#ffc8d9; color:#000; }
</style>
""", unsafe_allow_html=True)

st.title("🥩🍖 Microbe-Driven Spoilage Predictor 🍕🥘")
st.caption("Binary Safe/Unsafe uses **probability**.")

# ---------------- SIDEBAR SETTINGS ----------------
with st.sidebar:
    st.header("⚙️ Settings")
    default_artifacts = str(Path.home() / "Desktop" / "results_with_day")
    artifacts_dir_str = st.text_input(
        "Artifacts folder", value=default_artifacts,
        help="Needs: logit_model.pkl, reg_pipeline.joblib, model_meta.json"
    )
    prob_thr = st.slider("Unsafe probability threshold", 0.10, 0.90, 0.40, step=0.01)

    st.markdown(
        """
        ---
        **Disclaimer (RUO)**  
        This tool is for **Research Use Only** — not for clinical/regulatory decisions.  
        Models were trained on **poultry & pork** datasets and may not generalize to all products.  
        Exploring any cutoff is fine, but **0.40–0.50** is a sensible range for safety interpretation.
        """,
        unsafe_allow_html=True
    )

artifacts_dir = Path(artifacts_dir_str).expanduser()
need = ["logit_model.pkl", "reg_pipeline.joblib", "model_meta.json"]
if not artifacts_dir.exists() or any(not (artifacts_dir / f).exists() for f in need):
    st.warning(f"Pick a valid artifacts folder containing: {need}")
    st.stop()

@st.cache_resource
def load_artifacts(folder: Path):
    with open(folder / "logit_model.pkl", "rb") as f:
        clf = pickle.load(f)                           # classifier
    reg = joblib.load(folder / "reg_pipeline.joblib")  # CFU regression
    meta = json.loads((folder / "model_meta.json").read_text())
    return clf, reg, meta

try:
    clf, reg, meta = load_artifacts(artifacts_dir)
except Exception as e:
    st.error(f"Failed to load artifacts: {e}")
    st.stop()

feature_names = meta.get("feature_names", [])
cfu_thr = float(meta.get("threshold_cfu", 7.0))

# ---------------- DATA INPUT ----------------
csv_file = st.file_uploader("Upload CSV (same schema as training)", type=["csv"])
if csv_file is None:
    st.info("Upload a CSV to get predictions and microbiome views.")
    st.stop()

try:
    df_raw = pd.read_csv(csv_file)
except Exception as e:
    st.error(f"Could not read CSV: {e}")
    st.stop()

# ---------------- CORE PREDICTIONS (compute once) ----------------
# Align to training features; coerce numerics
X = df_raw.reindex(columns=feature_names, fill_value=0)
for c in X.columns:
    if not np.issubdtype(X[c].dtype, np.number):
        X[c] = pd.to_numeric(X[c], errors="coerce")
X = X.fillna(0.0)

# 1) Probability & label (always)
prob_unsafe = clf.predict_proba(X)[:, 1]
pred_class = np.where(prob_unsafe >= prob_thr, "Unsafe", "Safe")

# 2) Regressor outputs (optional)
regressor_ok = True
cfu_today = None
try:
    cfu_today = reg.predict(X)
except Exception:
    regressor_ok = False

# 3) Column mapping for display
ITEM_COL = "Sample_Name_y"
PTYPE_COL = "EnvType"
DAY_COL = "Day_numeric" if "Day_numeric" in df_raw.columns else ("Days_Numeric" if "Days_Numeric" in df_raw.columns else None)
if not DAY_COL:
    st.error("Neither 'Day_numeric' nor 'Days_Numeric' found in the CSV.")
    st.stop()

# Pastel palette for Plotly charts
PALETTE = ["#d5b4ab", "#c3b1e1", "#a3c1da", "#a8c69f", "#f4c2c2"]

# ---------------- TABS ----------------
tab_pred, tab_viz, tab_micro, tab_interp, tab_refs = st.tabs(
    ["🔮 Predictions", "📈 Visuals", "🧬 Microbiome", "📚 Interpretations", "📖 References"]
)

# ====== PREDICTIONS ======
with tab_pred:
    if not regressor_ok:
        st.info("CFU model disabled this run (feature mismatch). Showing probabilities only.")

    disp = pd.DataFrame({
        "Days in Refrigerator": df_raw[DAY_COL],
        "Item": df_raw[ITEM_COL] if ITEM_COL in df_raw.columns else "",
        "Product Type": df_raw[PTYPE_COL] if PTYPE_COL in df_raw.columns else "",
        "Safe/Unsafe": pred_class,
        "Probability_UnsafeToday": np.round(prob_unsafe, 3),
    })
    if regressor_ok and cfu_today is not None:
        disp["log CFU (pred today)"] = np.round(cfu_today, 2)

    st.markdown("### 📊 Predictions")
    st.dataframe(disp, use_container_width=True)

    st.download_button(
        "⬇️ Download predictions as CSV",
        disp.to_csv(index=False).encode("utf-8"),
        "spoilage_predictions.csv",
        "text/csv"
    )

# ====== VISUALS ======
with tab_viz:
    c1, c2 = st.columns(2)
    with c1:
        fig1, ax1 = plt.subplots(figsize=(5.8,3.8))
        ax1.hist(prob_unsafe, bins=20, edgecolor="white")
        ax1.set_xlabel("Probability Unsafe Today")
        ax1.set_ylabel("Count")
        ax1.set_title("Risk score distribution")
        st.pyplot(fig1)

    with c2:
        if regressor_ok and cfu_today is not None:
            fig2, ax2 = plt.subplots(figsize=(5.8,3.8))
            ax2.hist(cfu_today, bins=20, edgecolor="white")
            ax2.set_xlabel("log CFU (pred today)")
            ax2.set_ylabel("Count")
            ax2.set_title("Predicted microbial load")
            st.pyplot(fig2)
        else:
            st.write("CFU distribution not available.")

# ====== MICROBIOME (PLOTLY) ======
with tab_micro:
    st.markdown("### Relative abundance views")

    NON_MICROBE_COLS = {
        "Sample_Name_y","EnvType","EnvTyp",
        "Day_numeric","Days_Numeric","Day","Days",
        "SRR_ID","SAMN_ID","SamplingTime",
        "Total mesophilic aerobic flora (log10 CFU·g-1)",
        # predictions / outputs
        "Safe/Unsafe","Probability_UnsafeToday","Days Until Spoiled","log CFU (pred today)"
    }

    # Numeric columns not in NON_MICROBE_COLS are treated as microbes
    micro_cols = [
        c for c in df_raw.columns
        if c not in NON_MICROBE_COLS and np.issubdtype(df_raw[c].dtype, np.number)
    ]

    if not micro_cols:
        st.info("No numeric microbial columns detected to compute relative abundances.")
    else:
        micro_df = df_raw[micro_cols].copy()
        row_totals = micro_df.sum(axis=1).replace(0, np.nan)
        rel_abund = micro_df.div(row_totals, axis=0).fillna(0.0)  # 0..1

        # --- Single-sample Top-N (Plotly)
        left, right = st.columns([2, 1], gap="large")
        with left:
            sample_labels = df_raw["Sample_Name_y"].astype(str) if "Sample_Name_y" in df_raw.columns else df_raw.index.astype(str)
            default_idx = 0
            selected_label = st.selectbox("Select a sample", options=list(sample_labels), index=default_idx, key="sample_select")

            # resolve row index + title
            if "Sample_Name_y" in df_raw.columns:
                row_idx = df_raw.index[df_raw["Sample_Name_y"].astype(str) == selected_label][0]
                title_item = selected_label
            else:
                row_idx = int(selected_label)
                title_item = f"Row {row_idx}"

            topN = st.slider("Top microbes (N)", 5, 20, 10, 1, key="topN_single")
            series = (rel_abund.iloc[row_idx].sort_values(ascending=False).head(topN) * 100.0)
            df_top = series.rename("Relative_Abundance_%").reset_index().rename(columns={"index": "Microbe"})

            fig = px.bar(
                df_top, x="Microbe", y="Relative_Abundance_%",
                color="Microbe",
                color_discrete_sequence=(PALETTE * ((len(df_top)//len(PALETTE))+1))[:len(df_top)],
                title=f"Top {topN} microbes — {title_item}"
            )
            fig.update_layout(yaxis_title="Relative abundance (%)", xaxis_title="")
            st.plotly_chart(fig, use_container_width=True)

        with right:
            st.write("Top microbes (this sample)")
            st.dataframe(df_top, use_container_width=True, width=500)

        st.markdown("---")

        # ----- ⏱️ Time-course line chart (inside Microbiome tab only) -----
        with st.expander("⏱️ Genus time-course (mean relative abundance by day)"):
            if "Day_numeric" not in df_raw.columns and "Days_Numeric" not in df_raw.columns:
                st.info("No day column found (Day_numeric/Days_Numeric); cannot build time-course.")
            else:
                day_col = "Day_numeric" if "Day_numeric" in df_raw.columns else "Days_Numeric"
                rel_with_day = rel_abund.copy()
                rel_with_day[day_col] = pd.to_numeric(df_raw[day_col], errors="coerce").fillna(0).astype(int)

                overall_means = rel_with_day[micro_cols].mean().sort_values(ascending=False)
                taxa_all = list(overall_means.index)
                curated = ["Brochothrix","Lactobacillus","Leuconostoc","Pseudomonas",
                           "Carnobacterium","Shewanella","Photobacterium"]

                left_tc, right_tc = st.columns([2,1])
                with right_tc:
                    use_curated = st.checkbox("Use curated genus list", value=True,
                                              help="Brochothrix, Lactobacillus, Leuconostoc, Pseudomonas, Carnobacterium, Shewanella, Photobacterium",
                                              key="tc_use_curated")
                    k = st.slider("How many genera to show", 3, 12, 7, 1, key="tc_k")

                if use_curated:
                    taxa_top = [t for t in curated if t in rel_abund.columns][:k]
                    if not taxa_top:
                        taxa_top = taxa_all[:k]
                else:
                    taxa_top = taxa_all[:k]

                by_day = rel_with_day.groupby(day_col)[taxa_top].mean() * 100.0
                df_time = by_day.reset_index().melt(id_vars=day_col, var_name="Genus", value_name="Mean_RelAbund_%")
                df_time = df_time.sort_values([day_col, "Genus"])

                colors = (PALETTE * ((len(taxa_top)//len(PALETTE))+1))[:len(taxa_top)]
                fig_tc = px.line(
                    df_time, x=day_col, y="Mean_RelAbund_%", color="Genus",
                    markers=True,
                    color_discrete_sequence=colors,
                    title="Genus time-course (Days {})".format(", ".join(map(str, sorted(df_time[day_col].unique()))))
                )
                fig_tc.update_layout(
                    xaxis_title="Day",
                    yaxis_title="Mean relative abundance (%)",
                    legend_title_text="Genus"
                )
                with left_tc:
                    st.plotly_chart(fig_tc, use_container_width=True)

        # stash for other tabs if needed
        st.session_state["_rel_abund"] = rel_abund
        st.session_state["_micro_cols"] = micro_cols

# ====== INTERPRETATIONS ======
with tab_interp:
    st.markdown("### Model drivers & educational notes")

    # ---- 1) Top features (odds ratios) from artifacts
    top_or_csv = next(iter(sorted(artifacts_dir.glob("top*_coeffs_with_OR.csv"))), None)
    if top_or_csv and top_or_csv.exists():
        try:
            df_or = pd.read_csv(top_or_csv)

            # compute OR if not present
            if "odds_ratio_per_SD" not in df_or.columns and "Coefficient" in df_or.columns:
                df_or["odds_ratio_per_SD"] = np.exp(df_or["Coefficient"])

            # filter out day-related & unnamed features
            def _bad_feature(name: str) -> bool:
                n = str(name).lower()
                return n.startswith("day") or " day" in n or "unnamed:" in n
            df_or = df_or[~df_or["Feature"].astype(str).apply(_bad_feature)].copy()

            df_or["abs_coef"] = df_or["Coefficient"].abs()
            showN = st.slider("How many features to show", 5, 30, 20, 1, key="feat_n")
            show = df_or.sort_values("abs_coef", ascending=False).head(showN)

            # Positive vs Negative coloring
            show_plot = show.sort_values("odds_ratio_per_SD")
            show_plot["Direction"] = np.where(show_plot["Coefficient"] >= 0, "Positive", "Negative")

            fig3 = px.bar(
                show_plot,
                x="odds_ratio_per_SD", y="Feature", orientation="h",
                color="Direction",
                color_discrete_map={"Positive": "#a8c69f", "Negative": "#b91c1c"},
                title="Top model features by effect size (Odds Ratio per 1 SD)"
            )
            fig3.update_layout(xaxis_title="Odds Ratio (per 1 SD)", yaxis_title="")
            st.plotly_chart(fig3, use_container_width=True)

            st.dataframe(
                show[["Feature","Coefficient","odds_ratio_per_SD"]]
                .rename(columns={"odds_ratio_per_SD":"OR_per_SD"}),
                use_container_width=True, height=340
            )

            st.caption("Interpretation: OR > 1 increases odds of 'Unsafe'; OR < 1 decreases odds. Coefficients are pre-logit.")
        except Exception as e:
            st.info(f"Couldn't render feature odds ratios: {e}")
    else:
        st.info("No 'top*_coeffs_with_OR.csv' found in artifacts to visualize top features.")

    st.markdown("---")

    # ---- 2) Top taxa: late/high-risk vs early/low-risk (from relative abundance)
    st.markdown("#### Top taxa observed in spoilage")
    rel_abund = st.session_state.get("_rel_abund", None)
    micro_cols = st.session_state.get("_micro_cols", None)

    if rel_abund is None or not micro_cols:
        st.info("No numeric microbial columns detected; cannot propose top taxa lists.")
    else:
        threshold_for_lists = max(prob_thr, 0.50)
        high_mask = prob_unsafe >= threshold_for_lists
        low_mask  = prob_unsafe < threshold_for_lists

        top_high = rel_abund[high_mask].mean().sort_values(ascending=False).head(10) * 100.0 if high_mask.sum() > 0 else pd.Series(dtype=float)
        top_low  = rel_abund[low_mask].mean().sort_values(ascending=False).head(10) * 100.0 if low_mask.sum() > 0 else pd.Series(dtype=float)

        cA, cB = st.columns(2)
        with cA:
            st.markdown("**Top 10 taxa — Late / High-risk group**")
            if not top_high.empty:
                st.dataframe(
                    top_high.rename("Mean_RelAbund_%").reset_index().rename(columns={"index":"Taxon"}),
                    use_container_width=True, height=340
                )
            else:
                st.write("—")
        with cB:
            st.markdown("**Top 10 taxa — Early / Low-risk group**")
            if not top_low.empty:
                st.dataframe(
                    top_low.rename("Mean_RelAbund_%").reset_index().rename(columns={"index":"Taxon"}),
                    use_container_width=True, height=340
                )
            else:
                st.write("—")

        st.markdown("##### Interpretations (markdown supported)")
        default_text = "Add notes about early/late taxa here…"
        st.text_area("Notes:", value=default_text, height=220, key="notes_md")

# ====== REFERENCES (with download) ======
with tab_refs:
    st.markdown("### Key References")

    default_refs_md = """\
- **ICMSF** (2018). *Microorganisms in Foods 8: Use of Data for Assessing Process Control and Product Acceptance.*
- **Kraken2 / Bracken** methodology papers for metagenomic taxonomic profiling.
- **Pseudomonas & Brochothrix in chilled meats** — foundational studies on early vs late spoilage roles. *(add specific citations/DOIs)*
- **Predictive microbiology / logistic regression** — standard OR interpretation resources. *(add sources)*
- **Dataset sources** (BioProjects, DOIs, institutional datasets). *(list exact accessions)*
"""
    refs_md = st.text_area("References (Markdown):", value=default_refs_md, height=260, key="refs_md")

    include_notes = st.checkbox("Append 'Interpretations' notes to export", value=True, key="refs_append")
    export_buf = BytesIO()
    export_text = "# References\n\n" + refs_md.strip() + "\n"
    if include_notes and "notes_md" in st.session_state:
        export_text += "\n\n---\n\n# Interpretations Notes\n\n" + st.session_state["notes_md"].strip() + "\n"
    export_buf.write(export_text.encode("utf-8"))
    st.download_button(
        "⬇️ Download references (.md)",
        data=export_buf.getvalue(),
        file_name="references_and_notes.md",
        mime="text/markdown"
    )
