#!/usr/bin/env python3
import json, pickle
from pathlib import Path
from io import BytesIO

import joblib
from sklearn.metrics import confusion_matrix, classification_report # For Performance tab
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
st.caption("Classification is driven by your Random Forest model.") # <<< CHANGED

# ---------------- SIDEBAR SETTINGS ----------------
with st.sidebar:
    st.header("⚙️ Settings")
    # default_artifacts = str(Path.home() / "Desktop" / "results_with_day")
    default_artifacts = str(Path.home() / "OneDrive" / "Documents" / "UMGC" / "BIOT 670i" / "BIOT670i---Project-5" / "Streamlit" / "data")
    
    artifacts_dir_str = st.text_input(
        "Artifacts folder", value=default_artifacts,
        # <<< CHANGED help text for new RF files >>>
        help="Needs: rf_model_tuned.joblib, model_meta.json, rf_feature_importances.csv"
    )
    
    st.markdown(
        """
        ---
        **Disclaimer (RUO)** This tool is for **Research Use Only** — not for clinical/regulatory decisions.  
        Models were trained on **poultry & pork** datasets and may not generalize to all products.  
        A classification threshold of **0.50** is used.
        """,
        unsafe_allow_html=True
    )

artifacts_dir = Path(artifacts_dir_str).expanduser()

prob_thr = 0.50 # Hard-coded threshold
GT_THRESHOLD = 7.0  # Ground truth threshold for "not-safe"

# <<< CHANGED required file list for new RF files >>>
need = ["rf_model_tuned.joblib", "model_meta.json", "rf_feature_importances.csv"]
if not artifacts_dir.exists() or any(not (artifacts_dir / f).exists() for f in need):
    st.warning(f"Pick a valid artifacts folder containing: {need}")
    st.stop()

@st.cache_resource
def load_artifacts(folder: Path):
    # <<< CHANGED model filename >>>
    reg = joblib.load(folder / "rf_model_tuned.joblib")  # Your RF model
    meta = json.loads((folder / "model_meta.json").read_text())
    return reg, meta

try:
    reg, meta = load_artifacts(artifacts_dir)
except Exception as e:
    st.error(f"Failed to load artifacts: {e}")
    st.stop()

feature_names = meta.get("feature_names", [])
cfu_thr = float(meta.get("threshold_cfu", 7.0))
LOG_CFU_COL = "Total mesophilic aerobic flora (log10 CFU.g-1)" # With period


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


# 1) Get predictions from your Random Forest model
regressor_ok = True
pred_score = None 
try:
    # <<< CHANGED: Use .predict_proba() for classifiers >>>
    # This gets the probability for class 1 ("not-safe")
    pred_score = reg.predict_proba(X)[:, 1]
except Exception as e:
    st.error(f"Failed to get predictions from your model: {e}")
    regressor_ok = False
    st.stop()

# 2) Create "safe"/"not-safe" classification from that score
pred_class = np.where(pred_score >= prob_thr, "not-safe", "safe")

# 3) Calculate Confidence
safe_conf = (prob_thr - pred_score) / prob_thr
notsafe_conf = (pred_score - prob_thr) / (1.0 - prob_thr)
confidence_score = np.where(pred_class == "safe", safe_conf, notsafe_conf)
confidence_score = np.nan_to_num(confidence_score, nan=0.0)
confidence_score = np.clip(confidence_score, 0.0, 1.0)

# 4) Column mapping for display
ITEM_COL = "Sample_Name_y"
PTYPE_COL = "EnvType"
DAY_COL = "Day_numeric" if "Day_numeric" in df_raw.columns else ("Days_Numeric" if "Days_Numeric" in df_raw.columns else None)
if not DAY_COL:
    st.error("Neither 'Day_numeric' nor 'Days_Numeric' found in the CSV.")
    st.stop()

# Pastel palette for Plotly charts
PALETTE = ["#d5b4ab", "#c3b1e1", "#a3c1da", "#a8c69f", "#f4c2c2"]

def style_predictions(row):
    """Applies color to Prediction and Confidence cols based on the value."""
    styles = [''] * len(row) # default
    pred_val = row['Prediction']
    style_str = ''
    
    if pred_val == 'not-safe':
        style_str = 'background-color: #fde8e8; color: #9b1c1c;'
    elif pred_val == 'safe':
        style_str = 'background-color: #e8f5e9; color: #1b5e20;'

    try:
        pred_idx = row.index.get_loc('Prediction')
        conf_idx = row.index.get_loc('Confidence')
        styles[pred_idx] = style_str
        styles[conf_idx] = style_str
    except KeyError:
        pass 
    return styles

# ---------------- TABS ----------------
tab_pred, tab_viz, tab_perf, tab_micro, tab_interp, tab_refs = st.tabs(
    ["🔮 Predictions", "📈 Visuals", "📊 Performance", "🧬 Microbiome", "📚 Interpretations", "📖 References"]
)

# ====== PREDICTIONS ======
with tab_pred:
    st.markdown("### 📊 Predictions from your Model")

    disp = pd.DataFrame({
        "Days in Refrigerator": df_raw[DAY_COL],
        "Item": df_raw[ITEM_COL] if ITEM_COL in df_raw.columns else "",
        "Product Type": df_raw[PTYPE_COL] if PTYPE_COL in df_raw.columns else "",
        "Prediction": pred_class,
        "Confidence": confidence_score,
    })
    
    formatters = {"Confidence": "{:.1%}"}
    
    if LOG_CFU_COL in df_raw.columns:
        disp["Log CFU (Input)"] = df_raw[LOG_CFU_COL]
        formatters["Log CFU (Input)"] = "{:.2f}" 

    st.dataframe(
        disp.style
        .apply(style_predictions, axis=1)
        .format(formatters),
        use_container_width=True
    )

    st.download_button(
        "⬇️ Download predictions as CSV",
        disp.to_csv(index=False).encode("utf-8"),
        "spoilage_predictions.csv",
        "text/csv"
    )

# ====== VISUALS ======
with tab_viz:
    st.markdown("### Prediction Probability Distribution") # <<< CHANGED
    
    if regressor_ok and pred_score is not None:
        fig1, ax1 = plt.subplots(figsize=(7, 4))
        # Title and label now reflect "Probability"
        ax1.hist(pred_score, bins=25, edgecolor="white")
        ax1.set_xlabel("Predicted Probability of 'not-safe' (from your RF model)") # <<< CHANGED
        ax1.set_ylabel("Count")
        ax1.set_title("Distribution of Model Prediction Probabilities") # <<< CHANGED
        
        ax1.axvline(prob_thr, color='#e75480', linestyle='--', linewidth=2)
        ax1.text(prob_thr + 0.02, ax1.get_ylim()[1]*0.9, f'Threshold: {prob_thr:.2f}', color='#e75480')
        
        st.pyplot(fig1)
    else:
        st.info("Model scores not available to visualize.")

# ====== PERFORMANCE ======
with tab_perf:
    st.markdown("### 📊 Model Performance Metrics")
    st.caption(f"Compares model predictions (at {prob_thr} probability threshold) to ground truth (at {GT_THRESHOLD} log CFU threshold).")

    if LOG_CFU_COL not in df_raw.columns:
        st.warning(f"Cannot calculate performance. Input CSV must contain the ground truth column: **'{LOG_CFU_COL}'**")
    else:
        # Create true labels from the ground truth CFU column
        y_true = np.where(df_raw[LOG_CFU_COL] >= GT_THRESHOLD, "not-safe", "safe")
        # Predicted labels (already computed)
        y_pred = pred_class
        labels = ["safe", "not-safe"]

        st.markdown("---")
        st.markdown("#### Classification Report")
        try:
            report_dict = classification_report(y_true, y_pred, labels=labels, output_dict=True, zero_division=0)
            report_df = pd.DataFrame(report_dict).transpose()
            
            st.dataframe(
                report_df.style.format({
                    "precision": "{:.2f}",
                    "recall": "{:.2f}",
                    "f1-score": "{:.2f}",
                    "support": "{:.0f}"
                }),
                use_container_width=True
            )
        except Exception as e:
            st.error(f"Could not generate classification report: {e}")


        st.markdown("---")
        st.markdown("#### Confusion Matrix")
        try:
            cm = confusion_matrix(y_true, y_pred, labels=labels)
            
            fig_cm = px.imshow(cm, text_auto=True,
                               labels=dict(x="Predicted Label", y="True Label"),
                               x=labels, y=labels,
                               color_continuous_scale='Reds'
                              )
            fig_cm.update_layout(title="Confusion Matrix", coloraxis_showscale=False)
            st.plotly_chart(fig_cm, use_container_width=True)

        except Exception as e:
            st.error(f"Could not generate confusion matrix: {e}")

# ====== MICROBIOME (PLOTLY) ======
with tab_micro:
    st.markdown("### Relative abundance views")

    NON_MICROBE_COLS = {
        "Sample_Name_y","EnvType","EnvTyp",
        "Day_numeric","Days_Numeric","Day","Days",
        "SRR_ID","SAMN_ID","SamplingTime",
        "Total mesophilic aerobic flora (log10 CFU.g-1)",
        # predictions / outputs
        "Safe/Unsafe","Probability_UnsafeToday","Days Until Spoiled","log CFU (pred today)",
        "Prediction", "Confidence"
    }

    micro_cols = [
        c for c in df_raw.columns
        if c not in NON_MICROBE_COLS and np.issubdtype(df_raw[c].dtype, np.number)
    ]

    if not micro_cols:
        st.info("No numeric microbial columns detected to compute relative abundances.")
    else:
        micro_df = df_raw[micro_cols].copy()
        row_totals = micro_df.sum(axis=1).replace(0, np.nan)
        rel_abund = micro_df.div(row_totals, axis=0).fillna(0.0)

        # --- Single-sample Top-N (Plotly)
        left, right = st.columns([2, 1], gap="large")
        with left:
            sample_labels = df_raw["Sample_Name_y"].astype(str) if "Sample_Name_y" in df_raw.columns else df_raw.index.astype(str)
            default_idx = 0
            selected_label = st.selectbox("Select a sample", options=list(sample_labels), index=default_idx, key="sample_select")

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

        # ----- ⏱️ Time-course line chart -----
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

        st.session_state["_rel_abund"] = rel_abund
        st.session_state["_micro_cols"] = micro_cols

# ====== INTERPRETATIONS ======
with tab_interp:
    st.markdown("### Model drivers & educational notes")

    # <<< THIS ENTIRE BLOCK IS CHANGED FOR RANDOM FOREST >>>
    
    # ---- 1) Top features (importances) from your model
    coeff_csv = artifacts_dir / "rf_feature_importances.csv" # <<< CHANGED
    if coeff_csv.exists():
        try:
            df_import = pd.read_csv(coeff_csv) # <<< CHANGED
            
            def _bad_feature(name: str) -> bool:
                n = str(name).lower()
                return n.startswith("day") or " day" in n or "unnamed:" in n
            df_import = df_import[~df_import["Feature"].astype(str).apply(_bad_feature)].copy() # <<< CHANGED

            df_import["abs_coef"] = df_import["Importance"].abs() # <<< CHANGED
            showN = st.slider("How many features to show", 5, 30, 20, 1, key="feat_n")
            show = df_import.sort_values("abs_coef", ascending=False).head(showN) # <<< CHANGED

            # Importances are always positive, so no "Direction" needed
            show_plot = show.sort_values("Importance") # <<< CHANGED

            fig3 = px.bar(
                show_plot,
                x="Importance", y="Feature", orientation="h", # <<< CHANGED
                color_discrete_sequence=["#e75480"], # Use a single theme color
                title="Top Model Features by Importance (from your Random Forest model)" # <<< CHANGED
            )
            fig3.update_layout(xaxis_title="Feature Importance", yaxis_title="") # <<< CHANGED
            st.plotly_chart(fig3, use_container_width=True)

            st.dataframe(
                show[["Feature","Importance"]], # <<< CHANGED
                use_container_width=True, height=340
            )

            st.caption("Interpretation: Higher importance means the feature was used more by the model to make decisions.") # <<< CHANGED
        except Exception as e:
            st.info(f"Couldn't render feature importances: {e}")
    else:
        st.info("No 'rf_feature_importances.csv' found in artifacts to visualize top features.") # <<< CHANGED
    # <<< END OF CHANGED BLOCK >>>

    st.markdown("---")

    # ---- 2) Top taxa: late/high-risk vs early/low-risk (from relative abundance)
    st.markdown("#### Top taxa observed in spoilage")
    rel_abund = st.session_state.get("_rel_abund", None)
    micro_cols = st.session_state.get("_micro_cols", None)

    if rel_abund is None or not micro_cols or pred_score is None:
        st.info("No numeric microbial columns or model scores available; cannot propose top taxa lists.")
    else:
        threshold_for_lists = max(prob_thr, 0.50)
        
        high_mask = pred_score >= threshold_for_lists
        low_mask  = pred_score < threshold_for_lists

        top_high = rel_abund[high_mask].mean().sort_values(ascending=False).head(10) * 100.0 if high_mask.sum() > 0 else pd.Series(dtype=float)
        top_low  = rel_abund[low_mask].mean().sort_values(ascending=False).head(10) * 100.0 if low_mask.sum() > 0 else pd.Series(dtype=float)

        cA, cB = st.columns(2)
        with cA:
            st.markdown(f"**Top 10 taxa — High-risk (Score >= {threshold_for_lists:.2f})**")
            if not top_high.empty:
                st.dataframe(
                    top_high.rename("Mean_RelAbund_%").reset_index().rename(columns={"index":"Taxon"}),
                    use_container_width=True, height=340
                )
            else:
                st.write("—")
        with cB:
            st.markdown(f"**Top 10 taxa — Low-risk (Score < {threshold_for_lists:.2f})**")
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
- **Random Forest** - (Breiman, 2001). *(add sources)*
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