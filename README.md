# 🥩🍖 Microbe-Driven Spoilage Predictor

Streamlit app for exploring safe/unsafe spoilage predictions and microbiome visualizations.  
The app uses pre-trained models and saved artifacts — no retraining required.

---

## Files You Need
Clone or download the repo and make sure you have **all of these** in the same project folder (or point the app to them):

- `spoilage_app.py`  ← the Streamlit app  
- `results_with_day.csv`  ← your dataset with predictions (can also upload CSV manually in the app)  
- `fully_combined_project_data.csv`  
- `toy_spoilage.csv`  
- `spoilage_regr.py`  
- `clf_pipeline.joblib`  
- `coefficients_logit.csv`  
- `confusion_matrix.txt`  
- `logit_model.pkl`  
- `metrics.txt`  
- `model_meta.json`  
- `reg_pipeline.joblib`  
- `reg_r2.txt`  
- `roc_curve.png`  
- `threshold_sweep.csv`  
- `top20_coeffs_bar.png`  
- `top20_coeffs_with_OR.csv`  
- `top20_odds_ratios.png`  

These artifacts are required for the app to load models, display coefficients/odds ratios, and show reference plots.

---

## Installation

Make sure you’re in a clean virtual environment (conda or venv). Then:

```bash
pip install -r requirements.txt

streamlit run spoilage_app.py

```

In the sidebar, set the Artifacts folder to wherever the model files live (defaults to ~/Desktop/results_with_day).
Required inside the artifacts folder:
logit_model.pkl
reg_pipeline.joblib
model_meta.json
You can upload a CSV (results_with_day.csv or your own) with the same schema as the training data to generate predictions.
Adjust the “Unsafe probability threshold” slider to explore sensitivity.
Outputs include prediction tables, CFU regressions, risk histograms, microbiome relative abundances, and top model features.