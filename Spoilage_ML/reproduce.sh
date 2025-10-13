# 0) (optional) create a clean env
conda create -n spoilage_env python=3.11 -y
conda activate spoilage_env

# 1) install deps
pip install -r requirements.txt

# 2) run baseline
./spoilage_ml.py --csv fully_combined_project_data.csv --outdir results

# 3) class-balanced + threshold sweep/ROC + top20 OR plots (recommended headline)
./spoilage_ml.py --csv fully_combined_project_data.csv --outdir results_bal_thr40 --balanced --threshold 0.40 --topn 20

# 4) microbe-only variant (optional)
./spoilage_ml.py --csv fully_combined_project_data.csv --outdir results_microbe --microbe_only --balanced --threshold 0.40 --topn 20
