import pandas as pd
import numpy as np
import hyperopt as hp
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_val_score
from joblib import dump, load
import json
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials, space_eval

# (I've hidden the 'display' import as it's not used)
# from IPython.display import display 

spoilage_data = pd.read_csv(r"fully_combined_project_data.csv")

# Create the mapping
type_map = {'pork': 0, 'poultry': 1}
spoilage_data['meat_type'] = spoilage_data['EnvType'].map(type_map)


# Drop columns that contain nan values
spoilage_data = spoilage_data.dropna(axis=1)

# Select Columns that are numeric
num_spoilage_data = spoilage_data.select_dtypes(include=np.number)

num_spoilage_data['earlyvlatespoilage'] = (num_spoilage_data['Total mesophilic aerobic flora (log10 CFU.g-1)'] >= 7).astype(int)

# print(num_spoilage_data[['Total mesophilic aerobic flora (log10 CFU.g-1)', 'earlyvlatespoilage']])

X = num_spoilage_data.drop(['earlyvlatespoilage', 'Total mesophilic aerobic flora (log10 CFU.g-1)'], axis=1)
Y = num_spoilage_data['earlyvlatespoilage']

# X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.25 , random_state=100, shuffle=True)


# --- START RANDOM FOREST HYPEROPT ---

# Define the objective function to minimize (1 - ROC_AUC)
def objective(params):
    # Cast float values from hyperopt to int where required
    params['n_estimators'] = int(params['n_estimators'])
    params['max_depth'] = int(params['max_depth'])
    params['min_samples_leaf'] = int(params['min_samples_leaf'])

    model = RandomForestClassifier(
        **params,
        random_state=0,
        n_jobs=-1  # Use all available cores
    )

    # Perform cross-validation on the full dataset
    score = cross_val_score(model, X, Y, cv=5, scoring='roc_auc', n_jobs=-1)

    # We want to MINIMIZE loss, so we return (1 - mean_score)
    loss = 1 - np.mean(score)

    return {'loss': loss, 'status': STATUS_OK}

# Define the search space for Random Forest hyperparameters
space = {
    'n_estimators': hp.quniform('n_estimators', 100, 500, 50),
    'max_depth': hp.quniform('max_depth', 10, 50, 5),
    'min_samples_leaf': hp.quniform('min_samples_leaf', 1, 10, 1),
    'max_features': hp.choice('max_features', ['sqrt', 'log2', 0.75])
}

# Run the optimization
trials = Trials()
best = fmin(
    fn=objective,
    space=space,
    algo=tpe.suggest,
    max_evals=50,
    trials=trials
)

print("Best hyperparameters found (raw):", best)
print("Minimum loss (1-ROC_AUC) achieved:", trials.best_trial['result']['loss'])

# --- TRAIN FINAL MODEL ---

# Get the actual parameter values from the 'best' dictionary
final_params = space_eval(space, best)

# Cast numerical params to int for the final model
final_params['n_estimators'] = int(final_params['n_estimators'])
final_params['max_depth'] = int(final_params['max_depth'])
final_params['min_samples_leaf'] = int(final_params['min_samples_leaf'])

print("Training final model with best params:", final_params)

# Train the final classifier on the FULL dataset
clf = RandomForestClassifier(
    **final_params,
    random_state=0,
    n_jobs=-1,
    oob_score=True
).fit(X, Y)

print(f"\nModel OOB (out-of-bag) score on training data: {clf.oob_score_:.3f}")

# --- SAVE ARTIFACTS ---

# 1. Save Feature Importances
importances_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance': clf.feature_importances_  # <<< Use correct column name
})

print("\nTop 10 Features (by Importance):")
print(importances_df.sort_values(by='Importance', ascending=False).head(10))

# Save with a new, specific name
importances_df.to_csv("rf_feature_importances.csv", index=False) # <<< CHANGED
print("\nSuccessfully saved feature importances to rf_feature_importances.csv")

# 2. Save the Model
# Save with a new, specific name
filename = "rf_model_tuned.joblib" # <<< CHANGED
dump(clf, filename)
print(f"Successfully saved model to {filename}")

# 3. Save the meta.json file (this name is generic and can stay the same)
print("Creating model_meta.json...")
meta_data = {
    "feature_names": list(X.columns), 
    "threshold_cfu": 7.0 
}

with open("model_meta.json", "w") as f:
    json.dump(meta_data, f, indent=2)

print("Successfully created model_meta.json")