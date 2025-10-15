
import sys
from pathlib import Path

# --- Avoid import-name conflict with local 'dash' folder ---
REPO_ROOT = Path(__file__).parent.resolve()
if str(REPO_ROOT) in sys.path:
    sys.path.remove(str(REPO_ROOT))

from dash import Dash, html, dcc, Input, Output
import plotly.express as px
import pandas as pd
import numpy as np

# ---------- Relative data paths (portable) ----------
BASE_DIR = Path(__file__).parent
DATA_DIR = BASE_DIR / "dash" / "data"
META_PATH = DATA_DIR / "metadata.tsv"
BRACKEN_PATH = DATA_DIR / "bracken_genus_matrix.tsv"

# ---------- Utilities ----------
def _read_tsv_any(path: Path) -> pd.DataFrame:
    """Robust TSV reader that tries common encodings."""
    last_err = None
    for enc in ["utf-8", "utf-8-sig", "latin1", "cp1252", "utf-16", "utf-16le", "utf-16be"]:
        try:
            return pd.read_csv(path, sep="\t", encoding=enc)
        except Exception as e:
            last_err = e
    raise last_err

def _optional_load(path: Path):
    if not path.exists():
        return None
    try:
        return _read_tsv_any(path)
    except Exception:
        return None

# ---------- Load data (if present) ----------
meta = _optional_load(META_PATH)
bracken = _optional_load(BRACKEN_PATH)

# Clean metadata
PRODUCT_OPTIONS = []
sampling_time_values = []
if meta is not None:
    meta = meta.drop(columns=[c for c in meta.columns if str(c).startswith("Unnamed")], errors="ignore")
    if "Meat_Type" in meta.columns:
        meta["Meat_Type"] = meta["Meat_Type"].astype(str).str.strip().str.title()
    else:
        meta["Meat_Type"] = np.nan
    PRODUCT_OPTIONS = sorted([x for x in meta["Meat_Type"].dropna().unique().tolist() if x])
    if "Sampling_Time" in meta.columns:
        sampling_time_values = sorted(meta["Sampling_Time"].dropna().astype(str).str.strip().unique().tolist())

# Clean Bracken
GENUS_COL = "name"
BRACKEN_SAMPLES = []
if bracken is not None and GENUS_COL in bracken.columns:
    bracken.columns = [str(c) for c in bracken.columns]  # ensure strings
    BRACKEN_SAMPLES = [c for c in bracken.columns if c != GENUS_COL]

# Defaults if data missing
if not PRODUCT_OPTIONS:
    PRODUCT_OPTIONS = ["Pork", "Poultry"]
if not sampling_time_values:
    sampling_time_values = ["T0", "T1", "T2", "T3"]

# ---------- App constants ----------
DEFAULT_TIME_TO_DAYS = {"T0": 0, "T1": 3, "T2": 7, "T3": 10}
SPOILAGE_THRESHOLDS = {"Pork": 5, "Poultry": 3}

def predict_spoilage(product: str, days: int):
    product_norm = (product or "").strip().title()
    thr = SPOILAGE_THRESHOLDS.get(product_norm, 4)
    if days is None:
        return "Unknown", {}
    status = "Spoiled" if int(days) > thr else "Fresh"
    return status, {}

def composition_figure_from_sample(bracken_df: pd.DataFrame, sample_id: str, top_n: int = 8):
    if bracken_df is None or GENUS_COL not in bracken_df.columns:
        return None
    if not sample_id or sample_id not in bracken_df.columns:
        return None
    ser = bracken_df[[GENUS_COL, sample_id]].copy()
    ser.columns = [GENUS_COL, "Abundance"]
    ser = ser.replace([np.inf, -np.inf], np.nan).dropna(subset=["Abundance"])
    ser = ser[ser["Abundance"] > 0]
    if ser.empty:
        return None
    ser_sorted = ser.sort_values("Abundance", ascending=False)
    top = ser_sorted.head(top_n)
    other_sum = ser_sorted.iloc[top_n:]["Abundance"].sum()
    if other_sum > 0:
        top = pd.concat([top, pd.DataFrame({GENUS_COL: ["Other"], "Abundance": [other_sum]})], ignore_index=True)
    fig = px.pie(top, names=GENUS_COL, values="Abundance", title=f"Microbial Composition — {sample_id}")
    return fig

# ---------- Dash app ----------
app = Dash(
    __name__,
    meta_tags=[{"name": "viewport", "content": "width=device-width, initial-scale=1"}],
)

app.layout = html.Div([
    html.H1("Food Spoilage Predictor"),
    html.Div("Place your data files in ./dash/data: metadata.tsv and bracken_genus_matrix.tsv."),

    html.Label("Product Type:"),
    dcc.Dropdown(
        id="product-dropdown",
        options=[{"label": p, "value": p} for p in PRODUCT_OPTIONS],
        placeholder="Select a product type"
    ),

    html.Label("Storage Duration (days):"),
    dcc.Slider(id="days-slider", min=0, max=30, step=1, value=0,
               marks={0:"0", 3:"3", 5:"5", 7:"7", 10:"10", 14:"14", 21:"21", 28:"28"}),
    html.Div(id="days-readout", style={"marginBottom": 12}),

    html.Label("Sampling Time (optional):"),
    dcc.Dropdown(
        id="sampling-time-dropdown",
        options=[{"label": s, "value": s} for s in sampling_time_values],
        placeholder="T0/T1/T2/T3"
    ),

    html.Label("Sample ID (from Bracken):"),
    dcc.Dropdown(
        id="sample-dropdown",
        options=[{"label": s, "value": s} for s in BRACKEN_SAMPLES],
        placeholder="Select SRR/ERR sample…",
        searchable=True, clearable=True, style={"maxWidth": 520}
    ),
    html.Div(f"Loaded {len(BRACKEN_SAMPLES)} samples from dash/data/bracken_genus_matrix.tsv" if BRACKEN_SAMPLES else
             "No Bracken data found. Place bracken_genus_matrix.tsv in ./dash/data/", style={"color": "#666", "marginTop": 4}),

    html.Hr(),
    html.H3("Prediction Results:"),
    html.Div(id="spoilage-status", style={"fontWeight": "bold"}),
    dcc.Graph(id="microbial-composition-graph", style={"display": "none"}),
    html.Div(id="helper-message", style={"color": "#666", "marginTop": 8})
])

@app.callback(Output("days-readout", "children"), Input("days-slider", "value"))
def _show_days(v):
    return f"Days stored: {int(v)}"

@app.callback(Output("days-slider", "value"), Input("sampling-time-dropdown", "value"), prevent_initial_call=True)
def fill_days_from_sampling_time(sampling_code):
    if sampling_code:
        return DEFAULT_TIME_TO_DAYS.get(str(sampling_code).strip(), 0)
    return dash.no_update

@app.callback(
    [Output("spoilage-status", "children"),
     Output("microbial-composition-graph", "figure"),
     Output("microbial-composition-graph", "style"),
     Output("helper-message", "children")],
    [Input("product-dropdown", "value"),
     Input("days-slider", "value"),
     Input("sample-dropdown", "value")]
)
def update_prediction(product, days, sample_id):
    # Spoilage status
    if product is not None and days is not None:
        status, _ = predict_spoilage(product, int(days))
        status_text = f"The product is predicted to be: {status} (Product: {product}, Days: {int(days)})."
    elif product is not None:
        status_text = "Move the slider or pick a Sampling Time to set days."
    else:
        status_text = "Select a product to get a spoilage prediction."

    # Composition figure
    helper = ""
    fig = None
    style = {"display": "none"}
    if sample_id:
        fig = composition_figure_from_sample(bracken, sample_id, top_n=8)
        if fig is None:
            helper = ("Sample found, but abundances are zero/missing or Bracken file is absent. "
                      "Ensure bracken_genus_matrix.tsv is in ./dash/data and the sample has nonzero values.")
        else:
            style = {"display": "block"}
    return status_text, (fig or {}), style, helper

if __name__ == "__main__":
    app.run(debug=False)
