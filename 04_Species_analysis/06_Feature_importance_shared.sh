#!/bin/bash
#SBATCH --job-name=feature_importance_shared
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=36:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/feature_importance_shared_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/feature_importance_shared_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate


python3 << 'EOF'
# ============================================================
# FEATURE IMPORTANCE STABILITY PIPELINE — SHARED DATASET
# Models: GradientBoosting + RandomForest
# GLOBAL + PER-CLASS (Healthy, Adenoma, Cancer)
# CV=5, Repeats=3, Bootstrap=100
# PLOTS: Top 20 only, thick black borders
# COLORS:
#   - Global plot → viridis (as full dataset)
#   - Per-class → magma (as full dataset)
# ============================================================

import pandas as pd
import numpy as np
import os, joblib
from tqdm import tqdm
from sklearn.model_selection import StratifiedKFold
import seaborn as sns
import matplotlib.pyplot as plt


# ------------------ SETTINGS ------------------
N_FOLDS    = 5
N_REPEATS  = 3
N_BOOT     = 100
RAND       = 42
TOPN       = 20

# ------------------ PATHS ------------------
DATA_PATH = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/ML/data"

SHARED_FILE = f"{DATA_PATH}/CLR_Shared_CRC_Markers_Only.csv"
META_FILE   = f"{DATA_PATH}/CRC_metadata_aligned_bact.csv"

GB_SHARED = f"{DATA_PATH}/NestedCV_shared/GB_best_shared.pkl"
RF_SHARED = f"{DATA_PATH}/NestedCV_shared/RF_best_shared.pkl"

OUT     = "FeatureImportance_shared"
OUT_CLS = "FeatureImportance_shared/classes"

os.makedirs(OUT, exist_ok=True)
os.makedirs(OUT_CLS, exist_ok=True)


# ------------------ LOAD DATA ------------------
meta = pd.read_csv(META_FILE)
meta["Sample_ID"] = meta["Sample_ID"].astype(str).str.strip()
meta = meta.set_index("Sample_ID")

y = meta["Group"].astype(str).str.strip().values
CLASSES = ["Healthy", "Adenoma", "Cancer"]

X_df = pd.read_csv(SHARED_FILE, index_col=0).loc[meta.index]
X = X_df.values
species = X_df.columns.tolist()


# ------------------ LOAD MODELS ------------------
gb_model = joblib.load(GB_SHARED)
rf_model = joblib.load(RF_SHARED)

models = {
    "GB": gb_model,
    "RF": rf_model
}


# ------------------ IMPORTANCE FUNCTIONS ------------------
def gb_imp(model): return model.feature_importances_
def rf_imp(model): return model.feature_importances_


# ------------------ CONSENSUS TABLE ------------------
def consensus_table(imp_dict, species):
    df = pd.DataFrame({"Species": species})
    for name, vals in imp_dict.items():
        df[name] = np.mean(vals, axis=0)
        df[name + "_norm"] = df[name] / df[name].max()
    norm_cols = [c for c in df.columns if c.endswith("_norm")]
    df["Consensus"] = df[norm_cols].mean(axis=1)
    return df.sort_values("Consensus", ascending=False)


# ------------------ PLOTTING FUNCTION ------------------
def plot_top(df, title, outname, palette="magma"):
    plt.figure(figsize=(8,10))
    sns.barplot(
        data=df.head(TOPN),
        x="Consensus",
        y="Species",
        palette=palette,
        edgecolor="black",
        linewidth=2.2
    )
    plt.title(title, fontsize=15, weight="bold")
    plt.tight_layout()
    plt.savefig(outname, dpi=600)
    plt.close()


# ============================================================
# 1. GLOBAL FEATURE IMPORTANCE
# ============================================================

print("\n=== GLOBAL FEATURE IMPORTANCE (SHARED) ===")

cv = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=RAND)
global_imp = {k: [] for k in models}
boot_imp   = {k: [] for k in models}

# ----- CV × Repeats -----
for rep in range(N_REPEATS):
    for tr, te in cv.split(X, y):
        X_tr, y_tr = X[tr], y[tr]
        for name, model in models.items():
            model.fit(X_tr, y_tr)
            imp = gb_imp(model) if name=="GB" else rf_imp(model)
            global_imp[name].append(imp)

# ----- Bootstraps -----
for name, model in models.items():
    print(f"Bootstrapping (Global): {name}")
    for b in tqdm(range(N_BOOT)):
        idx = np.random.choice(len(y), len(y), replace=True)
        X_b, y_b = X[idx], y[idx]
        model.fit(X_b, y_b)
        imp = gb_imp(model) if name=="GB" else rf_imp(model)
        boot_imp[name].append(imp)

# ----- Consensus -----
global_cons = consensus_table(
    {f"{k}_cv":v for k,v in global_imp.items()} |
    {f"{k}_boot":v for k,v in boot_imp.items()},
    species
)

global_cons.to_csv(f"{OUT}/Global_Consensus.csv", index=False)

# GLOBAL COLOR = VIRIDIS
plot_top(global_cons,
         "Top 20 Global Consensus Features (Shared)",
         f"{OUT}/Global_Top20.png",
         palette="viridis")



# ============================================================
# 2. PER-CLASS FEATURE IMPORTANCE
# ============================================================

print("\n=== PER-CLASS FEATURE IMPORTANCE (SHARED) ===")

for cls in CLASSES:

    print(f"\n→ Processing class: {cls}")

    y_bin = np.where(y == cls, 1, 0)

    if y_bin.sum() < 3:
        print(f"Not enough positives for class {cls}, skipping…")
        continue

    cls_imp  = {k: [] for k in models}
    cls_boot = {k: [] for k in models}

    # ----- CV Repeats -----
    for rep in range(N_REPEATS):
        for tr, te in cv.split(X, y_bin):
            X_tr, y_tr = X[tr], y_bin[tr]
            for name, model in models.items():
                model.fit(X_tr, y_tr)
                imp = gb_imp(model) if name=="GB" else rf_imp(model)
                cls_imp[name].append(imp)

    # ----- Bootstraps -----
    for name, model in models.items():
        print(f"Bootstrapping {cls}: {name}")
        for b in tqdm(range(N_BOOT)):
            idx = np.random.choice(len(y), len(y), replace=True)
            X_b, y_b = X[idx], y_bin[idx]
            model.fit(X_b, y_b)
            imp = gb_imp(model) if name=="GB" else rf_imp(model)
            cls_boot[name].append(imp)

    # ----- Consensus -----
    cls_cons = consensus_table(
        {f"{k}_cv":v for k,v in cls_imp.items()} |
        {f"{k}_boot":v for k,v in cls_boot.items()},
        species
    )

    cls_cons.to_csv(f"{OUT_CLS}/{cls}_Consensus.csv", index=False)

    # PER-CLASS COLOR = MAGMA
    plot_top(
        cls_cons,
        f"Top 20 Consensus Features — {cls} (Shared)",
        f"{OUT_CLS}/{cls}_Top20.png",
        palette="magma"
    )

print("\n SHARED FEATURE IMPORTANCE COMPLETED SUCCESSFULLY")

EOF
