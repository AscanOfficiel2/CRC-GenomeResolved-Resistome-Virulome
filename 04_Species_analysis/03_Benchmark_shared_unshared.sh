#!/bin/bash
#SBATCH --job-name=benchmark_shared_unshared
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/benchmark_shared_unshared_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/benchmark_shared_unshared_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 << 'EOF'
# =====================================================================
# BENCHMARKING FOR SHARED (CANONICAL) & UNSHARED (NONCANONICAL)
# - Model-specific colors
# - BLACK borders on 3×3 benchmark plots
# - All output files labeled with _shared / _unshared
# - Correct .str.strip() to avoid Series errors
# =====================================================================

import pandas as pd
import numpy as np
import os

from sklearn.model_selection import StratifiedKFold, cross_validate, train_test_split
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    cohen_kappa_score, matthews_corrcoef, make_scorer
)
from sklearn.preprocessing import LabelBinarizer

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

import seaborn as sns
import matplotlib.pyplot as plt

# =====================================================================
# PATHS
# =====================================================================

DATA_PATH = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/ML/data"

SHARED_FILE   = f"{DATA_PATH}/CLR_Shared_CRC_Markers_Only.csv"
UNSHARED_FILE = f"{DATA_PATH}/CLR_Unshared_Species_Only.csv"
META_FILE     = f"{DATA_PATH}/CRC_metadata_aligned_bact.csv"


# =====================================================================
# LOAD METADATA (FIXED LINE)
# =====================================================================

meta = pd.read_csv(META_FILE)
meta["Sample_ID"] = meta["Sample_ID"].astype(str).str.strip()   # <-- FIXED
meta = meta.set_index("Sample_ID")
y = meta["Group"].values


# =====================================================================
# LOAD MATRICES
# =====================================================================

shared_df   = pd.read_csv(SHARED_FILE, index_col=0).loc[meta.index]
unshared_df = pd.read_csv(UNSHARED_FILE, index_col=0).loc[meta.index]


# =====================================================================
# MACRO-MCC
# =====================================================================

def macro_mcc(y_true, y_pred):
    lb = LabelBinarizer().fit(y_true)
    Y_true = lb.transform(y_true)
    Y_pred = lb.transform(y_pred)
    return np.mean([
        matthews_corrcoef(Y_true[:, i], Y_pred[:, i])
        for i in range(Y_true.shape[1])
    ])


# =====================================================================
# SCORERS
# =====================================================================

scorers = {
    "accuracy": make_scorer(accuracy_score),
    "balanced_accuracy": make_scorer(balanced_accuracy_score),
    "f1_macro": make_scorer(f1_score, average="macro"),
    "kappa": make_scorer(cohen_kappa_score),
    "mcc": make_scorer(matthews_corrcoef),
    "macro_mcc": make_scorer(macro_mcc)
}


# =====================================================================
# MODELS
# =====================================================================

models = {
    "DecisionTree": DecisionTreeClassifier(random_state=42, class_weight="balanced"),
    "RandomForest": RandomForestClassifier(random_state=42, n_jobs=-1, class_weight="balanced"),
    "LogisticRegression": LogisticRegression(max_iter=5000, random_state=42, class_weight="balanced"),
    "SVM": SVC(kernel="rbf", probability=True, random_state=42, class_weight="balanced"),
    "GradientBoosting": GradientBoostingClassifier(random_state=42)
}


# =====================================================================
# MODEL COLORS
# =====================================================================

MODEL_COLORS = {
    "DecisionTree": "#1f77b4",       # blue
    "RandomForest": "#ff7f0e",       # orange
    "LogisticRegression": "#2ca02c", # green
    "SVM": "#d62728",                # red
    "GradientBoosting": "#9467bd"    # purple
}


# =====================================================================
# BENCHMARKING FUNCTION
# =====================================================================

def run_benchmark(X_df, y, tag):

    label = tag.capitalize()
    outdir = f"Benchmark_{tag}"
    os.makedirs(outdir, exist_ok=True)

    print(f"\n=========== BENCHMARKING {label} DATASET ===========")

    X = X_df.values

    # Train/Test split
    X_train, _, y_train, _ = train_test_split(
        X, y, test_size=0.20, stratify=y, random_state=42
    )

    # 10-fold CV
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    results = []

    for name, model in models.items():
        cv_out = cross_validate(
            model, X_train, y_train,
            cv=cv, scoring=scorers,
            n_jobs=-1
        )
        for metric in scorers.keys():
            results.append({
                "Model": name,
                "Metric": metric,
                "Mean": np.mean(cv_out[f"test_{metric}"]),
                "SD": np.std(cv_out[f"test_{metric}"])
            })

    df_results = pd.DataFrame(results)
    df_results.to_csv(f"{outdir}/Benchmark_AllModels_{tag}.csv", index=False)


    # =====================================================================
    # 3×3 BENCHMARK FIGURE (MODEL COLORS + BLACK BORDERS)
    # =====================================================================

    df_plot = df_results.copy()
    df_plot["Metric"] = df_plot["Metric"].str.replace("_"," ").str.title()
    df_plot["Metric"] = df_plot["Metric"].replace({"Mcc":"MCC","Macro Mcc":"Macro-MCC"})

    metric_order = ["Accuracy","Balanced Accuracy","F1 Macro","Kappa","MCC","Macro-MCC"]
    model_order  = list(models.keys())

    fig, axes = plt.subplots(3,3, figsize=(16,14), dpi=350)
    axes = axes.flatten()

    for i, metric in enumerate(metric_order):
        ax = axes[i]

        sub = df_plot[df_plot["Metric"] == metric]
        sub = sub.set_index("Model").loc[model_order].reset_index()

        # model-colored points + error bars
        for _, row in sub.iterrows():
            m = row["Model"]
            mc = MODEL_COLORS[m]
            ax.errorbar(
                row["Mean"], m,
                xerr=row["SD"],
                fmt="o", color=mc,
                markersize=8, lw=1.8, capsize=4
            )

        ax.set_title(metric, fontsize=14, weight="bold")
        ax.set_xlim(0.2, 0.95)
        ax.grid(False)
        ax.set_ylabel("")

        # BLACK borders
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor("black")
            spine.set_linewidth(2.0)

    # hide empty plots
    for j in range(len(metric_order),9):
        axes[j].axis("off")

    plt.tight_layout()
    plt.savefig(f"{outdir}/Benchmark_3x3_{tag}.png", dpi=600)
    plt.close()


    # =====================================================================
    # RADAR PLOT
    # =====================================================================

    radar = df_results.pivot(index="Model", columns="Metric", values="Mean")
    radar = radar.loc[model_order, ["accuracy","balanced_accuracy","f1_macro","kappa","mcc","macro_mcc"]]

    radar_norm = (radar - radar.min())/(radar.max()-radar.min())

    labels = radar_norm.columns
    angles = np.linspace(0,2*np.pi,len(labels),endpoint=False)
    angles = np.concatenate((angles,[angles[0]]))

    plt.figure(figsize=(6,5), dpi=600)

    for model in radar_norm.index:
        vals = radar_norm.loc[model].tolist() + [radar_norm.loc[model].tolist()[0]]
        plt.polar(
            angles, vals,
            marker="o", linewidth=1.5,
            label=model, color=MODEL_COLORS[model]
        )

    plt.xticks(angles[:-1], labels, fontsize=7)
    plt.yticks(fontsize=6)
    plt.legend(loc="upper right", bbox_to_anchor=(1.7,1.1), fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{outdir}/Benchmark_RadarPlot_{tag}.png", dpi=600)
    plt.close()


    print(f"✓ Finished Benchmark: {label}\n")


# =====================================================================
# RUN BENCHMARKS
# =====================================================================

run_benchmark(shared_df,   y, tag="shared")
run_benchmark(unshared_df, y, tag="unshared")

print("\n🎉 BENCHMARKING COMPLETE FOR SHARED & UNSHARED DATASETS 🎉")

EOF
