#!/bin/bash
#SBATCH --job-name=nestedCV_top2_shared_unshared
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=06:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/nestedCV_top2_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/nestedCV_top2_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate


python3 << 'EOF'
# ============================================================
# NESTED CV + HOLDOUT EVALUATION FOR TOP TWO MODELS
# Shared → Gradient Boosting + Random Forest
# Unshared → Gradient Boosting + Logistic Regression
#
# ALL FILES SAVED WITH _shared / _unshared SUFFIX
# MODELS, PREDICTIONS, METRICS, ROC, CONFUSION MATRIX
# ============================================================

import pandas as pd
import numpy as np
import os
import joblib

import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, f1_score,
    cohen_kappa_score, matthews_corrcoef, confusion_matrix,
    roc_curve, auc
)
from sklearn.preprocessing import label_binarize, StandardScaler
from sklearn.pipeline import Pipeline

from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression


# ============================================================
# PATHS
# ============================================================

DATA_PATH = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/ML/data"

SHARED_FILE   = f"{DATA_PATH}/CLR_Shared_CRC_Markers_Only.csv"
UNSHARED_FILE = f"{DATA_PATH}/CLR_Unshared_Species_Only.csv"
META_FILE     = f"{DATA_PATH}/CRC_metadata_aligned_bact.csv"


# ============================================================
# LOAD METADATA
# ============================================================

meta = pd.read_csv(META_FILE)
meta["Sample_ID"] = meta["Sample_ID"].astype(str).str.strip()
meta = meta.set_index("Sample_ID")

y = meta["Group"].values
classes = np.unique(y)


# ============================================================
# LOAD MATRICES
# ============================================================

shared_df   = pd.read_csv(SHARED_FILE, index_col=0).loc[meta.index]
unshared_df = pd.read_csv(UNSHARED_FILE, index_col=0).loc[meta.index]


# ============================================================
# CORE FUNCTION — PERFORMS NESTED CV + HOLDOUT
# SAVES EVERYTHING WITH TAGS (_shared / _unshared)
# ============================================================

def nestedCV_model(X_df, y, tag, model_name, model, param_grid):

    print(f"\n========== Nested CV for {tag.upper()} — {model_name} ==========")

    outdir = f"NestedCV_{tag}"
    os.makedirs(outdir, exist_ok=True)

    X = X_df.values

    # ----------------------------
    # TRAIN / TEST SPLIT
    # ----------------------------
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.20, stratify=y, random_state=42
    )

    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    inner = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

    gs = GridSearchCV(
        model, param_grid,
        cv=inner, scoring="balanced_accuracy",
        n_jobs=-1
    )

    outer_scores = []

    # ----------------------------
    # OUTER LOOP
    # ----------------------------
    for fold, (tr_idx, te_idx) in enumerate(outer.split(X_train, y_train), 1):
        X_tr, X_te = X_train[tr_idx], X_train[te_idx]
        y_tr, y_te = y_train[tr_idx], y_train[te_idx]

        gs.fit(X_tr, y_tr)
        best_model = gs.best_estimator_

        y_pred = best_model.predict(X_te)
        bal_acc = balanced_accuracy_score(y_te, y_pred)
        outer_scores.append(bal_acc)

        print(f"Fold {fold} Balanced Accuracy = {bal_acc:.3f}")

    mean_outer = np.mean(outer_scores)
    print(f"\nMean Outer Balanced Accuracy = {mean_outer:.3f}")
    print("Best Params:", gs.best_params_)

    # Save best params
    pd.DataFrame([gs.best_params_]).to_csv(
        f"{outdir}/best_params_{model_name}_{tag}.csv",
        index=False
    )

    # ----------------------------
    # TRAIN FINAL BEST MODEL
    # ----------------------------
    best_model = gs.best_estimator_
    best_model.fit(X_train, y_train)

    # Save model
    joblib.dump(best_model, f"{outdir}/{model_name}_best_{tag}.pkl")
    print(f"Saved: {model_name}_best_{tag}.pkl")

    # ----------------------------
    # HOLDOUT EVALUATION
    # ----------------------------
    y_pred = best_model.predict(X_test)
    y_prob = best_model.predict_proba(X_test)

    metrics = {
        "Accuracy": accuracy_score(y_test, y_pred),
        "BalancedAcc": balanced_accuracy_score(y_test, y_pred),
        "F1_macro": f1_score(y_test, y_pred, average="macro"),
        "Kappa": cohen_kappa_score(y_test, y_pred),
        "MCC": matthews_corrcoef(y_test, y_pred)
    }

    pd.DataFrame([metrics]).to_csv(
        f"{outdir}/HoldoutMetrics_{model_name}_{tag}.csv",
        index=False
    )

    # Save predictions
    pd.DataFrame({
        "TrueLabel": y_test,
        "PredLabel": y_pred
    }).to_csv(f"{outdir}/Predictions_{model_name}_{tag}.csv", index=False)

    # Save probabilities
    pd.DataFrame(
        y_prob,
        columns=[f"{cls}_prob" for cls in classes]
    ).to_csv(f"{outdir}/Probabilities_{model_name}_{tag}.csv", index=False)

    # ----------------------------
    # CONFUSION MATRIX
    # ----------------------------
    cm = confusion_matrix(y_test, y_pred)
    cm_norm = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    plt.figure(figsize=(5,5))
    sns.heatmap(
        cm_norm, annot=True, fmt=".2f",
        cmap="crest", cbar=False,
        xticklabels=classes, yticklabels=classes,
        square=True
    )
    plt.title(f"{model_name} — {tag.upper()} Confusion Matrix")
    plt.tight_layout()
    plt.savefig(f"{outdir}/ConfusionMatrix_{model_name}_{tag}.png", dpi=600)
    plt.close()

    # ----------------------------
    # MULTICLASS ROC CURVES
    # ----------------------------
    y_bin = label_binarize(y_test, classes=classes)

    plt.figure(figsize=(6,5))
    for i, cls in enumerate(classes):
        fpr, tpr, _ = roc_curve(y_bin[:,i], y_prob[:,i])
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, lw=2, label=f"{cls} (AUC={roc_auc:.2f})")

    plt.plot([0,1],[0,1],"--",color="gray")
    plt.title(f"{model_name} — {tag.upper()} ROC Curves")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{outdir}/AUROC_{model_name}_{tag}.png", dpi=600)
    plt.close()

    print(f"✓ Completed Nested CV + Holdout for {model_name} ({tag})\n")


# ============================================================
# PARAMETER GRIDS
# ============================================================

gb_grid = {
    "n_estimators": [200, 300, 400],
    "learning_rate": [0.05, 0.1, 0.2],
    "max_depth": [2,3,4],
    "subsample": [0.8, 1.0]
}

rf_grid = {
    "n_estimators": [300, 500, 800],
    "max_depth": [None, 10, 20],
    "min_samples_leaf": [1, 3, 5]
}

lr_model = Pipeline([
    ("scaler", StandardScaler()),
    ("clf", LogisticRegression(max_iter=5000, class_weight="balanced", random_state=42))
])

lr_grid = {
    "clf__C": [0.1, 1, 10, 50]
}


# ============================================================
# RUN TOP-2 MODELS FOR SHARED
# ============================================================

print("\n=== SHARED DATASET: TOP 2 MODELS (GB, RF) ===\n")

nestedCV_model(shared_df, y,
               tag="shared",
               model_name="GB",
               model=GradientBoostingClassifier(random_state=42),
               param_grid=gb_grid)

nestedCV_model(shared_df, y,
               tag="shared",
               model_name="RF",
               model=RandomForestClassifier(random_state=42, n_jobs=-1, class_weight="balanced"),
               param_grid=rf_grid)


# ============================================================
# RUN TOP-2 MODELS FOR UNSHARED
# ============================================================

print("\n=== UNSHARED DATASET: TOP 2 MODELS (GB, LR) ===\n")

nestedCV_model(unshared_df, y,
               tag="unshared",
               model_name="GB",
               model=GradientBoostingClassifier(random_state=42),
               param_grid=gb_grid)

nestedCV_model(unshared_df, y,
               tag="unshared",
               model_name="LR",
               model=lr_model,
               param_grid=lr_grid)


print("\n ALL TOP-2 MODELS FOR SHARED & UNSHARED COMPLETED SUCCESSFULLY ")

EOF
