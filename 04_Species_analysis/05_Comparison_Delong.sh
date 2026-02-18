#!/bin/bash
#SBATCH --job-name=compare_GB_shared_vs_unshared
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/compare_GB_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/compare_GB_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 << 'EOF'
# ==============================================================
# COMPARISON OF SHARED-GB VS UNSHARED-GB PERFORMANCE
# Using DeLong’s Test (Correlated ROC Curves)
# ==============================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc
from scipy import stats

# =========================================================================================
# INPUT PATHS — UPDATED (THIS IS WHERE YOUR NESTED CV OUTPUTS ARE LOCATED)
# =========================================================================================

DATA_PATH = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/DATA/ML/data"

shared_prob_path   = f"{DATA_PATH}/Probabilities_GB_shared.csv"
unshared_prob_path = f"{DATA_PATH}/Probabilities_GB_unshared.csv"

shared_pred_path   = f"{DATA_PATH}/Predictions_GB_shared.csv"
unshared_pred_path = f"{DATA_PATH}/Predictions_GB_unshared.csv"

print("Using input data from:")
print(shared_prob_path)
print(unshared_prob_path)

# =========================================================================================
# LOAD PREDICTIONS & TRUE LABELS
# =========================================================================================

y_shared = pd.read_csv(shared_pred_path)["TrueLabel"].values
y_unshared = pd.read_csv(unshared_pred_path)["TrueLabel"].values

assert np.array_equal(y_shared, y_unshared), "Sample order mismatch between shared and unshared predictions!"

y_true = y_shared

# Binary mapping: Cancer = 1, Otherwise = 0
binary_y = (y_true == "Cancer").astype(int)

# Load predicted cancer probabilities
shared_prob   = pd.read_csv(shared_prob_path)["Cancer_prob"].values
unshared_prob = pd.read_csv(unshared_prob_path)["Cancer_prob"].values

print("\nLoaded GB predicted probabilities for Shared & Unshared.\n")

# =========================================================================================
# DELONG IMPLEMENTATION
# =========================================================================================

def compute_midrank(x):
    J = len(x)
    idx = np.argsort(x)
    sorted_x = x[idx]
    t = np.arange(J)
    midranks = t + 1
    return midranks[np.argsort(idx)]

def fastDeLong(predictions_sorted_transposed, label_1_count):
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m

    positive = predictions_sorted_transposed[:, :m]
    negative = predictions_sorted_transposed[:, m:]

    k = predictions_sorted_transposed.shape[0]
    tx = np.apply_along_axis(compute_midrank, 1, positive)
    ty = np.apply_along_axis(compute_midrank, 1, negative)
    tz = np.apply_along_axis(compute_midrank, 1, predictions_sorted_transposed)

    aucs = tz[:, :m].sum(axis=1) / (m*n) - (m+1)/(2*n)

    v01 = (tz[:, :m] - tx) / n
    v10 = (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    S = sx/m + sy/n

    return aucs, S

def delong_roc_test(labels, pred1, pred2):
    labels = np.array(labels)
    order = np.argsort(-(pred1 + pred2))
    labels = labels[order]
    pred1 = pred1[order]
    pred2 = pred2[order]

    label_1_count = sum(labels)

    preds = np.vstack([pred1, pred2])
    aucs, cov = fastDeLong(preds, label_1_count)

    auc_diff = aucs[0] - aucs[1]
    var = cov[0,0] + cov[1,1] - 2*cov[0,1]
    z = auc_diff / np.sqrt(var)
    p = 2 * (1 - stats.norm.cdf(abs(z)))

    return aucs[0], aucs[1], auc_diff, p, z

# =========================================================================================
# COMPUTE ROC CURVES AND AUCs
# =========================================================================================

fpr_s, tpr_s, _ = roc_curve(binary_y, shared_prob)
fpr_u, tpr_u, _ = roc_curve(binary_y, unshared_prob)

auc_s = auc(fpr_s, tpr_s)
auc_u = auc(fpr_u, tpr_u)

print(f"Shared GB AUC:   {auc_s:.4f}")
print(f"Unshared GB AUC: {auc_u:.4f}")

# =========================================================================================
# RUN DELONG TEST
# =========================================================================================

auc1, auc2, auc_diff, p_val, z = delong_roc_test(binary_y, shared_prob, unshared_prob)

print("\n================ DE-LONG TEST RESULTS ================")
print(f"AUC Shared-GB:      {auc1:.4f}")
print(f"AUC Unshared-GB:    {auc2:.4f}")
print(f"AUC Difference:     {auc_diff:.4f}")
print(f"Z-score:            {z:.4f}")
print(f"P-value:            {p_val:.6f}")
print("======================================================")

# Save results
pd.DataFrame({
    "AUC_Shared_GB":[auc1],
    "AUC_Unshared_GB":[auc2],
    "AUC_Difference":[auc_diff],
    "z_score":[z],
    "p_value":[p_val]
}).to_csv("AUC_Comparison_SharedGB_vs_UnsharedGB.csv", index=False)


# =========================================================================================
# ROC COMPARISON PLOT — PUBLICATION QUALITY
# =========================================================================================

plt.figure(figsize=(7,6), dpi=600)

plt.plot(fpr_s, tpr_s, lw=2.5, color="#1f77b4",
         label=f"Shared (GB) AUC = {auc_s:.3f}")

plt.plot(fpr_u, tpr_u, lw=2.5, color="#d62728",
         label=f"Unshared (GB) AUC = {auc_u:.3f}")

plt.plot([0,1],[0,1],"--",color="gray", lw=1)

sig_text = f"DeLong p = {p_val:.3e}"
plt.title(f"Shared-GB vs Unshared-GB\nROC Comparison\n{sig_text}",
          fontsize=13, weight="bold")

plt.xlabel("False Positive Rate", fontsize=12, weight="bold")
plt.ylabel("True Positive Rate", fontsize=12, weight="bold")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("ROC_Comparison_SharedGB_vs_UnsharedGB.png", dpi=600)
plt.close()


# =========================================================================================
# SCORE DISTRIBUTION PLOT (OPTIONAL)
# =========================================================================================

df_scores = pd.DataFrame({
    "Shared_GB": shared_prob,
    "Unshared_GB": unshared_prob,
    "Cancer": binary_y
})

plt.figure(figsize=(7,5), dpi=600)
sns.violinplot(data=df_scores.melt(id_vars="Cancer",
                                   value_vars=["Shared_GB","Unshared_GB"]),
               x="variable", y="value",
               palette=["#1f77b4","#d62728"])
plt.title("Distribution of Predicted Cancer Probabilities", weight="bold")
plt.ylabel("Cancer Probability", weight="bold")
plt.xlabel("")
plt.tight_layout()
plt.savefig("ScoreDistribution_SharedGB_vs_UnsharedGB.png", dpi=600)
plt.close()


print("\n All Done — AUC Comparison + DeLong Test + Plots Generated Successfully 🎉")

EOF
