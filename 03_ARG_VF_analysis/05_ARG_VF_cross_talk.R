### ============================================================
### Load required packages
### ============================================================
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

### ============================================================
### 1. Load data
### ============================================================
arg <- read.csv("arg_cleaned_final_with_mechanisms.csv")
vf  <- read.csv("vfdb_with_species_cleaned.csv")
meta <- read.csv("ARG_VF_metadata_final_matched.csv")

### ============================================================
### 2. Prepare ARG and VF tables
### ============================================================
arg2 <- arg %>% 
  select(Sample_ID, Gene, TPM, Species, ARG_Mechanism) %>%
  mutate(Type = "ARG") %>%
  rename(Mechanism = ARG_Mechanism)

vf2 <- vf %>%
  select(Sample_ID, Gene, TPM, Species, Mechanism) %>%
  mutate(Type = "VF")

combined <- bind_rows(arg2, vf2)

combined <- combined %>%
  left_join(meta %>% select(Sample_ID, Group), by="Sample_ID")

### ============================================================
### 3. Identify species with both ARGs and VFs
### ============================================================
species_ARG <- arg2 %>%
  group_by(Species) %>%
  summarise(ARG_total = sum(TPM)) %>%
  filter(ARG_total > 0)

species_VF <- vf2 %>%
  group_by(Species) %>%
  summarise(VF_total = sum(TPM)) %>%
  filter(VF_total > 0)

species_both <- species_ARG %>%
  inner_join(species_VF, by = "Species") %>%
  mutate(ARG_VF_ratio = ARG_total / VF_total)

print("Species with both ARG and VF:")
print(species_both)

### ============================================================
### 4. Species-level ARG–VF co-expression correlation
### ============================================================

### Build aligned ARG matrix
arg_species_mat <- arg2 %>%
  group_by(Species, Sample_ID) %>%
  summarise(ARG_TPM = sum(TPM), .groups="drop") %>%
  pivot_wider(names_from = Sample_ID, values_from = ARG_TPM, values_fill = 0)

### Build aligned VF matrix
vf_species_mat <- vf2 %>%
  group_by(Species, Sample_ID) %>%
  summarise(VF_TPM = sum(TPM), .groups="drop") %>%
  pivot_wider(names_from = Sample_ID, values_from = VF_TPM, values_fill = 0)

### Align species
common_species <- intersect(arg_species_mat$Species, vf_species_mat$Species)
arg_sub <- arg_species_mat %>% filter(Species %in% common_species)
vf_sub  <- vf_species_mat %>% filter(Species %in% common_species)

### Align sample IDs
all_samples <- union(colnames(arg_sub)[-1], colnames(vf_sub)[-1])

arg_m <- arg_sub %>%
  pivot_longer(-Species, names_to="Sample_ID", values_to="ARG_TPM") %>%
  complete(Species, Sample_ID = all_samples, fill=list(ARG_TPM=0)) %>%
  pivot_wider(names_from=Sample_ID, values_from=ARG_TPM) %>%
  arrange(match(Species, common_species)) %>%
  select(-Species) %>% as.matrix()

vf_m <- vf_sub %>%
  pivot_longer(-Species, names_to="Sample_ID", values_to="VF_TPM") %>%
  complete(Species, Sample_ID = all_samples, fill=list(VF_TPM=0)) %>%
  pivot_wider(names_from=Sample_ID, values_from=VF_TPM) %>%
  arrange(match(Species, common_species)) %>%
  select(-Species) %>% as.matrix()

### Compute species-level ARG–VF correlation
cor_scores <- sapply(1:length(common_species), function(i) {
  cor(arg_m[i, ], vf_m[i, ], method="spearman")
})

ARG_VF_correlation <- tibble(
  Species = common_species,
  ARG_VF_Correlation = cor_scores
)

print("ARG–VF correlation per species:")
print(ARG_VF_correlation)

### ============================================================
### 5. ARG/VF Ratio Table (pathogenicity balance)
### ============================================================
ARG_VF_ratio_table <- species_both %>%
  arrange(desc(ARG_VF_ratio))

print("ARG/VF Ratio Table:")
print(ARG_VF_ratio_table)

### ============================================================
### 6. ARG–VF Mechanism Crosstalk Matrix
### ============================================================

arg_mech_species <- arg2 %>%
  group_by(Species, Mechanism) %>%
  summarise(ARG_Mech_TPM = sum(TPM), .groups="drop") %>%
  pivot_wider(names_from=Mechanism, values_from=ARG_Mech_TPM, values_fill=0)

vf_mech_species <- vf2 %>%
  group_by(Species, Mechanism) %>%
  summarise(VF_Mech_TPM = sum(TPM), .groups="drop") %>%
  pivot_wider(names_from=Mechanism, values_from=VF_Mech_TPM, values_fill=0)

common_mech_species <- intersect(arg_mech_species$Species, vf_mech_species$Species)

arg_mech_mat <- arg_mech_species %>%
  filter(Species %in% common_mech_species) %>%
  select(-Species) %>% as.matrix()

vf_mech_mat <- vf_mech_species %>%
  filter(Species %in% common_mech_species) %>%
  select(-Species) %>% as.matrix()

### Remove zero-variance columns
arg_mech_mat <- arg_mech_mat[, apply(arg_mech_mat, 2, sd) > 0, drop=FALSE]
vf_mech_mat  <- vf_mech_mat[, apply(vf_mech_mat, 2, sd) > 0, drop=FALSE]

### Compute mechanism crosstalk correlation
crosstalk_matrix <- cor(arg_mech_mat, vf_mech_mat, method="spearman")
crosstalk_matrix[!is.finite(crosstalk_matrix)] <- 0

print("Mechanism Crosstalk Matrix:")
print(crosstalk_matrix)

### ============================================================
### 7. High-Risk Pathobiont Score
### ============================================================
high_risk <- species_both %>%
  left_join(ARG_VF_correlation, by="Species") %>%
  mutate(
    risk_score = scale(ARG_total) +
      scale(VF_total) +
      scale(ARG_VF_Correlation)
  ) %>%
  arrange(desc(risk_score))

print("High-risk ARG+VF Pathobionts:")
print(high_risk)

### ============================================================
### SAVE OUTPUT TABLES
### ============================================================

write.csv(species_both, "Species_with_ARG_and_VF.csv", row.names = FALSE)
write.csv(ARG_VF_ratio_table, "ARG_VF_Ratio_Table.csv", row.names = FALSE)
write.csv(high_risk, "High_Risk_Pathobionts.csv", row.names = FALSE)

cat("Saved: Species_with_ARG_and_VF.csv\n")
cat("Saved: ARG_VF_Ratio_Table.csv\n")
cat("Saved: High_Risk_Pathobionts.csv\n")

### ============================================================
### pheatmap for Crosstalk (clean version)
### ============================================================

crosstalk_clean <- crosstalk_matrix
crosstalk_clean[!is.finite(crosstalk_clean)] <- 0

breaks <- seq(-1, 1, length.out = 200)
colors <- colorRampPalette(c("blue", "white", "red"))(199)
colors <- viridisLite::plasma(200)
colors <- viridisLite::magma(200)


# ----------------------------------------
# ARG–VF Mechanism Crosstalk Heatmap
# ----------------------------------------

png(
  "ARG_VF_Crosstalk_pheatmap.png",
  width = 3500,
  height = 3800,
  res = 600
)

pheatmap(
  crosstalk_clean,
  color = viridisLite::magma(200),
  breaks = seq(-1, 1, length.out = 200),
  border_color = NA,
  clustering_method = "complete",
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 9,
  main = "ARG–VF Mechanism Crosstalk",
  angle_col = "45"   # FIXED
)

dev.off()

cat("Saved: ARG_VF_Crosstalk_pheatmap.png (600 dpi)\n")

################################################################################
###############################################################################
### OPTIONAL — Save Long-Format Spearman Correlation Table
###############################################################################

crosstalk_long <- as.data.frame(as.table(crosstalk_clean)) %>%
  rename(
    ARG_Mechanism = Var1,
    VF_Mechanism  = Var2,
    Spearman_r    = Freq
  ) %>%
  arrange(desc(abs(Spearman_r)))

write.csv(
  crosstalk_long,
  "ARG_VF_Crosstalk_Spearman_LongFormat.csv",
  row.names = FALSE
)

cat("Saved: ARG_VF_Crosstalk_Spearman_LongFormat.csv\n")
###################################################################################
