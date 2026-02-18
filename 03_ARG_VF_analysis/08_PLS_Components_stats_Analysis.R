############################################################
# PLS → Mechanism Analysis (ARG + VF)
# PLS Components: PLS_1 to PLS_10 combined
# Stats: Kruskal–Wallis (raw p) + Dunn post-hoc
############################################################

library(dplyr)
library(tidyr)
library(stringr)
library(FSA)

# -----------------------------
# Load Files
# -----------------------------
PLS  <- read.csv("PLS_TopGenes_Ranked_with_ComponentStats.csv", stringsAsFactors = FALSE)
ARG  <- read.csv("arg_cleaned_final_with_mechanisms.csv", stringsAsFactors = FALSE)
VF   <- read.delim("vfdb_with_species_cleaned.csv", stringsAsFactors = FALSE)
TPM  <- read.csv("ARG_VF_TPM_raw_filtered_genes.csv", stringsAsFactors = FALSE, check.names = FALSE)
META <- read.csv("ARG_VF_metadata_final_matched.csv", stringsAsFactors = FALSE)

# -----------------------------
# Fix TPM first column
# -----------------------------
colnames(TPM)[1] <- "Sample_ID"
TPM$Sample_ID <- as.character(TPM$Sample_ID)

# -----------------------------
# Filter PLS Genes (PLS_1 to PLS_10)
# -----------------------------
PLS <- PLS %>%
  dplyr::filter(Component %in% paste0("PLS_", 1:10)) %>%
  dplyr::mutate(Gene_clean = str_remove(Gene, "^(ARG|VF)_"))

PLS_genes <- unique(PLS$Gene)

# -----------------------------
# Metadata
# -----------------------------
META <- META %>%
  dplyr::select(Sample_ID, Group) %>%
  dplyr::mutate(Group = factor(Group, levels = c("Healthy", "Adenoma", "Cancer")))

# -----------------------------
# ARG Mechanism Map
# -----------------------------
arg_mech <- ARG %>%
  dplyr::select(Gene, ARG_Mechanism) %>%
  dplyr::rename(Mechanism = ARG_Mechanism) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    Gene_clean = stringr::str_remove(Gene, "^(ARG|VF)_"),
    Source = "ARG"
  )

# -----------------------------
# VF Mechanism Map
# -----------------------------
colnames(VF) <- trimws(colnames(VF))

vf_mech <- VF %>%
  dplyr::select(Gene, Mechanism) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(Mechanism)) %>%
  dplyr::mutate(
    Gene_clean = stringr::str_remove(Gene, "^(ARG|VF)_"),
    Source = "VF"
  )

# -----------------------------
# Combine Mechanisms
# -----------------------------
mech_map <- bind_rows(arg_mech, vf_mech)

# -----------------------------
# Filter TPM and melt
# -----------------------------
TPM <- TPM %>%
  dplyr::semi_join(META, by = "Sample_ID") %>%
  dplyr::select(Sample_ID, dplyr::any_of(PLS_genes))

tpm_long <- TPM %>%
  tidyr::pivot_longer(
    cols = -Sample_ID,
    names_to = "Gene",
    values_to = "Abundance"
  ) %>%
  dplyr::mutate(
    Gene_clean = stringr::str_remove(Gene, "^(ARG|VF)_"),
    Abundance = as.numeric(Abundance)
  )

# -----------------------------
# Merge all
# -----------------------------
merged <- tpm_long %>%
  dplyr::left_join(
    mech_map,
    by = "Gene_clean",
    relationship = "many-to-many"
  ) %>%
  dplyr::left_join(META, by = "Sample_ID") %>%
  dplyr::filter(!is.na(Mechanism))

# -----------------------------
# Sum TPM per Sample × Mechanism
# -----------------------------
sample_mech <- merged %>%
  dplyr::group_by(Sample_ID, Group, Mechanism) %>%
  dplyr::summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  sample_mech,
  "PLS_1_to_10_Mechanism_Abundance_PerSample.csv",
  row.names = FALSE
)

# -----------------------------
# Kruskal–Wallis per Mechanism
# -----------------------------
kw_results <- sample_mech %>%
  dplyr::group_by(Mechanism) %>%
  dplyr::summarise(
    KW_pval = tryCatch(
      kruskal.test(Abundance ~ Group)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  )

write.csv(
  kw_results,
  "PLS_1_to_10_Mechanism_KW.csv",
  row.names = FALSE
)

# -----------------------------
# Dunn post-hoc for significant KW (p < 0.05)
# -----------------------------
sig_mechs <- kw_results %>%
  dplyr::filter(KW_pval < 0.05) %>%
  dplyr::pull(Mechanism)

dunn_all <- list()

for (mech in sig_mechs) {
  df <- sample_mech %>% dplyr::filter(Mechanism == mech)
  
  res <- tryCatch({
    d <- FSA::dunnTest(Abundance ~ Group, data = df, method = "bh")$res
    d$Mechanism <- mech
    d
  }, error = function(e) NULL)
  
  if (!is.null(res)) dunn_all[[mech]] <- res
}

if (length(dunn_all) > 0) {
  dunn_df <- bind_rows(dunn_all)
  write.csv(
    dunn_df,
    "PLS_1_to_10_Mechanism_Dunn.csv",
    row.names = FALSE
  )
}

cat("\n SUCCESS: PLS_1 to PLS_10 mechanism analysis complete.\n")
######################################################################################
#####################################################################################
############################################################
# Mechanism × Group Heatmap (PLS_1–PLS_10)
# Rows = Mechanisms | Columns = Healthy, Adenoma, Cancer
# Output: Publication-ready PNG
############################################################
library(pheatmap)
library(dplyr)
library(tidyr)
library(readr)
library(RColorBrewer)
library(tibble)

############################################################
# Files & Parameters
############################################################

abund_file <- "PLS_1_to_10_Mechanism_Abundance_PerSample.csv"
meta_file  <- "ARG_VF_metadata_final_matched.csv"
output_png <- "PLS_1_to_10_Mechanism_Group_Heatmap.png"

# TRUE = row-wise Z-score (recommended for main figure)
# FALSE = raw mean abundance
scale_rows <- TRUE

# Desired column order
group_order <- c("Healthy", "Adenoma", "Cancer")

############################################################
# Load Data
############################################################

abund <- read.csv(abund_file, stringsAsFactors = FALSE)
meta  <- read.csv(meta_file,  stringsAsFactors = FALSE)

############################################################
# ROBUST METADATA HANDLING (Group vs group)
############################################################

if ("Group" %in% colnames(meta)) {
  meta_clean <- meta[, c("Sample_ID", "Group")]
} else if ("group" %in% colnames(meta)) {
  meta_clean <- meta
  colnames(meta_clean)[colnames(meta_clean) == "group"] <- "Group"
  meta_clean <- meta_clean[, c("Sample_ID", "Group")]
} else {
  stop("Metadata must contain a column named 'Group' or 'group'")
}

############################################################
# CRITICAL FIX: REMOVE Group FROM ABUND BEFORE JOIN
############################################################

# This prevents Group.x / Group.y collisions
if ("Group" %in% colnames(abund)) {
  abund <- abund %>% dplyr::select(-Group)
}

# Ensure Sample_IDs match exactly
abund$Sample_ID <- trimws(abund$Sample_ID)
meta_clean$Sample_ID <- trimws(meta_clean$Sample_ID)

############################################################
# Merge Group Information (NOW SAFE)
############################################################

abund <- abund %>%
  dplyr::left_join(meta_clean, by = "Sample_ID") %>%
  dplyr::filter(!is.na(Group)) %>%
  dplyr::mutate(
    Group = factor(Group, levels = group_order)
  )

############################################################
# Compute Mean Abundance per Mechanism × Group
############################################################

mech_group_matrix <- abund %>%
  dplyr::group_by(Mechanism, Group) %>%
  dplyr::summarise(
    Mean_Abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = Group,
    values_from = Mean_Abundance,
    values_fill = 0
  ) %>%
  tibble::column_to_rownames("Mechanism") %>%
  as.matrix()

############################################################
# Enforce Column Order (Healthy → Adenoma → Cancer)
############################################################

mech_group_matrix <- mech_group_matrix[, group_order]

############################################################
# Optional Row-wise Z-score Scaling
############################################################

if (scale_rows) {
  mech_group_matrix <- t(scale(t(mech_group_matrix)))
  mech_group_matrix[is.na(mech_group_matrix)] <- 0
}

############################################################
# Color Palette (Clean, Journal-Safe)
############################################################

heat_colors <- colorRampPalette(
  rev(brewer.pal(n = 11, name = "RdBu"))
)(100)

############################################################
# Plot and Save Heatmap
############################################################

png(
  filename = output_png,
  width    = 1600,
  height   = 1200,
  res      = 300
)

pheatmap(
  mat = mech_group_matrix,
  color = heat_colors,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",        # scaling already handled
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 11,
  angle_col = 45,  
  main = "Mechanism Abundance by Group (PLS_1–PLS_10)",
  display_numbers = FALSE,
  labels_col = group_order
)

dev.off()

cat("\n Heatmap successfully saved as:\n", output_png, "\n")
###############################################################################
############################################################
# Link PLS_1 and PLS_9 Genes to Species using TPM data
# Output: Table of Gene × Component × Species
############################################################

library(dplyr)
library(readr)
library(stringr)

# -----------------------------
# Load Files
# -----------------------------
pls <- read.csv("PLS_TopGenes_Ranked_with_ComponentStats.csv", stringsAsFactors = FALSE)
arg <- read.csv("arg_cleaned_final_with_mechanisms.csv", stringsAsFactors = FALSE)
vf  <- read.delim("vfdb_with_species_cleaned.csv", stringsAsFactors = FALSE)  # TAB-delimited

# -----------------------------
# Filter PLS_1 and PLS_9 Genes
# -----------------------------
pls <- pls %>%
  filter(Component %in% c("PLS_1", "PLS_9")) %>%
  mutate(Gene_clean = str_remove(Gene, "^(ARG|VF)_"))

# -----------------------------
# Extract ARG + VF Gene TPM + Species Data
# -----------------------------

# Fix column names
colnames(arg) <- trimws(colnames(arg))
colnames(vf)  <- trimws(colnames(vf))

# Ensure both have Sample_ID, Gene, TPM, Species
arg_sub <- arg %>%
  select(Sample_ID, Gene, TPM, Species) %>%
  mutate(Source = "ARG")

vf_sub <- vf %>%
  select(Sample_ID, Gene, TPM, Species) %>%
  mutate(Source = "VF")

tpm_long <- bind_rows(arg_sub, vf_sub) %>%
  mutate(Gene_clean = str_remove(Gene, "^(ARG|VF)_"))

# -----------------------------
# Join PLS info
# -----------------------------
linked <- tpm_long %>%
  inner_join(pls %>% select(Gene_clean, Component), by = "Gene_clean")

# -----------------------------
# Summarize: Gene × Component × Species
# -----------------------------
summary_table <- linked %>%
  group_by(Gene_clean, Component, Species) %>%
  summarise(
    Sample_Count = n_distinct(Sample_ID),
    Mean_TPM = mean(TPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Component, desc(Mean_TPM))

# -----------------------------
# Save Output
# -----------------------------
write.csv(summary_table, "PLS_1_9_Gene_to_Species_Link.csv", row.names = FALSE)

cat("\n Gene-to-species mapping saved as: PLS_1_9_Gene_to_Species_Link.csv\n")
#############################################################################################################
############################################################
# Link All PLS_1 to PLS_10 Genes to Species Using TPM Data
# Output: Gene × Component × Species mapping with expression stats
############################################################

library(dplyr)
library(stringr)
library(readr)

# -----------------------------
# Load PLS Gene–Component Table
# -----------------------------
pls <- read.csv("PLS_TopGenes_Ranked_with_ComponentStats.csv", stringsAsFactors = FALSE)

# Filter for PLS_1 to PLS_10 only
pls <- pls %>%
  filter(Component %in% paste0("PLS_", 1:10)) %>%
  mutate(Gene_clean = str_remove(Gene, "^(ARG|VF)_"))

# -----------------------------
# Load ARG + VF files with TPM + Species
# -----------------------------
arg <- read.csv("arg_cleaned_final_with_mechanisms.csv", stringsAsFactors = FALSE)
vf  <- read.delim("vfdb_with_species_cleaned.csv", stringsAsFactors = FALSE)  # tab-delimited

# Fix column names
colnames(arg) <- trimws(colnames(arg))
colnames(vf)  <- trimws(colnames(vf))

# Subset needed columns (ensure all exist)
arg_sub <- arg %>%
  select(Sample_ID, Gene, TPM, Species) %>%
  mutate(Source = "ARG")

vf_sub <- vf %>%
  select(Sample_ID, Gene, TPM, Species) %>%
  mutate(Source = "VF")

# -----------------------------
# Combine ARG + VF into one long-format TPM table
# -----------------------------
tpm_long <- bind_rows(arg_sub, vf_sub) %>%
  mutate(Gene_clean = str_remove(Gene, "^(ARG|VF)_"))

# -----------------------------
# Join with PLS Component Info
# -----------------------------
linked <- tpm_long %>%
  inner_join(pls %>% select(Gene_clean, Component), by = "Gene_clean")

# -----------------------------
# Summarize: Gene × Component × Species
# -----------------------------
summary_table <- linked %>%
  group_by(Gene_clean, Component, Species) %>%
  summarise(
    Sample_Count = n_distinct(Sample_ID),
    Mean_TPM     = mean(TPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Component, desc(Mean_TPM))

# -----------------------------
# Save Output
# -----------------------------
write.csv(
  summary_table,
  "PLS_1_to_10_Gene_to_Species_Link.csv",
  row.names = FALSE
)

cat("\n Gene-to-species mapping saved as: PLS_1_to_10_Gene_to_Species_Link.csv\n")

