# ============================================================
# Mechanism Enrichment & Activity Analysis
# PLS₁ vs PLS₉ (ARG + VF)
# FINAL — REVIEWER-SAFE VERSION
# ============================================================
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(scales)

# -----------------------------
# FILE PATHS
# -----------------------------
tpm_file        <- "ARG_VF_TPM_raw_filtered_genes.csv"
meta_file       <- "ARG_VF_metadata_final_matched.csv"
pls_file        <- "PLS_TopGenes_Ranked_with_ComponentStats.csv"
arg_mech_file   <- "arg_cleaned_final_with_mechanisms.csv"
vf_mech_file    <- "vfdb_with_species_cleaned.csv"

# -----------------------------
# HELPER: sanitize gene names (JOIN ONLY)
# -----------------------------
sanitize_gene <- function(x) {
  x %>%
    tolower() %>%
    str_replace("^(arg_|vf_)", "") %>%
    str_replace_all("[\\(\\)\\-\\.]", "_") %>%
    str_replace_all("__+", "_") %>%
    trimws()
}

# ============================================================
# LOAD TPM MATRIX
# ============================================================
tpm <- read_csv(tpm_file, show_col_types = FALSE) %>% as.data.frame()
rownames(tpm) <- tpm[[1]]
tpm <- tpm[, -1, drop = FALSE]
cat("✔ TPM matrix:", dim(tpm), "\n")

# ============================================================
# LOAD METADATA
# ============================================================
meta <- read_csv(meta_file, show_col_types = FALSE) %>%
  as.data.frame() %>%
  filter(Sample_ID %in% rownames(tpm))

rownames(meta) <- meta$Sample_ID
meta$Group <- factor(meta$Group, levels = c("Healthy", "Adenoma", "Cancer"))

# Align samples
common_samples <- intersect(rownames(tpm), meta$Sample_ID)
tpm  <- tpm[common_samples, , drop = FALSE]
meta <- meta[common_samples, ]
stopifnot(all(rownames(tpm) == meta$Sample_ID))
cat(" Samples used:", length(common_samples), "\n")

# ============================================================
# LOAD PLS TOP GENES (AUTHORITATIVE)
# ============================================================
pls_genes <- read_csv(pls_file, show_col_types = FALSE) %>%
  filter(Component %in% c("PLS_1", "PLS_9")) %>%
  transmute(
    Component,
    Gene_full = Gene,
    Gene_core = sanitize_gene(Gene)
  ) %>%
  distinct()

cat(" PLS genes:", nrow(pls_genes), "\n")

# ============================================================
# LOAD ARG MECHANISMS
# ============================================================
arg_mech <- read_csv(arg_mech_file, show_col_types = FALSE) %>%
  transmute(
    Gene_core = sanitize_gene(Gene),
    Mechanism = ARG_Mechanism,
    Gene_type = "ARG"
  ) %>%
  distinct()

cat("✔ ARG mechanisms:", nrow(arg_mech), "\n")

# ============================================================
# LOAD VF MECHANISMS  CRITICAL FIX
# (Your VF file is TAB-separated but named .csv)
# ============================================================
vf_raw <- read_tsv(vf_mech_file, show_col_types = FALSE)

stopifnot(all(c("Gene", "Mechanism") %in% colnames(vf_raw)))

vf_mech <- vf_raw %>%
  transmute(
    Gene_core = sanitize_gene(Gene),
    Mechanism = Mechanism,
    Gene_type = "VF"
  ) %>%
  distinct()

cat(" VF mechanisms:", nrow(vf_mech), "\n")

# ============================================================
# COMBINE MECHANISM TABLES
# ============================================================
mech_map <- bind_rows(arg_mech, vf_mech)

# ============================================================
# MAP PLS GENES → MECHANISMS
#  AUTHORITATIVE CONSTRAINT ENFORCED HERE
# ============================================================
pls_mech <- pls_genes %>%
  left_join(mech_map, by = "Gene_core") %>%
  filter(
    !(Component == "PLS_1" & Gene_type == "VF")
  )

cat(" PLS genes with mechanisms:",
    sum(!is.na(pls_mech$Mechanism)), "\n")

# ============================================================
# FILTER TPM TO PLS GENES
# ============================================================
tpm_pls <- tpm[, intersect(colnames(tpm), pls_mech$Gene_full), drop = FALSE]
cat(" TPM genes used:", ncol(tpm_pls), "\n")

# ============================================================
# LONG FORMAT TPM
# ============================================================
tpm_long <- tpm_pls %>%
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(-Sample_ID, names_to = "Gene_full", values_to = "TPM") %>%
  left_join(meta, by = "Sample_ID") %>%
  left_join(pls_mech, by = "Gene_full")

# ============================================================
# A) MECHANISM ENRICHMENT (COMPOSITION)
# ============================================================
mech_enrichment <- pls_mech %>%
  filter(!is.na(Mechanism)) %>%
  count(Component, Gene_type, Mechanism)

write_csv(
  mech_enrichment,
  "PLS1_PLS9_Mechanism_Enrichment.csv"
)

# ============================================================
# B) MECHANISM ACTIVITY ACROSS GROUPS
# ============================================================
mech_activity <- tpm_long %>%
  filter(!is.na(Mechanism)) %>%
  group_by(Component, Gene_type, Mechanism, Group) %>%
  summarise(
    Mean_TPM = mean(TPM, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  mech_activity,
  "PLS1_PLS9_Mechanism_Activity_ByGroup.csv"
)

cat("Mechanism enrichment & activity analysis COMPLETE\n")

#################################################################################
library(dplyr)
library(readr)
library(ggplot2)
library(scales)

# -----------------------------
# INPUT FILE
# -----------------------------
mech_file <- "PLS1_PLS9_Mechanism_Enrichment.csv"

# -----------------------------
# LOAD DATA
# -----------------------------
df <- read_csv(mech_file, show_col_types = FALSE)

# -----------------------------
# PREPARE DATA
# -----------------------------
df_plot <- df %>%
  mutate(
    Mechanism_Label = paste0(Gene_type, " – ", Mechanism)
  ) %>%
  group_by(Component) %>%
  mutate(
    total = sum(n),
    prop  = n / total
  ) %>%
  ungroup()

# -----------------------------
# DEFINE COLORS (MATCH FILE CONTENTS)
# -----------------------------
mech_colors <- c(
  # -------- ARG --------
  "ARG – antibiotic efflux"            = "#4C72B0",
  "ARG – antibiotic inactivation"      = "#55A868",
  "ARG – target alteration"            = "#C44E52",
  "ARG – target protection"            = "#E17C7C",
  "ARG – target replacement"           = "#B07AA1",
  "ARG – reduced permeability"         = "#6A8EC7",
  "ARG – Efflux; reduced permeability" = "#3B5B92",
  
  # -------- VF --------
  "VF – Adherence"           = "#DD8452",
  "VF – Immune modulation"   = "#937860",
  "VF – Invasion"            = "#DA8BC3",
  "VF – Motility"            = "#8C8C8C",
  "VF – Regulation"          = "#9D9D9D",
  "VF – Exoenzyme"           = "#CC79A7"
)

# -----------------------------
# FUNCTION: BARPLOT PER COMPONENT
# -----------------------------
plot_mech_bar <- function(component_name, outfile) {
  
  df_sub <- df_plot %>%
    filter(Component == component_name) %>%
    arrange(prop)
  
  p <- ggplot(
    df_sub,
    aes(
      x = prop,
      y = reorder(Mechanism_Label, prop),
      fill = Mechanism_Label
    )
  ) +
    geom_col(
      width = 0.7,
      color = "black",
      linewidth = 0.4
    ) +
    scale_fill_manual(values = mech_colors, drop = FALSE) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste0(component_name, " Mechanism Enrichment"),
      subtitle = "ARG and VF mechanisms (gene proportions)",
      x = "Proportion of genes",
      y = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.title.x = element_text(face = "bold"),
      legend.position = "none",
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 1
      )
    )
  
  # ONLY CHANGE IS HERE
  ggsave(
    filename = outfile,
    plot = p,
    width = 5,
    height = 5,
    dpi = 600,
    bg = "white"
  )
  
  print(p)
}

# -----------------------------
# DRAW BAR PLOTS
# -----------------------------
plot_mech_bar(
  component_name = "PLS_1",
  outfile = "PLS1_Mechanism_Enrichment_Barplot.png"
)

plot_mech_bar(
  component_name = "PLS_9",
  outfile = "PLS9_Mechanism_Enrichment_Barplot.png"
)

cat("PLS₁ and PLS₉ mechanism enrichment bar plots saved (5 × 5 inches)\n")
######################################################################################

