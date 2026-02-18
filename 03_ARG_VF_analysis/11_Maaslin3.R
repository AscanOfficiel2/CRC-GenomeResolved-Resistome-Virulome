# ============================================================
# MaAsLin3 Differential Abundance Analysis
# PLS₁ + PLS₉ Species (Species TPM-based)
# ============================================================
library(maaslin3)
library(dplyr)
library(readr)
library(tibble)

# -----------------------------
# FILE PATHS
# -----------------------------
species_tpm_file <- "TPM_raw_filtered_bact.csv"
metadata_file    <- "CRC_metadata_aligned_bact.csv"
pls_species_file <- "PLS_1_9_Gene_to_Species_Link.csv"
output_dir       <- "Maaslin3_PLS1_PLS9"

dir.create(output_dir, showWarnings = FALSE)

# -----------------------------
# LOAD SPECIES TPM MATRIX
# (rows = species, cols = samples)
# -----------------------------
species_tpm_raw <- read_csv(
  species_tpm_file,
  show_col_types = FALSE
) %>% as.data.frame()

rownames(species_tpm_raw) <- species_tpm_raw[[1]]
species_tpm_raw <- species_tpm_raw[, -1]

# -----------------------------
# TRANSPOSE → samples as rows
# -----------------------------
species_tpm <- as.data.frame(t(species_tpm_raw))
species_tpm$Sample_ID <- rownames(species_tpm)
rownames(species_tpm) <- species_tpm$Sample_ID
species_tpm$Sample_ID <- NULL

cat("✔ Species TPM matrix shape:", dim(species_tpm), "\n")

# -----------------------------
# LOAD METADATA
# -----------------------------
meta <- read_csv(metadata_file, show_col_types = FALSE) %>%
  as.data.frame()

rownames(meta) <- meta$Sample_ID

# -----------------------------
# LOAD PLS₁ + PLS₉ SPECIES
# -----------------------------
pls_species <- read_csv(pls_species_file, show_col_types = FALSE) %>%
  distinct(Species, Component)

# -----------------------------
# FILTER TPM MATRIX
# -----------------------------
keep_species <- intersect(
  colnames(species_tpm),
  pls_species$Species
)

cat(" Number of PLS₁ + PLS₉ species:", length(keep_species), "\n")

species_tpm <- species_tpm[, keep_species, drop = FALSE]

# -----------------------------
# ALIGN SAMPLES
# -----------------------------
common_samples <- intersect(
  rownames(species_tpm),
  rownames(meta)
)

species_tpm <- species_tpm[common_samples, , drop = FALSE]
meta <- meta[common_samples, ]

cat("Samples used:", length(common_samples), "\n")

# -----------------------------
# FACTOR SETUP
# -----------------------------
meta$Group <- factor(
  meta$Group,
  levels = c("Healthy", "Adenoma", "Cancer")
)

# -----------------------------
# RUN MAASLIN3
# -----------------------------
fit <- maaslin3(
  input_data     = species_tpm,
  input_metadata = meta,
  output         = output_dir,
  
  formula        = "~ Group",
  
  normalization  = "NONE",
  transform      = "LOG",
  correction     = "BH",
  
  standardize    = FALSE,
  augment        = TRUE,
  
  median_comparison_abundance  = FALSE,
  median_comparison_prevalence = FALSE,
  
  max_significance = 0.25,
  cores = 1
)

cat("MAsLin3 completed successfully\n")
###############################################################################

# ============================================================
# MaAsLin3 Volcano Plots
# Separate PLS₁ and PLS₉
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# Canonical species function
# ------------------------------------------------------------
canonical_species <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all("[\\._-]", " ") %>%
    str_replace_all(" +", " ") %>%
    trimws()
}

# ------------------------------------------------------------
# FILE PATHS
# ------------------------------------------------------------
maaslin_file <- "Maaslin3_PLS1_PLS9/significant_results.tsv"
pls_file     <- "PLS_1_9_Gene_to_Species_Link.csv"
out_table    <- "Maaslin3_Significant_With_PLS_Component.csv"

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
maaslin <- read_tsv(maaslin_file, show_col_types = FALSE)
pls <- read_csv(pls_file, show_col_types = FALSE)

# ------------------------------------------------------------
# KEEP ONLY ABUNDANCE MODELS
# ------------------------------------------------------------
maaslin <- maaslin %>%
  filter(model == "abundance")

# ------------------------------------------------------------
# PREPARE MaAsLin3 RESULTS
# ------------------------------------------------------------
maaslin <- maaslin %>%
  rename(Species_raw = feature) %>%
  mutate(
    Species_canon = canonical_species(Species_raw),
    neg_log10_q = -log10(qval_individual)
  )

# ------------------------------------------------------------
# PREPARE PLS TABLE
# ------------------------------------------------------------
pls_species <- pls %>%
  distinct(Species, Component) %>%
  mutate(
    Species_canon = canonical_species(Species)
  )

# ------------------------------------------------------------
# MERGE
# ------------------------------------------------------------
merged <- maaslin %>%
  left_join(
    pls_species %>% select(Species_canon, Component),
    by = "Species_canon"
  )

write_csv(merged, out_table)

# ============================================================
# VOLCANO FUNCTION
# ============================================================
plot_component_volcano <- function(df, component, color, n_labels, outfile) {
  
  df_sub <- df %>% filter(Component == component)
  
  label_df <- df_sub %>%
    slice_max(
      order_by = abs(coef) * neg_log10_q,
      n = n_labels,
      with_ties = FALSE
    )
  
  p <- ggplot(df_sub, aes(x = coef, y = neg_log10_q)) +
    geom_point(
      fill = color,
      shape = 21,
      size = 3,
      stroke = 1,
      color = "black",
      alpha = 0.9
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_text_repel(
      data = label_df,
      aes(label = Species_raw),
      size = 2,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey50"
    ) +
    labs(
      title = paste0("Differentially Abundant Species — ", component),
      subtitle = "MaAsLin3 (TPM-based)",
      x = "Effect Size (log fold-change)",
      y = expression(-log[10]("FDR q-value"))
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, linewidth = 1.5)
    )
  
  ggsave(outfile, p, width = 5, height = 5, dpi = 600)
  print(p)
}

# ============================================================
# MAKE PLOTS
# ============================================================

plot_component_volcano(
  merged,
  "PLS_1",
  "#E74C3C",
  4,
  "Maaslin3_Volcano_PLS1_Annotated.png"
)

plot_component_volcano(
  merged,
  "PLS_9",
  "#2ECC71",
  4,
  "Maaslin3_Volcano_PLS9_Annotated.png"
)
###############################################################################

-----------------------------
  # LOAD LIBRARIES
  # -----------------------------
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -----------------------------
# LOAD DATA
# -----------------------------
merged <- read_csv(
  "Maaslin3_Significant_With_PLS_Component.csv",
  show_col_types = FALSE
)

# -----------------------------
# FILTER SIGNIFICANT ABUNDANCE RESULTS
# -----------------------------
df_sig <- merged %>%
  filter(
    model == "abundance",
    qval_individual <= 0.05,
    !is.na(coef)
  )

# -----------------------------
# SELECT TOP SPECIES PER PLS AXIS
# -----------------------------
top_pls1 <- df_sig %>%
  filter(Component == "PLS_1") %>%
  group_by(Species_raw) %>%
  summarise(max_abs = max(abs(coef)), .groups = "drop") %>%
  slice_max(max_abs, n = 40) %>%
  pull(Species_raw)

top_pls9 <- df_sig %>%
  filter(Component == "PLS_9") %>%
  group_by(Species_raw) %>%
  summarise(max_abs = max(abs(coef)), .groups = "drop") %>%
  slice_max(max_abs, n = 25) %>%
  pull(Species_raw)

df_keep <- df_sig %>%
  filter(
    (Component == "PLS_1" & Species_raw %in% top_pls1) |
      (Component == "PLS_9" & Species_raw %in% top_pls9)
  )

# -----------------------------
# COLLAPSE TO ONE ROW PER SPECIES × PLS × CONTRAST
# -----------------------------
heat_df <- df_keep %>%
  mutate(
    Contrast = case_when(
      value == "Cancer"  ~ "Cancer",
      value == "Adenoma" ~ "Adenoma"
    )
  ) %>%
  group_by(Species_raw, Component, Contrast) %>%
  summarise(
    coef = coef[which.max(abs(coef))],
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from  = Contrast,
    values_from = coef
  )

# -----------------------------
# ORDER SPECIES WITHIN EACH PLS BLOCK
# -----------------------------
heat_df <- heat_df %>%
  group_by(Component) %>%
  mutate(
    max_abs = pmax(
      abs(Cancer),
      abs(Adenoma),
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  arrange(Component, desc(max_abs))

# -----------------------------
# CREATE UNIQUE ROW IDS (PLS + SPECIES)
# -----------------------------
heat_df <- heat_df %>%
  mutate(
    row_id = paste(Component, Species_raw, sep = " | ")
  )

# -----------------------------
# BUILD HEATMAP MATRIX
# -----------------------------
heat_mat <- heat_df %>%
  select(row_id, Cancer, Adenoma) %>%
  column_to_rownames("row_id") %>%
  as.matrix()

# Row annotations
row_pls     <- heat_df$Component
row_species <- heat_df$Species_raw

# -----------------------------
# LEFT ROW ANNOTATION (PLS AXIS)
# -----------------------------
row_anno <- rowAnnotation(
  PLS = row_pls,
  col = list(
    PLS = c(
      "PLS_1" = "#D7191C",  # red
      "PLS_9" = "#1B9E77"   # green
    )
  ),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold"),
  annotation_legend_param = list(
    title = "PLS axis",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# -----------------------------
# COLOR SCALE (SHARED, SYMMETRIC)
# -----------------------------
lim <- max(abs(heat_mat), na.rm = TRUE)

col_fun <- colorRamp2(
  c(-lim, 0, lim),
  c("#2C7BB6", "white", "#D7191C")
)

# -----------------------------
# MAIN HEATMAP
# -----------------------------
ht <- Heatmap(
  heat_mat,
  name               = "β",
  col                = col_fun,
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  show_row_dend      = FALSE,
  show_column_dend   = FALSE,
  na_col             = "white",
  rect_gp            = gpar(col = "grey40", lwd = 0.6),
  row_labels         = row_species,
  row_names_gp       = gpar(fontsize = 7, fontface = "italic"),
  column_names_gp    = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title     = expression(beta),
    at        = c(-round(lim, 1), 0, round(lim, 1)),
    labels    = c("↓ Healthy", "0", "↑ Disease"),
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

# -----------------------------
# SAVE AS HIGH-RES PNG (600 DPI)
# -----------------------------
png(
  filename = "Maaslin3_PLS1_PLS9_Annotated_Heatmap.png",
  width    = 2600,
  height   = 4200,
  res      = 600
)

draw(
  row_anno + ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "left",
  padding = unit(c(6, 6, 6, 6), "mm")
)

dev.off()

