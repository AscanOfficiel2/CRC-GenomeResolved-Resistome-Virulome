# ============================================================
# Pairwise log2 Fold Change for PLS₁ and PLS₉ Genes
# Using PLS Top Genes
# DESCRIPTIVE ONLY
# ============================================================

# -----------------------------
# LIBRARIES
# -----------------------------
library(dplyr)
library(readr)
library(tibble)
library(pheatmap)
library(stringr)

# -----------------------------
# FILE PATHS
# -----------------------------
gene_tpm_file <- "ARG_VF_TPM_raw_filtered_genes.csv"
meta_file     <- "ARG_VF_metadata_final_matched.csv"
pls_file      <- "PLS_TopGenes_Ranked_with_ComponentStats.csv"

# -----------------------------
# PARAMETERS
# -----------------------------
pseudocount <- 1e-6

# ============================================================
# LOAD GENE TPM MATRIX
# ============================================================
gene_tpm <- read_csv(gene_tpm_file, show_col_types = FALSE) %>%
  as.data.frame()

rownames(gene_tpm) <- gene_tpm[[1]]
gene_tpm <- gene_tpm[, -1, drop = FALSE]

cat("✔ Gene TPM matrix:", dim(gene_tpm), "\n")

# ============================================================
# LOAD METADATA
# ============================================================
meta <- read_csv(meta_file, show_col_types = FALSE) %>%
  as.data.frame() %>%
  filter(Sample_ID %in% rownames(gene_tpm))

rownames(meta) <- meta$Sample_ID

meta$Group <- factor(
  meta$Group,
  levels = c("Healthy", "Adenoma", "Cancer")
)

# ============================================================
# ALIGN SAMPLES
# ============================================================
common_samples <- intersect(rownames(gene_tpm), meta$Sample_ID)

gene_tpm <- gene_tpm[common_samples, , drop = FALSE]
meta <- meta[common_samples, ]

stopifnot(all(rownames(gene_tpm) == meta$Sample_ID))
cat("✔ Samples used:", length(common_samples), "\n")

# ============================================================
# LOAD PLS TOP GENES
# ============================================================
pls_genes <- read_csv(pls_file, show_col_types = FALSE) %>%
  filter(Component %in% c("PLS_1", "PLS_9")) %>%
  select(Component, Gene) %>%
  distinct()

# ============================================================
# FILTER TPM MATRIX TO PLS GENES
# ============================================================
pls_gene_list <- intersect(colnames(gene_tpm), pls_genes$Gene)

if (length(pls_gene_list) == 0) {
  stop("No PLS genes matched TPM gene names")
}

gene_tpm_pls <- gene_tpm[, pls_gene_list, drop = FALSE]
cat("✔ PLS genes used:", ncol(gene_tpm_pls), "\n")

# ============================================================
# COMPUTE GROUP MEANS
# ============================================================
group_means <- gene_tpm_pls %>%
  as.data.frame() %>%
  mutate(Group = meta$Group) %>%
  group_by(Group) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  column_to_rownames("Group")

group_means <- group_means + pseudocount

# ============================================================
# PAIRWISE LOG2 FOLD CHANGES
# ============================================================
log2fc <- tibble(
  Gene = colnames(group_means),
  log2FC_Adenoma_vs_Healthy =
    as.numeric(log2(group_means["Adenoma", ] / group_means["Healthy", ])),
  log2FC_Cancer_vs_Healthy =
    as.numeric(log2(group_means["Cancer", ] / group_means["Healthy", ])),
  log2FC_Cancer_vs_Adenoma =
    as.numeric(log2(group_means["Cancer", ] / group_means["Adenoma", ]))
)

# ============================================================
# ADD COMPONENT + GENE TYPE
# ============================================================
log2fc <- log2fc %>%
  left_join(pls_genes, by = "Gene") %>%
  mutate(
    Gene_type = ifelse(startsWith(Gene, "ARG_"), "ARG", "VF")
  )

write_csv(
  log2fc,
  "PLS1_PLS9_Genes_Pairwise_log2FC.csv"
)

cat("✔ Saved: PLS1_PLS9_Genes_Pairwise_log2FC.csv\n")

# ============================================================
# HEATMAP FUNCTION (NO CAPPING)
# ============================================================
plot_pls_heatmap <- function(component, top_n, outfile) {
  
  df_sub <- log2fc %>%
    filter(Component == component) %>%
    mutate(
      max_abs_fc = pmax(
        abs(log2FC_Adenoma_vs_Healthy),
        abs(log2FC_Cancer_vs_Healthy),
        abs(log2FC_Cancer_vs_Adenoma)
      )
    ) %>%
    slice_max(max_abs_fc, n = top_n)
  
  if (nrow(df_sub) < 2) {
    warning(paste("Skipping", component, "- not enough genes"))
    return(NULL)
  }
  
  mat <- df_sub %>%
    select(
      Gene,
      log2FC_Adenoma_vs_Healthy,
      log2FC_Cancer_vs_Healthy,
      log2FC_Cancer_vs_Adenoma
    ) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  colnames(mat) <- c(
    "Adenoma vs Healthy",
    "Cancer vs Healthy",
    "Cancer vs Adenoma"
  )
  
  row_annot <- df_sub %>%
    select(Gene, Gene_type) %>%
    column_to_rownames("Gene")
  
  ann_colors <- list(
    Gene_type = c(
      ARG = "#4C72B0",
      VF  = "#DD8452"
    )
  )
  
  pheatmap(
    mat,
    scale = "none",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_row = row_annot,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    border_color = "black",
    fontsize_row = 7,
    fontsize_col = 11,
    angle_col = 45,
    main = paste0(component, " Genes — Pairwise log2 Fold Change"),
    filename = outfile,
    width = 6.5,
    height = max(6, nrow(mat) * 0.22)
  )
}

# ============================================================
# DRAW HEATMAPS
# ============================================================
plot_pls_heatmap(
  component = "PLS_1",
  top_n = 50,
  outfile = "Heatmap_PLS1_Genes_log2FC_RAW.png"
)

plot_pls_heatmap(
  component = "PLS_9",
  top_n = 50,
  outfile = "Heatmap_PLS9_Genes_log2FC_RAW.png"
)

cat("RAW log2FC PLS₁ / PLS₉ heatmaps generated\n")
#############################################################################
# ============================================================
# Supplementary Figure:
# |log2FC| Distribution for PLS₁ vs PLS₉ Genes
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# INPUT FILE
# (from previous analysis)
# ------------------------------------------------------------
log2fc_file <- "PLS1_PLS9_Genes_Pairwise_log2FC.csv"

# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------
df <- read_csv(log2fc_file, show_col_types = FALSE)

# ------------------------------------------------------------
# LONG FORMAT + ABS(log2FC)
# ------------------------------------------------------------
df_long <- df %>%
  filter(Component %in% c("PLS_1", "PLS_9")) %>%
  pivot_longer(
    cols = starts_with("log2FC_"),
    names_to = "Comparison",
    values_to = "log2FC"
  ) %>%
  mutate(
    abs_log2FC = abs(log2FC),
    Comparison = recode(
      Comparison,
      "log2FC_Adenoma_vs_Healthy" = "Adenoma vs Healthy",
      "log2FC_Cancer_vs_Healthy"  = "Cancer vs Healthy",
      "log2FC_Cancer_vs_Adenoma"  = "Cancer vs Adenoma"
    ),
    Component = factor(Component, levels = c("PLS_1", "PLS_9"))
  )

# ------------------------------------------------------------
# BOXPLOT (PUBLICATION STYLE)
# ------------------------------------------------------------
p <- ggplot(
  df_long,
  aes(x = Component, y = abs_log2FC, fill = Component)
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.2,
    alpha = 0.35
  ) +
  facet_wrap(
    ~ Comparison,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "PLS_1" = "#4C72B0",
      "PLS_9" = "#DD8452"
    ),
    guide = "none"
  ) +
  labs(
    x = "PLS Component",
    y = expression("|log"[2]*" Fold Change|")
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1.3
    )
  )

# ------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------
ggsave(
  filename = "Supplementary_Boxplot_Abs_log2FC_PLS1_vs_PLS9.png",
  plot = p,
  width = 7,
  height = 4.5,
  dpi = 600
)

print(p)
cat("Supplementary |log2FC| boxplot saved\n")
#################################################################################

