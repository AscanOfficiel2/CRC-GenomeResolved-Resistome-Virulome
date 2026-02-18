###############################################################################
## ARG + VF CLR + Batch Correction Pipeline
## Mirrors the species pipeline but corrected for ARG/VF matrix structure.
## Includes: CLR, zero-var filter, PCA/UMAP/t-SNE before/after,
## PERMANOVA (Euclidean), LIMMA batch correction with proper diagnostics.
##
## Final publication-ready version for AbdulAziz Ascandari
###############################################################################

set.seed(43)

# --- Load libraries ---
library(limma)
library(compositions)
library(vegan)
library(ggplot2)
library(gridExtra)
library(umap)
library(Rtsne)
library(dplyr)
library(shadowtext)

###############################################################################
# THEME + COLORS
###############################################################################

pub_theme <- theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.4),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

group_colors <- c(
  "Healthy" = "#ffcc00",
  "Adenoma" = "#984ea3",
  "Cancer"  = "#4daf4a"
)

###############################################################################
## STEP 1 — TPM–METADATA ALIGNMENT (Correct for ARG+VF matrix)
###############################################################################

tpm <- read.csv("Combined_ARG_VF_TPM_Matrix.csv", row.names = 1, check.names = FALSE)
meta_raw <- read.csv("metadata_MAGs.csv", check.names = FALSE)

# Clean whitespace
rownames(tpm) <- trimws(rownames(tpm))
meta_raw$Sample_ID <- trimws(meta_raw$Sample_ID)

# Identify samples in both TPM and metadata
common_ids <- intersect(rownames(tpm), meta_raw$Sample_ID)
if (length(common_ids) == 0) stop("No overlapping sample IDs between TPM and metadata.")

# Filter TPM & align metadata
tpm_filtered <- tpm[common_ids, , drop = FALSE]
meta_clean  <- meta_raw[match(common_ids, meta_raw$Sample_ID), ]

write.csv(meta_clean, "ARG_VF_metadata_final_matched.csv", row.names = FALSE)
write.csv(tpm_filtered, "ARG_VF_TPM_raw_filtered.csv")
cat("✔ TPM–metadata alignment complete.\n")

###############################################################################
# STEP 2 — Replace NAs
###############################################################################

tpm_filtered[is.na(tpm_filtered)] <- 0
write.csv(tpm_filtered, "ARG_VF_TPM_raw_filtered.csv")

###############################################################################
# STEP 3 — CLR TRANSFORMATION
###############################################################################

pseudocount <- 1e-6
tpm_clr <- t(apply(tpm_filtered + pseudocount, 2, clr))
counts <- t(tpm_clr)   # samples = rows

write.csv(counts, "ARG_VF_TPM_clr.csv")
cat("✔ CLR transformation complete.\n")

###############################################################################
# STEP 4 — REMOVE ZERO-VARIANCE FEATURES
###############################################################################

var_ok <- apply(counts, 2, function(x) var(x, na.rm = TRUE) > 1e-12)
counts_filt <- counts[, var_ok]

write.csv(counts_filt, "ARG_VF_TPM_clr_filtered_zero_var.csv")
cat("✔ Zero-variance filtering complete.\n")

###############################################################################
# PLOTTING UTILITY (updated: no deprecated aes_string)
###############################################################################

plot_points <- function(df, xvar, yvar, title){
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(aes(fill = Group), 
               color = "black", shape = 21, size = 3.5, stroke = 0.7, alpha = 0.85) +
    scale_fill_manual(values = group_colors) +
    labs(title = title) +
    pub_theme
}

###############################################################################
# STEP 5 — PCA / UMAP / t-SNE BEFORE BATCH CORRECTION
###############################################################################

## PCA BEFORE
pca <- prcomp(counts_filt, scale. = TRUE)
pca_df <- data.frame(pca$x, Group = meta_clean$Group)
ggsave("ARG_VF_PCA_Before.png",
       plot_points(pca_df, "PC1","PC2","PCA Before Batch Correction"),
       width=5, height=4, dpi=600)

## UMAP BEFORE
umap_res <- umap(counts_filt)
umap_df <- data.frame(UMAP1=umap_res$layout[,1],
                      UMAP2=umap_res$layout[,2],
                      Group=meta_clean$Group)
ggsave("ARG_VF_UMAP_Before.png",
       plot_points(umap_df,"UMAP1","UMAP2","UMAP Before Batch Correction"),
       width=5, height=5, dpi=600)

## t-SNE BEFORE (remove duplicates)
uniq_idx <- !duplicated(counts_filt)
counts_tsne <- counts_filt[uniq_idx, ]
meta_tsne <- meta_clean[uniq_idx, ]

tsne <- Rtsne(counts_tsne, dims=2, perplexity=30)
tsne_df <- data.frame(tSNE1 = tsne$Y[,1],
                      tSNE2 = tsne$Y[,2],
                      Group = meta_tsne$Group)
ggsave("ARG_VF_tSNE_Before.png",
       plot_points(tsne_df,"tSNE1","tSNE2","t-SNE Before Batch Correction"),
       width=5, height=5, dpi=600)

###############################################################################
# STEP 5b — PERMANOVA BEFORE CORRECTION
###############################################################################
meta_perma <- meta_clean  # make a safe copy

# Continuous variables: replace NA with median
meta_perma$Age <- as.numeric(meta_clean$Age)
meta_perma$Age[is.na(meta_perma$Age)] <- median(meta_perma$Age, na.rm = TRUE)

meta_perma$BMI <- as.numeric(meta_clean$BMI)
meta_perma$BMI[is.na(meta_perma$BMI)] <- median(meta_perma$BMI, na.rm = TRUE)

# Categorical variables: replace NA with "Unknown"
fact_cols <- c("Group","Project","Continent","Country","Center_Name","Instrument","Sex")
meta_perma[fact_cols] <- lapply(meta_clean[fact_cols], function(x){
  x <- as.character(x)
  x[is.na(x)] <- "Unknown"
  as.factor(x)
})

per_before <- list(
  Group       = adonis2(counts_filt ~ Group,        data = meta_perma, method = "euclidean"),
  Project     = adonis2(counts_filt ~ Project,      data = meta_perma, method = "euclidean"),
  Continent   = adonis2(counts_filt ~ Continent,    data = meta_perma, method = "euclidean"),
  Country     = adonis2(counts_filt ~ Country,      data = meta_perma, method = "euclidean"),
  Center      = adonis2(counts_filt ~ Center_Name,  data = meta_perma, method = "euclidean"),
  Instrument  = adonis2(counts_filt ~ Instrument,   data = meta_perma, method = "euclidean"),
  Age         = adonis2(counts_filt ~ Age,          data = meta_perma, method = "euclidean"),
  Sex         = adonis2(counts_filt ~ Sex,          data = meta_perma, method = "euclidean"),
  BMI         = adonis2(counts_filt ~ BMI,          data = meta_perma, method = "euclidean")
)

per_before_df <- bind_rows(
  lapply(names(per_before), function(x){
    df <- as.data.frame(per_before[[x]])
    df$Factor <- x
    df
  })
) %>%
  filter(!is.na(F)) %>%
  group_by(Factor) %>%
  slice(1) %>%
  ungroup()

write.csv(per_before_df, "ARG_VF_PERMANOVA_before_cleaned.csv", row.names=FALSE)
cat("PERMANOVA Before Correction saved.\n")

###############################################################################
# STEP 6 — LIMMA BATCH CORRECTION (v1: covariates, v2: covariates + project)
###############################################################################

design <- model.matrix(~ Group, data=meta_clean)
covars <- model.matrix(~ Instrument + Center_Name, data=meta_clean)[,-1]

## Correction Model 1
corrected_v1 <- t(removeBatchEffect(
  t(counts_filt),
  covariates = covars,
  design = design
))
write.csv(corrected_v1, "ARG_VF_TPM_clr_batch_corrected_v1.csv")

## Correction Model 2
corrected_v2 <- t(removeBatchEffect(
  t(counts_filt),
  batch = meta_clean$Project,
  covariates = covars,
  design = design
))
write.csv(corrected_v2, "ARG_VF_TPM_clr_batch_corrected_v2.csv")

###############################################################################
# STEP 7 — DIAGNOSTICS AFTER CORRECTION
###############################################################################

diagnostics <- function(mat, meta, prefix){
  
  ## PCA AFTER
  pca <- prcomp(mat, scale.=TRUE)
  pca_df <- data.frame(pca$x, Group = meta$Group)
  ggsave(paste0("ARG_VF_PCA_after_",prefix,".png"),
         plot_points(pca_df,"PC1","PC2",paste("PCA After Correction",prefix)),
         width=5, height=4, dpi=600)
  
  ## PERMANOVA AFTER
  per <- list(
    Group     = adonis2(mat ~ Group,       data=meta, method="euclidean"),
    Instrument= adonis2(mat ~ Instrument,  data=meta, method="euclidean"),
    Center    = adonis2(mat ~ Center_Name, data=meta, method="euclidean"),
    Project   = adonis2(mat ~ Project,     data=meta, method="euclidean")
  )
  per_df <- bind_rows(lapply(names(per), function(x){
    df <- as.data.frame(per[[x]])
    df$Factor <- x
    df
  })) %>% filter(!is.na(F)) %>%
    group_by(Factor) %>% slice(1) %>% ungroup()
  
  write.csv(per_df, paste0("ARG_VF_PERMANOVA_cleaned_",prefix,".csv"), row.names=FALSE)
  
  ## UMAP AFTER
  um_res <- umap(mat)
  um_df <- data.frame(UMAP1=um_res$layout[,1],
                      UMAP2=um_res$layout[,2],
                      Group=meta$Group)
  ggsave(paste0("ARG_VF_UMAP_after_",prefix,".png"),
         plot_points(um_df,"UMAP1","UMAP2",paste("UMAP After Correction",prefix)),
         width=5, height=5, dpi=600)
  
  ## t-SNE AFTER
  uniq_idx <- !duplicated(mat)
  mat_tsne <- mat[uniq_idx, ]
  meta_tsne <- meta[uniq_idx, ]
  
  ts <- Rtsne(mat_tsne, dims=2, perplexity=30)
  ts_df <- data.frame(tSNE1=ts$Y[,1],
                      tSNE2=ts$Y[,2],
                      Group=meta_tsne$Group)
  
  ggsave(paste0("ARG_VF_tSNE_after_",prefix,".png"),
         plot_points(ts_df,"tSNE1","tSNE2",paste("t-SNE After Correction",prefix)),
         width=5, height=5, dpi=600)
}

diagnostics(corrected_v1, meta_clean, "v1")
diagnostics(corrected_v2, meta_clean, "v2")

cat(" ARG+VF CLR + Batch Correction Completed Successfully\n")
