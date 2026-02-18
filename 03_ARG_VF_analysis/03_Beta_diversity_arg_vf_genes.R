###############################################################################
# ARG + VF — BETA DIVERSITY USING CLR+BATCH CORRECTED MATRIX (v2)
###############################################################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(tibble)

set.seed(42)

###############################################################################
# STEP 1 — LOAD CLR MATRIX & METADATA
###############################################################################

cat("Loading CLR batch-corrected matrix...\n")

clr_full <- read.csv("ARG_VF_TPM_clr_batch_corrected_v2.csv",
                     row.names = 1, check.names = FALSE)

meta <- read.csv("ARG_VF_metadata_final_matched.csv",
                 row.names = 1, check.names = FALSE)

# Align samples
common <- intersect(rownames(clr_full), rownames(meta))
clr_full <- clr_full[common, , drop=FALSE]
meta <- meta[common, , drop=FALSE]

cat("Samples aligned:", nrow(meta), "samples\n")

###############################################################################
# REMOVE SAMPLES WITH MISSING METADATA
###############################################################################

cat("Filtering samples with missing (Group, Age, BMI, Sex)...\n")

meta <- meta %>% drop_na(Group, Age, BMI, Sex)
clr_full <- clr_full[rownames(meta), , drop=FALSE]

cat("Remaining samples:", nrow(meta), "\n\n")

###############################################################################
# STEP 2 — SPLIT INTO ARG and VF CLR MATRICES + FILTER FEATURES
###############################################################################

arg_cols <- grep("^ARG_", colnames(clr_full), value=TRUE)
vf_cols  <- grep("^VF_",  colnames(clr_full), value=TRUE)

clr_ARG <- clr_full[, arg_cols, drop=FALSE]
clr_VF  <- clr_full[, vf_cols, drop=FALSE]

# Filter features present in ≥5% of samples
min_freq <- 0.05
n <- nrow(clr_full)

arg_keep <- names(which(colSums(clr_ARG != 0)/n >= min_freq))
vf_keep  <- names(which(colSums(clr_VF != 0)/n >= min_freq))

clr_ARG <- clr_ARG[, arg_keep, drop=FALSE]
clr_VF  <- clr_VF[, vf_keep, drop=FALSE]

cat("ARG CLR features kept:", length(arg_keep), "\n")
cat("VF CLR features kept :", length(vf_keep), "\n\n")

###############################################################################
# STEP 3 — COMPUTE AITCHISON DISTANCE
###############################################################################

dist_ARG <- dist(clr_ARG, method="euclidean")
dist_VF  <- dist(clr_VF,  method="euclidean")

###############################################################################
# STEP 4 — NMDS
###############################################################################
library(parallel)


NMDS_ARG <- metaMDS(clr_ARG, distance="euclidean", k=2, trymax=50, parallel = detectCores() - 1)
NMDS_VF  <- metaMDS(clr_VF,  distance="euclidean", k=2, trymax=50, parallel = detectCores() - 1)

# Save stress values
write.csv(data.frame(Stress = NMDS_ARG$stress),
          "ARG_NMDS_Stress.csv", row.names=FALSE)
write.csv(data.frame(Stress = NMDS_VF$stress),
          "VF_NMDS_Stress.csv", row.names=FALSE)

###############################################################################
# Fix metadata for joining (Sample_ID must be a column)
###############################################################################
meta <- meta %>% rownames_to_column("Sample_ID")

###############################################################################
# STEP 5 — NMDS PLOTS
###############################################################################

COLORS <- c("Healthy"="#ffcc00",
            "Adenoma"="#984ea3",
            "Cancer" ="#4daf4a")

plot_NMDS <- function(N, meta, prefix){
  
  df <- as.data.frame(N$points) %>%
    rownames_to_column("Sample_ID") %>%
    left_join(meta, by="Sample_ID")
  
  p <- ggplot(df, aes(MDS1, MDS2, fill=Group)) +
    geom_point(shape=21, size=3, alpha=0.9, color="black") +
    scale_fill_manual(values=COLORS) +
    theme_classic(base_size=13) +
    theme(panel.border = element_rect(color="black", fill=NA, linewidth=1.2),
          legend.position="bottom") +
    labs(title=paste(prefix,"NMDS (Aitchison)"))
  
  ggsave(paste0(prefix,"_NMDS.png"), p, width=5, height=4.8, dpi=600)
}

plot_NMDS(NMDS_ARG, meta, "ARG")
plot_NMDS(NMDS_VF,  meta, "VF")

###############################################################################
# STEP 6 — PERMANOVA
###############################################################################

perm_ARG <- adonis2(dist_ARG ~ Group + Age + BMI + Sex, data=meta, permutations=999, parallel = detectCores() - 1)
perm_VF  <- adonis2(dist_VF  ~ Group + Age + BMI + Sex, data=meta, permutations=999, parallel = detectCores() - 1)

write.csv(as.data.frame(perm_ARG), "ARG_PERMANOVA.csv")
write.csv(as.data.frame(perm_VF),  "VF_PERMANOVA.csv")

###############################################################################
# STEP 7 — PERMDISP
###############################################################################

disp_ARG <- betadisper(dist_ARG, meta$Group)
disp_VF  <- betadisper(dist_VF,  meta$Group)

write.csv(as.data.frame(anova(disp_ARG)), "ARG_PERMDISP.csv")
write.csv(as.data.frame(anova(disp_VF)),  "VF_PERMDISP.csv")

###############################################################################
# STEP 8 — ENVFIT
###############################################################################

env_data <- meta %>% select(Age, BMI, Sex)

fit_ARG <- envfit(NMDS_ARG, env_data, permutations=999, parallel = detectCores() - 1)
fit_VF  <- envfit(NMDS_VF,  env_data, permutations=999, parallel = detectCores() - 1)

sink("ARG_ENVFIT.txt"); print(fit_ARG); sink()
sink("VF_ENVFIT.txt");  print(fit_VF);  sink()

###############################################################################
# STEP 9 — Pairwise PERMANOVA
###############################################################################

pairwise_permanova <- function(distmat, group, prefix){
  lev <- unique(group)
  pairs <- combn(lev, 2, simplify=FALSE)
  out <- list()
  
  for (p in pairs){
    idx <- group %in% p
    dsub <- as.dist(as.matrix(distmat)[idx, idx])
    metasub <- data.frame(Group = factor(group[idx]))
    ad <- adonis2(dsub ~ Group, data=metasub, permutations=999, parallel = detectCores() - 1)
    
    out[[paste0(p[1],"_vs_",p[2])]] <- tibble(
      Comparison = paste(p, collapse=" vs "),
      R2 = ad$R2[1],
      F  = ad$F[1],
      p  = ad$`Pr(>F)`[1]
    )
  }
  
  out_df <- bind_rows(out)
  write.csv(out_df, paste0(prefix,"_Pairwise_PERMANOVA.csv"), row.names=FALSE)
}

pairwise_permanova(dist_ARG, meta$Group, "ARG")
pairwise_permanova(dist_VF,  meta$Group, "VF")

###############################################################################
# STEP 11 — ARG–VF Resistome Similarity (Mantel + Procrustes)
###############################################################################

cat("🔗 Computing ARG–VF resistome similarity…\n")

mantel_res <- mantel(dist_ARG, dist_VF, method="spearman", permutations=999)
proc <- protest(NMDS_ARG, NMDS_VF, permutations=999)

similarity_df <- tibble(
  Mantel_r = mantel_res$statistic,
  Mantel_p = mantel_res$signif,
  Procrustes_r = proc$t0,
  Procrustes_p = proc$signif
)

write.csv(similarity_df, "ARG_VF_ResistomeSimilarity.csv", row.names=FALSE)

cat("\n BETA diversity + BIOENV + ARG–VF similarity COMPLETED.\n")

###############################################################################
# STEP 12 — PLOT: Procrustes Alignment (ARG → VF)
###############################################################################

cat(" Creating Procrustes alignment plot…\n")

# Extract Procrustes coordinates
proc_fit <- procrustes(NMDS_ARG, NMDS_VF)

# Extract coordinates into a tidy dataframe
proc_df <- data.frame(
  Sample_ID = rownames(proc_fit$Yrot),
  ARG_MDS1  = proc_fit$X[,1],
  ARG_MDS2  = proc_fit$X[,2],
  VF_MDS1   = proc_fit$Yrot[,1],
  VF_MDS2   = proc_fit$Yrot[,2]
)

# Add metadata
proc_df <- proc_df %>% 
  left_join(meta, by = "Sample_ID")

# Colors
COLORS <- c(
  "Healthy" = "#ffcc00",
  "Adenoma" = "#984ea3",
  "Cancer"  = "#4daf4a"
)

# Plot
p_proc <- ggplot(proc_df) +
  geom_segment(aes(x = ARG_MDS1, y = ARG_MDS2,
                   xend = VF_MDS1, yend = VF_MDS2,
                   color = Group),
               arrow = arrow(length = unit(0.15, "cm")),
               alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(ARG_MDS1, ARG_MDS2, fill = Group),
             shape = 21, size = 3, color = "black") +
  geom_point(aes(VF_MDS1, VF_MDS2, fill = Group),
             shape = 24, size = 3, color = "black") +
  scale_fill_manual(values = COLORS) +
  scale_color_manual(values = COLORS) +
  labs(title = "Procrustes Alignment: ARG → VF",
       subtitle = paste0("Mantel r = ", round(mantel_res$statistic, 3),
                         " (p = ", mantel_res$signif, ");  ",
                         "Procrustes r = ", round(proc$t0, 3),
                         " (p = ", proc$signif, ")"),
       x = "ARG NMDS Axis 1",
       y = "ARG NMDS Axis 2") +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Save figure
ggsave("ARG_VF_ProcrustesAlignment.png",
       p_proc, width = 6, height = 5, dpi = 600)

cat("Procrustes plot saved: ARG_VF_ProcrustesAlignment.png\n")
###############################################################################
# STEP 12 — Stylized Procrustes Alignment Plot (Upgraded Publication-Ready)
###############################################################################

cat("Creating improved Procrustes alignment plot…\n")

# Compute Procrustes fit
proc_fit <- procrustes(NMDS_ARG, NMDS_VF)

# Extract coordinates
proc_df <- data.frame(
  Sample_ID = rownames(proc_fit$Yrot),
  ARG_MDS1  = proc_fit$X[,1],
  ARG_MDS2  = proc_fit$X[,2],
  VF_MDS1   = proc_fit$Yrot[,1],
  VF_MDS2   = proc_fit$Yrot[,2]
)

# Residual lengths = how far VF points must move to align to ARG
proc_df$Residual <- sqrt((proc_fit$X[,1] - proc_fit$Yrot[,1])^2 +
                           (proc_fit$X[,2] - proc_fit$Yrot[,2])^2)

# Attach metadata
proc_df <- proc_df %>% left_join(meta, by="Sample_ID")

# Colors
COLORS <- c("Healthy"="#ffcc00",
            "Adenoma"="#984ea3",
            "Cancer" ="#4daf4a")

# Optional hull function — helps visualize group separation
library(ggforce)
compute_hull <- function(df) df[chull(df$ARG_MDS1, df$ARG_MDS2), ]

hulls <- proc_df %>%
  group_by(Group) %>%
  do(compute_hull(.)) %>%
  ungroup()

# Plot
p_proc <- ggplot(proc_df) +
  
  # Group convex hulls (optional but very helpful)
  geom_polygon(data=hulls, aes(ARG_MDS1, ARG_MDS2, fill=Group),
               alpha=0.12, color=NA) +
  
  # Segments colored by residual distance
  geom_segment(aes(x=ARG_MDS1, y=ARG_MDS2,
                   xend=VF_MDS1, yend=VF_MDS2,
                   color=Residual),
               alpha=0.55, linewidth=0.55) +
  
  # ARG positions (circles)
  geom_point(aes(ARG_MDS1, ARG_MDS2, fill=Group),
             shape=21, size=3, color="black", alpha=0.95) +
  
  # VF positions (triangles)
  geom_point(aes(VF_MDS1, VF_MDS2, fill=Group),
             shape=24, size=3, color="black", alpha=0.95) +
  
  # Color scales
  scale_fill_manual(values=COLORS) +
  scale_color_gradient(low="grey70", high="darkred",
                       name="Residual Distance") +
  
  # Title + subtitle
  labs(
    title = "Procrustes Alignment: ARG → VF",
    subtitle = paste0(
      "Mantel r = ", round(mantel_res$statistic, 3),
      " (p = ", mantel_res$signif, ");  ",
      "Procrustes r = ", round(proc$t0, 3),
      " (p = ", proc$signif, ")"
    ),
    x = "ARG NMDS Axis 1",
    y = "ARG NMDS Axis 2"
  ) +
  
  # Clean theme
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1.2),
    legend.position = "bottom",
    legend.title = element_text(size=10),
    legend.key.width = unit(1.4, "cm")
  ) +
  
  coord_equal()   # Prevent distortion

# Save figure
ggsave("ARG_VF_ProcrustesAlignment11.png",
       p_proc, width=7, height=6, dpi=600)

cat("Improved Procrustes plot saved: ARG_VF_ProcrustesAlignment.png\n")
###############################################################################################
###############################################################################
### GENE-LEVEL BETA DIVERSITY (ARG + VF)
### Using ARG_VF_TPM_clr_batch_corrected_v2.csv
###############################################################################

library(tidyverse)
library(vegan)
library(ggplot2)
library(stringr)

###############################################################################
### 1. LOAD CLR GENE MATRIX & METADATA
###############################################################################

# Load gene-level CLR matrix
clr_gene <- read.csv("ARG_VF_TPM_clr_batch_corrected_v2.csv",
                     row.names = 1, check.names = FALSE)

# Load metadata
meta <- read.csv("ARG_VF_metadata_final_matched.csv")

# Clean metadata
meta <- meta %>% drop_na(Sample_ID, Group, Age, BMI, Sex)
meta$Group <- factor(meta$Group, levels = c("Healthy","Adenoma","Cancer"))

# Align samples
common <- intersect(rownames(clr_gene), meta$Sample_ID)
clr_gene <- clr_gene[common, , drop=FALSE]
meta <- meta[match(common, meta$Sample_ID), ]

cat("Gene-level samples aligned:", nrow(meta), "samples\n")

###############################################################################
### 2. SPLIT INTO ARG GENES AND VF GENES
###############################################################################

arg_gene_cols <- grep("^ARG_", colnames(clr_gene), value=TRUE)
vf_gene_cols  <- grep("^VF_",  colnames(clr_gene), value=TRUE)

clr_ARG_gene <- clr_gene[, arg_gene_cols, drop=FALSE]
clr_VF_gene  <- clr_gene[, vf_gene_cols,  drop=FALSE]

cat("ARG genes:", ncol(clr_ARG_gene), "\n")
cat("VF genes :", ncol(clr_VF_gene), "\n")

###############################################################################
### 3. COMPUTE AITCHISON DISTANCE (Euclidean in CLR space)
###############################################################################

dist_ARG_gene <- dist(clr_ARG_gene, method="euclidean")
dist_VF_gene  <- dist(clr_VF_gene,  method="euclidean")

###############################################################################
### 4. PCoA — ARG GENES
###############################################################################

pcoa_ARG_gene <- cmdscale(dist_ARG_gene, eig=TRUE, k=3)
var_ARG_gene <- round(100 * pcoa_ARG_gene$eig[1:2] / sum(pcoa_ARG_gene$eig), 2)

df_pcoa_ARG <- data.frame(
  PCoA1 = pcoa_ARG_gene$points[,1],
  PCoA2 = pcoa_ARG_gene$points[,2],
  Group = meta$Group
)

p_ARG_gene <- ggplot(df_pcoa_ARG, aes(PCoA1, PCoA2, fill = Group)) +
  geom_point(shape=21, size=3.5, color="black", alpha=0.85) +
  scale_fill_manual(values = c("Healthy"="#ffcc00","Adenoma"="#984ea3","Cancer"="#4daf4a")) +
  coord_equal() +
  theme_classic(base_size=14) +
  labs(
    title = "ARG Gene-Level PCoA (Aitchison)",
    x = paste0("PCoA1 (", var_ARG_gene[1], "%)"),
    y = paste0("PCoA2 (", var_ARG_gene[2], "%)")
  )

ggsave("ARG_Genes_PCoA2D.png", p_ARG_gene, width=5, height=5, dpi=600)

###############################################################################
### 5. PCoA — VF GENES
###############################################################################

pcoa_VF_gene <- cmdscale(dist_VF_gene, eig=TRUE, k=3)
var_VF_gene <- round(100 * pcoa_VF_gene$eig[1:2] / sum(pcoa_VF_gene$eig), 2)

df_pcoa_VF <- data.frame(
  PCoA1 = pcoa_VF_gene$points[,1],
  PCoA2 = pcoa_VF_gene$points[,2],
  Group = meta$Group
)

p_VF_gene <- ggplot(df_pcoa_VF, aes(PCoA1, PCoA2, fill = Group)) +
  geom_point(shape=21, size=3.5, color="black", alpha=0.85) +
  scale_fill_manual(values = c("Healthy"="#ffcc00","Adenoma"="#984ea3","Cancer"="#4daf4a")) +
  coord_equal() +
  theme_classic(base_size=14) +
  labs(
    title = "VF Gene-Level PCoA (Aitchison)",
    x = paste0("PCoA1 (", var_VF_gene[1], "%)"),
    y = paste0("PCoA2 (", var_VF_gene[2], "%)")
  )

ggsave("VF_Genes_PCoA2D.png", p_VF_gene, width=5, height=5, dpi=600)

###############################################################################
### 6. CAP ANALYSIS — ARG GENES
###############################################################################

cap_ARG_gene <- capscale(dist_ARG_gene ~ Group, data = meta, distance="euclidean")
cap_ARG_gene_test <- anova(cap_ARG_gene, permutations = 999)

arg_scores <- scores(cap_ARG_gene, display="sites")

df_CAP_ARG <- data.frame(
  CAP1 = arg_scores[,1],
  CAP2 = arg_scores[,2],
  Group = meta$Group
)

eig_ARG_gene <- cap_ARG_gene$CCA$eig
var_CAP1_ARG <- round(100 * eig_ARG_gene[1] / sum(eig_ARG_gene), 2)
var_CAP2_ARG <- round(100 * eig_ARG_gene[2] / sum(eig_ARG_gene), 2)

p_CAP_ARG <- ggplot(df_CAP_ARG, aes(CAP1, CAP2, fill = Group)) +
  geom_point(shape=21, size=3.5, color="black") +
  scale_fill_manual(values = c("Healthy"="#ffcc00","Adenoma"="#984ea3","Cancer"="#4daf4a")) +
  coord_equal() +
  theme_classic(base_size=14) +
  labs(
    title="ARG Gene-Level CAP (Group-Constrained)",
    x = paste0("CAP1 (", var_CAP1_ARG, "%)"),
    y = paste0("CAP2 (", var_CAP2_ARG, "%)")
  )

ggsave("ARG_Genes_CAP.png", p_CAP_ARG, width=5, height=5, dpi=600)

write.csv(as.data.frame(cap_ARG_gene_test), "ARG_Genes_CAP_ANOVA.csv", row.names=FALSE)
write.csv(df_CAP_ARG, "ARG_Genes_CAP_SiteScores.csv", row.names=FALSE)

###############################################################################
### 7. CAP ANALYSIS — VF GENES
###############################################################################

cap_VF_gene <- capscale(dist_VF_gene ~ Group, data = meta, distance="euclidean")
cap_VF_gene_test <- anova(cap_VF_gene, permutations = 999)

vf_scores <- scores(cap_VF_gene, display="sites")

df_CAP_VF <- data.frame(
  CAP1 = vf_scores[,1],
  CAP2 = vf_scores[,2],
  Group = meta$Group
)

eig_VF_gene <- cap_VF_gene$CCA$eig
var_CAP1_VF <- round(100 * eig_VF_gene[1] / sum(eig_VF_gene), 2)
var_CAP2_VF <- round(100 * eig_VF_gene[2] / sum(eig_VF_gene), 2)

p_CAP_VF <- ggplot(df_CAP_VF, aes(CAP1, CAP2, fill = Group)) +
  geom_point(shape=21, size=3.5, color="black") +
  scale_fill_manual(values = c("Healthy"="#ffcc00","Adenoma"="#984ea3","Cancer"="#4daf4a")) +
  coord_equal() +
  theme_classic(base_size=14) +
  labs(
    title="VF Gene-Level CAP (Group-Constrained)",
    x = paste0("CAP1 (", var_CAP1_VF, "%)"),
    y = paste0("CAP2 (", var_CAP2_VF, "%)")
  )

ggsave("VF_Genes_CAP.png", p_CAP_VF, width=5, height=5, dpi=600)

write.csv(as.data.frame(cap_VF_gene_test), "VF_Genes_CAP_ANOVA.csv", row.names=FALSE)
write.csv(df_CAP_VF, "VF_Genes_CAP_SiteScores.csv", row.names=FALSE)

cat("\n Gene-level PCoA + CAP analysis for ARG & VF COMPLETED!\n")
#####################################################################################
###############################################################################
### SAVE ARG GENE-LEVEL CAP SUMMARY
###############################################################################

sink("ARG_Genes_CAP_Summary.txt")
print(summary(cap_ARG_gene))
sink()

###############################################################################
### SAVE VF GENE-LEVEL CAP SUMMARY
###############################################################################

sink("VF_Genes_CAP_Summary.txt")
print(summary(cap_VF_gene))
sink()

############################################