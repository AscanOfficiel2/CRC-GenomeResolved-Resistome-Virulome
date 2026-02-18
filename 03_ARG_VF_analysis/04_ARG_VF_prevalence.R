###############################################################################
# ARG FUNCTIONAL PROFILING — FINAL FIXED VERSION
###############################################################################
library(tidyverse)
library(ggplot2)
library(pheatmap)

set.seed(42)

###############################################################################
# STEP 1 — Load data
###############################################################################
arg <- read.csv("arg_cleaned_final_with_mechanisms.csv", check.names = FALSE)
meta <- read.csv("ARG_VF_metadata_final_matched.csv", check.names = FALSE)

meta$Sample_ID <- trimws(meta$Sample_ID)

cat("Loaded ARG annotation rows:", nrow(arg), "\n\n")

###############################################################################
# STEP 2 — Detect actual mechanism column name
###############################################################################

mech_col <- grep("mechanism", colnames(arg), ignore.case = TRUE, value = TRUE)

cat("Detected mechanism column(s):\n")
print(mech_col)

if(length(mech_col) != 1){
  stop("Could not uniquely identify mechanism column. Found: ",
       paste(mech_col, collapse=", "))
}

###############################################################################
# STEP 3 — Prepare long table
###############################################################################

arg_long <- arg %>%
  transmute(
    Sample    = Sample_ID,
    TPM       = TPM,
    Mechanism = .data[[mech_col]]
  ) %>%
  filter(!is.na(Mechanism)) %>%
  left_join(meta, by = c("Sample" = "Sample_ID"))

# >>> FIX GROUP ORDER <<<
arg_long$Group <- factor(arg_long$Group,
                         levels = c("Healthy", "Adenoma", "Cancer"))

cat("\nColumns in arg_long:\n")
print(colnames(arg_long))


###############################################################################
# STEP 4 — Compute mechanism mean TPM per group
###############################################################################

arg_mech_group <- arg_long %>%
  group_by(Group, Mechanism) %>%
  summarise(MeanTPM = mean(TPM), .groups = "drop")


cat("\nColumns in arg_mech_group:\n")
print(colnames(arg_mech_group))

###############################################################################
# STEP 5 — Top 10 mechanisms
###############################################################################

top10_mech <- arg_mech_group %>%
  group_by(Mechanism) %>%
  summarise(Total = sum(MeanTPM)) %>%
  arrange(desc(Total)) %>%
  slice(1:10) %>%
  pull(Mechanism)

###############################################################################
# STEP 6 — Prepare plot table
###############################################################################

arg_top10 <- arg_mech_group %>%
  mutate(
    Mechanism = ifelse(Mechanism %in% top10_mech, Mechanism, "Other")
  )

###############################################################################
# STEP 7 — Barplot
###############################################################################

p_major <- ggplot(arg_top10, aes(Group, MeanTPM, fill = Mechanism)) +
  geom_col(color="black", width = 0.7) +
  scale_fill_brewer(palette="Set3") +
  theme_classic(base_size=15) +
  theme(
    panel.border = element_rect(color="black", fill=NA, size=1.2),
    axis.text.x = element_text(angle=45, hjust=1, size=13)
  ) +
  labs(
    title = "Top ARG Mechanisms Across CRC Groups",
    y = "Mean TPM", x = ""
  )

ggsave("ARG_TopMechanisms_Group.png", p_major, width=7, height=5, dpi=600)

cat("\n✔ Saved: ARG_TopMechanisms_Group.png\n")
cat("ARG Functional Mechanism Profiling Completed Successfully!\n")
###########################################################################################
#################################################################################
# ARG ANTIBIOTIC CLASS PROFILING → WHO AWaRe CATEGORY (Final Master Script)
###############################################################################
library(tidyverse)
library(ggplot2)

set.seed(42)

###############################################################################
# 1. LOAD DATA
###############################################################################
arg  <- read.csv("arg_cleaned_final_with_mechanisms.csv", check.names = FALSE)
meta <- read.csv("ARG_VF_metadata_final_matched.csv", check.names = FALSE)
aware <- read.csv("Aware.csv", check.names = FALSE)

arg$Sample_ID  <- trimws(arg$Sample_ID)
meta$Sample_ID <- trimws(meta$Sample_ID)

cat("✔ ARG rows:", nrow(arg), "\n")
cat("✔ AWARE entries:", nrow(aware), "\n")

###############################################################################
# 2. CLEAN AWaRe TABLE
###############################################################################
aware_clean <- aware %>%
  transmute(
    Class = trimws(Class),
    Category = trimws(Category)   # Access / Watch / Reserve
  ) %>%
  distinct()

###############################################################################
# 3. PREPARE ARG TABLE (split multi-class entries)
###############################################################################
arg_class <- arg %>%
  select(Sample_ID, TPM, ARG_Drug_Class) %>%
  filter(!is.na(ARG_Drug_Class)) %>%
  separate_rows(ARG_Drug_Class, sep = ";") %>%
  mutate(ARG_Drug_Class = trimws(ARG_Drug_Class))

###############################################################################
# 4. SHORTEN DRUG CLASS NAMES FOR CONSISTENT MAPPING
###############################################################################
shorten <- c(
  "penicillin beta-lactam" = "Penicillins",
  "macrolide antibiotic" = "Macrolides",
  "lincosamide antibiotic" = "Lincosamides",
  "aminoglycoside antibiotic" = "Aminoglycosides",
  "fluoroquinolone antibiotic" = "Fluoroquinolones",
  "carbapenem" = "Carbapenems",
  "cephalosporin" = "Cephalosporins",
  "tetracycline antibiotic" = "Tetracyclines",
  "glycopeptide antibiotic" = "Glycopeptides",
  "phenicol antibiotic" = "Phenicols",
  "rifamycin antibiotic" = "Rifamycins",
  "phosphonic acid antibiotic" = "Fosfomycin",
  "diaminopyrimidine antibiotic" = "Trimethoprim",
  "sulfonamide antibiotic" = "Sulfonamides",
  "peptide antibiotic" = "Peptides",
  "monobactam" = "Monobactams"
)

arg_class$Drug_Class_Short <- plyr::mapvalues(
  arg_class$ARG_Drug_Class,
  from = names(shorten),
  to   = shorten,
  warn_missing = FALSE
)

###############################################################################
# 5. MAP → AWARE CATEGORY
###############################################################################

arg_mapped <- arg_class %>%
  left_join(aware_clean, by = c("Drug_Class_Short" = "Class"))

# Anything not matched = "Unclassified"
arg_mapped$Category[is.na(arg_mapped$Category)] <- "Unclassified"

###############################################################################
# 6. MERGE SAMPLE METADATA + FIX GROUP ORDER
###############################################################################
arg_mapped <- arg_mapped %>%
  left_join(meta, by = "Sample_ID")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# FIXED ORDER: Healthy → Adenoma → Cancer
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
arg_mapped$Group <- factor(arg_mapped$Group,
                           levels = c("Healthy", "Adenoma", "Cancer"))

###############################################################################
# 7. TPM SUM PER GROUP × AWaRe CATEGORY
###############################################################################
aware_group <- arg_mapped %>%
  group_by(Group, Category) %>%
  summarise(TotalTPM = sum(TPM), .groups = "drop")

write.csv(aware_group, "ARG_AWARE_GroupMeans.csv", row.names = FALSE)
cat("✔ Saved: ARG_AWARE_GroupMeans.csv\n")

###############################################################################
# 8. BARPLOT — WHO AWaRe CATEGORIES
###############################################################################

col_aware <- c(
  "Access"       = "#1f78b4",
  "Watch"        = "#e31a1c",
  "Reserve"      = "#33a02c",
  "Unclassified" = "grey60"
)

p_aware <- ggplot(aware_group, aes(Group, TotalTPM, fill = Category)) +
  geom_col(color="black", width=0.75) +
  scale_fill_manual(values = col_aware) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1.3),
    axis.text.x = element_text(angle=45, hjust=1, size=14),
    legend.title = element_text(face="bold", size=13)
  ) +
  labs(
    title = "ARG Antibiotic AWaRe Classification Across CRC Groups",
    y = "Total TPM",
    x = "",
    fill = "AWaRe Category"
  )

ggsave("ARG_AWARE_Barplot.png", p_aware, width = 7, height = 5, dpi = 600)
cat("✔ Saved: ARG_AWARE_Barplot.png\n")

###############################################################################
cat("AWaRe Antibiotic Class Profiling Completed Successfully!\n")
###############################################################################
###############################################################################
# ARG FUNCTIONAL PROFILING — Drug Class Level
####################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

set.seed(42)

###############################################################################
# STEP 1 — Load ARG annotation file + metadata
###############################################################################

arg <- read.csv("arg_cleaned_final_with_mechanisms.csv", check.names = FALSE)
meta <- read.csv("ARG_VF_metadata_final_matched.csv", check.names = FALSE)

# Clean sample IDs
arg$Sample_ID <- trimws(arg$Sample_ID)
meta$Sample_ID <- trimws(meta$Sample_ID)

###############################################################################
# STEP 2 — Keep relevant columns
###############################################################################

arg_small <- arg %>%
  select(Sample_ID, TPM, ARG_Drug_Class) %>%
  filter(!is.na(ARG_Drug_Class))

###############################################################################
# Important: split multiple classes
###############################################################################

arg_small <- arg_small %>%
  separate_rows(ARG_Drug_Class, sep = ";") %>%
  mutate(ARG_Drug_Class = trimws(ARG_Drug_Class))

###############################################################################
# STEP 3 — Merge metadata + FIX GROUP ORDER
###############################################################################

arg_small <- arg_small %>%
  left_join(meta, by = c("Sample_ID" = "Sample_ID")) %>%
  filter(!is.na(Group)) %>%
  mutate(Group = factor(Group, levels = c("Healthy", "Adenoma", "Cancer")))  ## ★ FIXED ORDER

###############################################################################
# STEP 4 — Compute mean TPM per Drug Class per Group
###############################################################################

drug_group <- arg_small %>%
  group_by(Group, ARG_Drug_Class) %>%
  summarise(MeanTPM = sum(TPM), .groups = "drop")

###############################################################################
# STEP 5 — Identify Top 15 Drug Classes
###############################################################################

top_drugs <- drug_group %>%
  group_by(ARG_Drug_Class) %>%
  summarise(Total = sum(MeanTPM)) %>%
  slice_max(Total, n = 15) %>%
  pull(ARG_Drug_Class)

drug_group2 <- drug_group %>%
  mutate(Class_Collapsed = ifelse(ARG_Drug_Class %in% top_drugs,
                                  ARG_Drug_Class, "Other"))

###############################################################################
# STEP 6 — Final aggregated dataframe
###############################################################################

drug_plot_df <- drug_group2 %>%
  group_by(Group, Class_Collapsed) %>%
  summarise(MeanTPM = sum(MeanTPM), .groups = "drop")

write.csv(drug_plot_df,
          "ARG_DrugClass_GroupMeans.csv",
          row.names = FALSE)

###############################################################################
# STEP 7 — Barplot
###############################################################################

palette_drugs <- c(
  RColorBrewer::brewer.pal(12, "Paired"),
  RColorBrewer::brewer.pal(8, "Set2")
)

drug_colors <- setNames(palette_drugs[1:length(unique(drug_plot_df$Class_Collapsed))],
                        unique(drug_plot_df$Class_Collapsed))

p <- ggplot(drug_plot_df, aes(x = Group, y = MeanTPM, fill = Class_Collapsed)) +
  geom_col(position = "stack", color = "black", width = 0.7) +
  scale_fill_manual(values = drug_colors) +
  labs(
    title = "ARG Drug Class Composition Across CRC Groups",
    y = "Total TPM",
    x = "",
    fill = "Drug Class"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=1.2),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave("ARG_DrugClass_Barplot.png",
       p, width = 9, height = 6, dpi = 600)

cat("ARG Drug Class Barplot saved: ARG_DrugClass_Barplot.png\n")


###################################################################################

###############################################################################
# VF FUNCTIONAL PROFILING — FINAL PUBLICATION VERSION
# Distinct muted colors for 11 mechanisms
##############################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(scales)

set.seed(42)

# --------------------------------------------------
# LOAD VF DATA (TSV inside .csv)
# --------------------------------------------------
vf <- read.delim(
  "vfdb_with_species_cleaned.csv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("VF file loaded successfully\n")
print(colnames(vf))
cat("VF rows:", nrow(vf), "\n\n")

# --------------------------------------------------
# LOAD METADATA
# --------------------------------------------------
meta <- read.csv(
  "ARG_VF_metadata_final_matched.csv",
  stringsAsFactors = FALSE
)

meta$Sample_ID <- trimws(meta$Sample_ID)

# --------------------------------------------------
# LONG FORMAT + MERGE GROUP
# --------------------------------------------------
vf_long <- vf %>%
  select(Sample_ID, TPM, Mechanism) %>%
  filter(!is.na(Mechanism)) %>%
  left_join(meta, by = "Sample_ID") %>%
  filter(!is.na(Group))

cat("VF rows after metadata join:", nrow(vf_long), "\n\n")

# --------------------------------------------------
# MEAN TPM PER GROUP × MECHANISM
# --------------------------------------------------
vf_mech_group <- vf_long %>%
  group_by(Group, Mechanism) %>%
  summarise(
    MeanTPM = mean(TPM, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------------------------------
# ENFORCE BIOLOGICAL GROUP ORDER
# --------------------------------------------------
vf_mech_group$Group <- factor(
  vf_mech_group$Group,
  levels = c("Healthy", "Adenoma", "Cancer")
)

# --------------------------------------------------
# TOP 10 MECHANISMS (GLOBAL)
# --------------------------------------------------
top10_mech <- vf_mech_group %>%
  group_by(Mechanism) %>%
  summarise(Total = sum(MeanTPM)) %>%
  slice_max(Total, n = 10) %>%
  pull(Mechanism)

vf_top10 <- vf_mech_group %>%
  mutate(
    Mechanism = ifelse(Mechanism %in% top10_mech, Mechanism, "Other")
  )

# --------------------------------------------------
# DISTINCT, MUTED, PUBLICATION-SAFE COLORS (11)
# --------------------------------------------------
mech_colors <- c(
  "Adherence"                    = "#0072B2",  # blue
  "Biofilm"                      = "#56B4E9",  # light blue
  "Effector delivery system"     = "#009E73",  # green
  "Exoenzyme"                    = "#E69F00",  # orange
  "Exotoxin"                     = "#D55E00",  # vermillion
  "Immune modulation"            = "#CC79A7",  # purple
  "Invasion"                     = "#F0E442",  # yellow
  "Motility"                     = "#999999",  # grey
  "Nutritional Metabolic factor" = "#8DD3C7",  # teal
  "Stress survival"              = "#5E5A8A",  # lavender
  "Other"                        = "#C9A227"   # muted gold (replaces grey)
)


# --------------------------------------------------
# PLOT: STACKED BARPLOT
# --------------------------------------------------
p_vf <- ggplot(
  vf_top10,
  aes(
    x = Group,
    y = MeanTPM,
    fill = Mechanism
  )
) +
  geom_col(
    width = 0.75,
    color = "black",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = mech_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = comma_format()
  ) +
  labs(
    title = "Top VF Mechanisms Across CRC Groups",
    subtitle = "Mean TPM per group",
    y = "Mean TPM",
    x = ""
  ) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1.2
    ),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )

# --------------------------------------------------
# SAVE FIGURE
# --------------------------------------------------
ggsave(
  "VF_TopMechanisms_Group.png",
  plot = p_vf,
  width = 6,
  height = 5,
  dpi = 600,
  bg = "white"
)

print(p_vf)
cat("Saved: VF_TopMechanisms_Group.png (distinct muted colors)\n")

#############################################################################################
################################################################################
# Top ARG + VF Heatmap (Species × Genes)
# Robust parsing + publication-ready ComplexHeatmap
################################################################################

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

set.seed(42)

# ============================================================
# 1. LOAD FILES (ROBUST)
# ============================================================

# ARG file (true CSV)
arg <- read.csv(
  "arg_cleaned_final_with_mechanisms.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# VF file (TSV inside .csv → MUST use read.delim)
vf <- read.delim(
  "vfdb_with_species_cleaned.csv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Metadata
meta <- read.csv(
  "ARG_VF_metadata_final_matched.csv",
  stringsAsFactors = FALSE
)

meta$Sample_ID <- trimws(meta$Sample_ID)

cat("Files loaded\n")
cat("ARG rows:", nrow(arg), "\n")
cat("VF rows :", nrow(vf), "\n\n")

# ============================================================
# 2. STANDARDIZE & COMBINE ARG + VF
# ============================================================

arg2 <- arg %>%
  select(Sample_ID, Gene, TPM, Species) %>%
  mutate(Type = "ARG")

vf2 <- vf %>%
  select(Sample_ID, Gene, TPM, Species) %>%
  mutate(Type = "VF")

combined <- bind_rows(arg2, vf2)

# Merge group info
combined <- combined %>%
  left_join(meta %>% select(Sample_ID, Group), by = "Sample_ID") %>%
  filter(!is.na(Group))

cat("Combined rows after merge:", nrow(combined), "\n\n")

# ============================================================
# 3. SELECT TOP N GENES (GLOBAL)
# ============================================================

topN <- 50

top_genes <- combined %>%
  group_by(Gene) %>%
  summarise(total = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = topN) %>%
  pull(Gene)

combined_top <- combined %>%
  filter(Gene %in% top_genes)

cat("✔ Top genes selected:", length(top_genes), "\n\n")

# ============================================================
# 4. AGGREGATE TPM (Species × Gene)
# ============================================================

agg <- combined_top %>%
  group_by(Species, Gene) %>%
  summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")

# ============================================================
# 5. BUILD MATRIX
# ============================================================

species_gene_matrix <- agg %>%
  pivot_wider(
    names_from = Gene,
    values_from = TPM,
    values_fill = 0
  )

mat <- as.matrix(species_gene_matrix[,-1])
rownames(mat) <- species_gene_matrix$Species

# ============================================================
# 6. FILTER TOP SPECIES
# ============================================================

top_speciesN <- 30

species_totals <- rowSums(mat)
top_species <- names(sort(species_totals, decreasing = TRUE))[1:top_speciesN]

mat <- mat[top_species, , drop = FALSE]

cat("Matrix size:", dim(mat), "\n\n")

# ============================================================
# 7. COLUMN ANNOTATIONS
# ============================================================

# ARG / VF annotation
gene_info <- combined_top %>%
  distinct(Gene, Type)

type_vec <- gene_info$Type[match(colnames(mat), gene_info$Gene)]

type_colors <- c(
  ARG = "#4C72B0",
  VF  = "#DD8452"
)

# Dominant group per gene
group_assign <- combined_top %>%
  group_by(Gene, Group) %>%
  summarise(total = sum(TPM), .groups = "drop") %>%
  group_by(Gene) %>%
  slice_max(total, n = 1) %>%
  ungroup()

group_vec <- group_assign$Group[match(colnames(mat), group_assign$Gene)]
group_vec <- factor(group_vec, levels = c("Healthy", "Adenoma", "Cancer"))

group_colors <- c(
  Healthy="#ffcc00",
  Adenoma="#984ea3",
  Cancer ="#4daf4a"
)
col_ha <- HeatmapAnnotation(
  Type = type_vec,
  Group = group_vec,
  col = list(
    Type  = type_colors,
    Group = group_colors
  ),
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold")
)

# ============================================================
# 8. DRAW HEATMAP
# ============================================================

ht <- Heatmap(
  log1p(mat),
  name = "log(TPM + 1)",
  top_annotation = col_ha,
  
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  
  show_column_names = TRUE,
  column_names_rot = 75,
  column_names_gp = gpar(fontsize = 9),
  
  heatmap_legend_param = list(
    title = "Expression",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  ),
  
  border = TRUE,
  row_title = "Top Species",
  row_title_gp = gpar(fontsize = 13, fontface = "bold"),
  column_title = "Top ARG + VF Genes",
  column_title_gp = gpar(fontsize = 13, fontface = "bold")
)

# ============================================================
# 9. SAVE HIGH-RES PNG
# ============================================================

png(
  filename = "Top_ARG_VF_heatmap.png",
  width = 7000,
  height = 4500,
  res = 600
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(10, 10, 10, 10), "mm")
)

dev.off()

cat("heatmap saved: Top_ARG_VF_heatmap.png\n")

