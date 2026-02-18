# ============================================================
# Integrated CRC Marker Reproducibility and Network Pipeline
# ============================================================
set.seed(42)

# ---- 1. Load packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(igraph)
  library(ggraph)
  library(ggvenn)
  library(ggplot2)
})

# ---- 2. Load data ----
tpm <- read.csv("TPM_raw_filtered_bact.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("CRC_metadata_aligned_bact.csv", stringsAsFactors = FALSE)
ref <- read.csv("CRC_reference_marker_panel.csv")

stopifnot(all(c("Sample_ID", "Group") %in% colnames(meta)))
rownames(meta) <- meta$Sample_ID
tpm[is.na(tpm)] <- 0

# ---- 3. Clean species names ----
rownames(tpm) <- gsub("_", " ", rownames(tpm))
ref$clean_species <- gsub("_", " ", tolower(ref$clean_species))
species_detected <- tolower(rownames(tpm))

# ---- 4. Identify overlap ----
ref_species <- unique(ref$clean_species)
shared_markers <- intersect(species_detected, ref_species)
write.csv(data.frame(Species = shared_markers), "CRC_Reproducible_Markers.csv", row.names = FALSE)
cat("✅ Shared CRC markers found:", length(shared_markers), "\n")

# ---- 5. Venn diagram (Your data vs Literature markers) ----

venn_plot <- ggvenn(
  list(
    "Your Dataset" = species_detected,
    "Reference CRC Markers" = ref_species
  ),
  fill_color = c("blue", "gold"),   # optional custom colors
  stroke_size = 1.2,
  set_name_size = 5,
  text_size = 5
)

# Save with white background (non-transparent)
ggsave("CRC_Marker_Overlap_Venn.png",
       plot = venn_plot,
       width = 5, height = 4, dpi = 600,
       bg = "white")     # <-- solid background


# ============================================================
# 6. Heatmap of Reproducible CRC Markers
# ============================================================
library(ComplexHeatmap)
library(circlize)
library(grid)

groups <- c("Healthy", "Adenoma", "Cancer")

# ---- Map marker categories ----
marker_map <- setNames(rep("Other", length(species_detected)), species_detected)
marker_map[species_detected %in% tolower(ref$clean_species[ref$source == "Wirbel_core"])] <- "Core (Wirbel)"
marker_map[species_detected %in% tolower(ref$clean_species[ref$source == "Wirbel_extended"])] <- "Extended"
marker_map[species_detected %in% tolower(ref$clean_species[ref$source == "Piccino_stage"])] <- "Stage-associated"

# ---- Keep only microbial markers ----
tpm_markers <- tpm[tolower(rownames(tpm)) %in% names(marker_map[marker_map != "Other"]), ]
meta <- meta[meta$Sample_ID %in% colnames(tpm_markers), ]
rownames(meta) <- meta$Sample_ID

# ---- Calculate mean abundance by group ----
log_tpm <- log10(tpm_markers + 1)
mean_by_group <- sapply(groups, function(g) {
  idx <- which(meta$Group == g)
  if (length(idx) > 0) rowMeans(log_tpm[, idx, drop = FALSE], na.rm = TRUE)
  else rep(NA, nrow(log_tpm))
})
mean_by_group <- as.data.frame(mean_by_group)
rownames(mean_by_group) <- rownames(log_tpm)

# ---- Match rows to their category ----
row_categories <- marker_map[tolower(rownames(mean_by_group))]
row_split_vec <- factor(row_categories,
                        levels = c("Core (Wirbel)", "Extended", "Stage-associated"))

# ---- Remove undefined rows ----
valid_rows <- !is.na(row_split_vec)
mean_by_group <- mean_by_group[valid_rows, ]
row_split_vec <- droplevels(row_split_vec[valid_rows])

# ============================================================
#   Z-score transform (per species)
# ============================================================
scaled_matrix <- t(scale(t(mean_by_group)))     # center and scale each species
scaled_matrix[is.na(scaled_matrix)] <- 0        # replace NAs

# ============================================================
#   Color palette and annotations
# ============================================================

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#01665E", "white", "#F46D43")
)


group_cols <- c("Healthy"="#ffcc00", "Adenoma"="#984ea3", "Cancer"="#4daf4a")

present_groups <- intersect(groups, colnames(mean_by_group))

top_annotation <- HeatmapAnnotation(
  Group = colnames(mean_by_group),
  col = list(Group = group_cols[present_groups]),
  annotation_legend_param = list(title = "Group")
)

# ============================================================
#  Draw Z-score heatmap
# ============================================================
tiff("CRC_Marker_Heatmap_Zscore.tiff", width = 7, height = 7, units = "in", res = 600)
Heatmap(
  as.matrix(scaled_matrix),
  name = "Z-score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_split = row_split_vec,
  row_title = levels(row_split_vec),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 7),
  top_annotation = top_annotation,
  column_names_gp = gpar(fontsize = 10, fontface = "bold")
)
dev.off()

cat("CRC marker Z-score heatmap saved as 'CRC_Marker_Heatmap_Zscore.tiff'\n")


# ============================================================
# 7. Co-occurrence Network of Reproducible Markers
# ============================================================

# ============================================================
# Informative CRC Marker Networks (Core + Hubs + Communities)
# ============================================================

library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)

# ---- 1. Subset reproducible markers ----
tpm_shared <- tpm[tolower(rownames(tpm)) %in% shared_markers, ]
tpm_shared <- t(tpm_shared)

# ---- 2. Correlation ----
cor_matrix <- cor(tpm_shared, method = "spearman")
diag(cor_matrix) <- 0

# keep moderate associations to reveal modules
threshold <- 0.5
adj_matrix <- ifelse(abs(cor_matrix) > threshold, cor_matrix, 0)
diag(adj_matrix) <- 0
network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE)

# ---- 3. Add attributes ----
V(network)$degree <- degree(network)
V(network)$source <- ifelse(
  tolower(V(network)$name) %in% tolower(ref$clean_species[ref$source == "Wirbel_core"]), "Core (Wirbel)",
  ifelse(tolower(V(network)$name) %in% tolower(ref$clean_species[ref$source == "Wirbel_extended"]), "Extended",
         ifelse(tolower(V(network)$name) %in% tolower(ref$clean_species[ref$source == "Piccino_stage"]), "Stage-associated", "Other"))
)
E(network)$sign <- ifelse(E(network)$weight > 0, "Positive", "Negative")

# ---- 4. Remove isolates (species with no edges) ----
network <- delete.vertices(network, degree(network) == 0)

# ---- 5. Community detection ----
comm <- cluster_louvain(network)
V(network)$community <- membership(comm)

# ---- Color palettes ----
marker_cols <- c("Core (Wirbel)"="#d62728",
                 "Extended"="#ff7f0e",
                 "Stage-associated"="#2ca02c",
                 "Other"="grey80")
edge_cols <- c("Positive"="#6baed6", "Negative"="#fb6a4a")
comm_cols <- RColorBrewer::brewer.pal(max(V(network)$community), "Set3")

# ============================================================
# A. Core Co-occurrence Subnetwork (main cluster)
# ============================================================
core_comm <- which.max(table(V(network)$community))
core_net <- induced_subgraph(network, vids = V(network)[community == core_comm])

tiff("CRC_Marker_Network_CoreCluster.png", width = 6, height = 6, units = "in", res = 600)
ggraph(core_net, layout = "fr") +
  geom_edge_link(aes(color = sign, alpha = abs(weight)), show.legend = TRUE) +
  scale_edge_color_manual(values = edge_cols) +
  scale_edge_alpha(range = c(0.2, 0.8)) +
  geom_node_point(aes(size = degree, color = source), alpha = 0.9) +
  scale_color_manual(values = marker_cols) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() +
  labs(title = "Core Co-occurrence Subnetwork of CRC Markers")
dev.off()

# ============================================================
# B. Hub-Focused Network (top 15)
# ============================================================
hub_nodes <- names(sort(degree(network), decreasing = TRUE))[1:20]
hub_net <- induced_subgraph(network, vids = V(network)[name %in% hub_nodes])

tiff("CRC_Marker_Network_TopHubs.png", width = 5, height = 5, units = "in", res = 600)
ggraph(hub_net, layout = "fr") +
  geom_edge_link(aes(color = sign, alpha = abs(weight)), show.legend = FALSE) +
  geom_node_point(aes(size = degree, color = source), alpha = 1.0) +
  scale_color_manual(values = marker_cols) +
  geom_node_label(aes(label = name), size = 2, fill = "white", repel = TRUE) +
  theme_void() +
  labs(title = "Top 15 CRC Marker Hubs (Most Connected Species)")
dev.off()

# ============================================================
# C. Full Community-Colored Network
# ============================================================
tiff("CRC_Marker_Network_Communities.png", width = 7, height = 7, units = "in", res = 600)
ggraph(network, layout = "fr") +
  geom_edge_link(aes(alpha = abs(weight)), color = "grey80") +
  geom_node_point(aes(size = degree, color = as.factor(community)), alpha = 0.9) +
  scale_color_manual(values = comm_cols, name = "Community") +
  geom_node_text(aes(label = ifelse(degree > quantile(degree, 0.75), name, "")),
                 repel = TRUE, size = 2, fontface = "bold") +
  theme_void() +
  labs(title = "Louvain Community Structure of CRC Marker Co-occurrence Network",
       size = "Degree")
dev.off()

cat("Generated: CoreCluster, TopHubs, and Community network figures\n")
# ============================================================
# Network Statistics for Reproducible CRC Marker Network
# ============================================================

library(igraph)
library(dplyr)

# Assuming you already created 'network' from your co-occurrence analysis
# (the one used for the figures)

# ---- 1.  Basic Network-Level Statistics ----
num_nodes <- gorder(network)
num_edges <- gsize(network)
density_val <- edge_density(network)
avg_degree <- mean(degree(network))
avg_path <- mean_distance(network, directed = FALSE)
modularity_val <- modularity(cluster_louvain(network))

net_summary <- data.frame(
  Metric = c("Nodes", "Edges", "Density", "Average Degree",
             "Average Path Length", "Modularity"),
  Value = round(c(num_nodes, num_edges, density_val,
                  avg_degree, avg_path, modularity_val), 4)
)
print(net_summary)

write.csv(net_summary, "CRC_Marker_Network_Summary.csv", row.names = FALSE)
cat("Network-level summary saved as 'CRC_Marker_Network_Summary.csv'\n")


# ---- 2.  Node-Level (Centrality) Statistics ----
node_stats <- data.frame(
  Species = V(network)$name,
  Source = V(network)$source,
  Degree = degree(network),
  Betweenness = betweenness(network, normalized = TRUE),
  Closeness = closeness(network, normalized = TRUE),
  Eigenvector = evcent(network)$vector,
  ClusteringCoeff = transitivity(network, type = "local", isolates = "zero")
)

# Sort by Eigenvector centrality (most influential)
node_stats <- node_stats %>%
  arrange(desc(Eigenvector))

head(node_stats, 10)   # preview top influential species

write.csv(node_stats, "CRC_Marker_Network_Node_Statistics.csv", row.names = FALSE)
cat("Node-level centrality statistics saved as 'CRC_Marker_Network_Node_Statistics.csv'\n")


# ---- 3. Identify Topological Hubs ----
top_hubs <- node_stats %>%
  arrange(desc(Degree)) %>%
  slice(1:10)

write.csv(top_hubs, "CRC_Marker_Top10_Hubs.csv", row.names = FALSE)
cat("top 10 most connected hub species saved as 'CRC_Marker_Top10_Hubs.csv'\n")


# ---- 4. Global Network Descriptors ----
global_measures <- data.frame(
  Network_Measure = c("Average Degree", "Graph Density", "Average Clustering Coefficient",
                      "Network Diameter", "Modularity (Louvain)"),
  Value = c(mean(node_stats$Degree),
            density_val,
            transitivity(network, type = "global"),
            diameter(network, directed = FALSE),
            modularity_val)
)
write.csv(global_measures, "CRC_Marker_Global_Metrics.csv", row.names = FALSE)
cat("Gobal metrics saved as 'CRC_Marker_Global_Metrics.csv'\n")

############
#################################################################################

######################################################################
# STEP — FILTER TPM + CLR MATRICES FOR SHARED CRC MARKERS
######################################################################

cat("\n==============================\n")
cat("   FILTERING TPM & CLR FOR SHARED CRC MARKERS\n")
cat("==============================\n\n")

# ----- 1.  TPM FILTERING -----

# get lowercase species names for TPM rows
tpm_clean <- tolower(rownames(tpm))

# mapping table (original → cleaned)
map_tpm <- data.frame(
  Original = rownames(tpm),
  Clean = tpm_clean,
  stringsAsFactors = FALSE
)

# select rows whose cleaned name is in shared_markers
tpm_shared_only <- tpm[map_tpm$Clean %in% shared_markers, , drop = FALSE]

# restore original rownames
rownames(tpm_shared_only) <- map_tpm$Original[map_tpm$Clean %in% shared_markers]

cat("✔ TPM shared markers dimensions:",
    nrow(tpm_shared_only), "species ×", ncol(tpm_shared_only), "samples\n")

# save TPM filtered matrix
write.csv(tpm_shared_only,
          "TPM_Shared_CRC_Markers_Only.csv",
          row.names = TRUE)
cat("Saved: TPM_Shared_CRC_Markers_Only.csv\n\n")


# ----- 2.  CLR FILTERING (TPM_clr_batch_corrected_v2.csv) -----

cat("Loading CLR batch-corrected matrix (v2)...\n")

species_clr_v2 <- read.csv("TPM_clr_batch_corrected_v2.csv",
                           row.names = 1, check.names = FALSE)

# Determine orientation automatically
# If first row looks like samples (DRR/ERR/SRR), transpose
first_row_clr <- rownames(species_clr_v2)[1]

if (!any(grepl("DRR|SRR|ERR", first_row_clr))) {
  cat("CLR matrix detected as species × samples → transposing...\n")
  species_clr_v2 <- t(species_clr_v2)
} else {
  cat("CLR matrix detected as samples × species.\n")
}

# Now species = columns
clr_species_clean <- tolower(colnames(species_clr_v2))

map_clr <- data.frame(
  Original = colnames(species_clr_v2),
  Clean = clr_species_clean,
  stringsAsFactors = FALSE
)

# keep only shared marker species
clr_keep <- map_clr$Original[map_clr$Clean %in% shared_markers]

clr_shared_only <- species_clr_v2[, clr_keep, drop = FALSE]

cat("CLR shared markers dimensions:",
    nrow(clr_shared_only), "samples ×", ncol(clr_shared_only), "species\n")

# reorder rows to match metadata (optional but recommended)
if ("Sample_ID" %in% colnames(meta)) {
  clr_shared_only <- clr_shared_only[meta$Sample_ID, , drop = FALSE]
}

# save CLR filtered matrix
write.csv(clr_shared_only,
          "CLR_Shared_CRC_Markers_Only.csv",
          row.names = TRUE)
cat("Saved: CLR_Shared_CRC_Markers_Only.csv\n\n")


# ----- 3. Save metadata subset (optional) -----

meta_shared <- meta[meta$Sample_ID %in% colnames(tpm_shared_only), ]
write.csv(meta_shared,
          "Metadata_Shared_CRC_Markers_Samples.csv",
          row.names = FALSE)

cat("Saved: Metadata_Shared_CRC_Markers_Samples.csv\n")
cat("Shared CRC marker matrices successfully generated!\n")
###############################################################################
######################################################################
# STEP — FILTER TPM + CLR FOR UNSHARED SPECIES (Not CRC Markers)
######################################################################

cat("\n==============================\n")
cat("   FILTERING UNSHARED MICROBIOME SPECIES\n")
cat("==============================\n\n")

# ============================================================
# 1. PREP: Clean species names again (TPM + CLR)
# ============================================================

# TPM cleaning mapping
tpm_clean <- tolower(rownames(tpm))
map_tpm <- data.frame(
  Original = rownames(tpm),
  Clean = tpm_clean,
  stringsAsFactors = FALSE
)

# CLR cleaning mapping
species_clr_v2 <- read.csv("TPM_clr_batch_corrected_v2.csv",
                           row.names = 1, check.names = FALSE)

# Detect orientation
first_row_clr <- rownames(species_clr_v2)[1]
if (!any(grepl("DRR|SRR|ERR", first_row_clr))) {
  species_clr_v2 <- t(species_clr_v2)
}

clr_clean <- tolower(colnames(species_clr_v2))
map_clr <- data.frame(
  Original = colnames(species_clr_v2),
  Clean = clr_clean,
  stringsAsFactors = FALSE
)


######################################################################
# 2. IDENTIFY UNSHARED SPECIES
######################################################################

unshared_species <- setdiff(tpm_clean, shared_markers)

cat("✔ Total TPM species:", length(tpm_clean), "\n")
cat("✔ Shared CRC markers:", length(shared_markers), "\n")
cat("✔ Unshared species:", length(unshared_species), "\n\n")


######################################################################
# 3. FILTER TPM: UNSHARED SPECIES ONLY
######################################################################

tpm_unshared <- tpm[map_tpm$Clean %in% unshared_species, , drop = FALSE]
rownames(tpm_unshared) <- map_tpm$Original[map_tpm$Clean %in% unshared_species]

write.csv(tpm_unshared, "TPM_Unshared_Species_Only.csv", row.names = TRUE)
cat("Saved: TPM_Unshared_Species_Only.csv\n")


######################################################################
# 4️. FILTER CLR: UNSHARED SPECIES ONLY
######################################################################

clr_keep_unshared <- map_clr$Original[map_clr$Clean %in% unshared_species]

clr_unshared <- species_clr_v2[, clr_keep_unshared, drop = FALSE]

# Optionally reorder rows to metadata order
if ("Sample_ID" %in% colnames(meta)) {
  clr_unshared <- clr_unshared[meta$Sample_ID, , drop = FALSE]
}

write.csv(clr_unshared,
          "CLR_Unshared_Species_Only.csv",
          row.names = TRUE)

cat("Saved: CLR_Unshared_Species_Only.csv\n")


######################################################################
# 5️⃣ SAVE SUBSETTED METADATA (same samples)
######################################################################

meta_unshared <- meta[meta$Sample_ID %in% colnames(tpm_unshared), ]
write.csv(meta_unshared,
          "Metadata_Unshared_Species_Samples.csv",
          row.names = FALSE)

cat("Saved: Metadata_Unshared_Species_Samples.csv\n")
cat("Unshared species matrices successfully generated!\n")

