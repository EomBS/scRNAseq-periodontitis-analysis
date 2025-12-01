# ============================================================
# 04_Trajectory.R
# Template for trajectory inference with Monocle3 (non-executable)
# ============================================================

library(Seurat)
library(monocle3)
library(ggplot2)
# library(viridis) # optional for coloring

# ------------------------------------------------------------
# 1. Convert Seurat object to Monocle3 cell_data_set
# ------------------------------------------------------------
# 'seu' should be a Seurat object containing the subset of interest
# e.g., B cell / plasma cell subset integrated across conditions.
# seu <- readRDS("PATH_TO_SEURAT_SUBSET.rds")

cds <- as.cell_data_set(seu)

# ------------------------------------------------------------
# 2. Preprocess and reduce dimensions (UMAP)
# ------------------------------------------------------------
# Adjust num_dim and UMAP parameters as needed for your dataset.
cds <- preprocess_cds(cds, num_dim = 50)

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  umap.n_neighbors = 30,
  umap.min_dist    = 0.3
)

# ------------------------------------------------------------
# 3. Clustering within Monocle3
# ------------------------------------------------------------
# Resolution and clustering behavior can be tuned.
cds <- cluster_cells(cds, resolution = 1e-3)

# ------------------------------------------------------------
# 4. Learn trajectory graph
# ------------------------------------------------------------
cds <- learn_graph(
  cds,
  use_partition = FALSE,
  learn_graph_control = list(ncenter = 250)
)

# ------------------------------------------------------------
# 5. Define root cells and order cells in pseudotime
# ------------------------------------------------------------
# Replace "Annotation" and "ROOT_LABEL" with your metadata and root population.
# For example, you may choose naive B cells, early precursors, etc.
root_cells <- colnames(cds)[cds@colData$Annotation == "ROOT_LABEL"]

cds <- order_cells(
  cds,
  root_cells = root_cells
)

# ------------------------------------------------------------
# 6. Basic pseudotime visualization in Monocle3
# ------------------------------------------------------------
# Plot pseudotime on UMAP with trajectory graph
# plot_cells(
#   cds,
#   color_cells_by       = "pseudotime",
#   label_branch_points  = TRUE,
#   label_leaves         = TRUE,
#   cell_size            = 1
# )

# Plot pseudotime without trajectory graph overlay
# plot_cells(
#   cds,
#   color_cells_by       = "pseudotime",
#   show_trajectory_graph = FALSE,
#   cell_size            = 1.5
# )
# ------------------------------------------------------------
# End of template script.
# All metadata fields (Annotation, group), gene names,
# and paths must be adapted to the user's own dataset.
