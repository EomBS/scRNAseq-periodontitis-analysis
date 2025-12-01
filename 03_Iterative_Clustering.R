# ============================================================
# 03_Iterative_Clustering.R
# Template for subset-level reclustering (non-executable)
# ============================================================

library(Seurat)
library(ggplot2)

# 'immune' should be an integrated immune Seurat object
# with an annotation column (e.g., "Annotation") and
# optional group label (e.g., "disease", "group").

# Example:
# immune <- readRDS("PATH_TO_INTEGRATED_IMMUNE_OBJECT.rds")

# ------------------------------------------------------------
# 1. Subset a major immune population (e.g., NK/T cells)
# ------------------------------------------------------------
# Replace "Annotation" and "NK/T" with your own metadata field and label.
subset_cells <- subset(
  immune,
  cells = rownames(
    immune@meta.data[
      immune@meta.data$Annotation == "NK/T",
    ]
  )
)

# ------------------------------------------------------------
# 2. Re-normalize and recluster within the subset
# ------------------------------------------------------------
subset_cells <- SCTransform(subset_cells)

# Use Harmony reduction if already computed, or PCA otherwise.
# Adjust dims and reduction according to your analysis.
subset_cells <- FindNeighbors(
  subset_cells,
  dims      = 1:N_PCS,
  reduction = "harmony"   # or "pca"
)

# Optional: inspect elbow plot for PCs
# ElbowPlot(subset_cells)

subset_cells <- FindClusters(
  subset_cells,
  resolution = SUBSET_RESOLUTION
)

subset_cells <- RunUMAP(
  subset_cells,
  dims      = 1:N_PCS,
  reduction = "harmony",
  verbose   = FALSE
)

# ------------------------------------------------------------
# 3. Visualization (examples)
# ------------------------------------------------------------
# Basic UMAP with cluster labels
# DimPlot(subset_cells, label = TRUE) + NoLegend()

# Split by condition / group (if available)
# DimPlot(subset_cells, label = TRUE, split.by = "group") + NoLegend()

# ------------------------------------------------------------
# 4. Marker gene identification for subset clusters
# ------------------------------------------------------------
# Adjust thresholds as needed.
subset_markers <- FindAllMarkers(
  subset_cells,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

# Example: export marker table (path is a placeholder)
# write.csv(subset_markers, file = "PATH_TO_OUTPUT/subset_markers.csv")

# ------------------------------------------------------------
# End of template script.
# All column names (Annotation, group), labels (e.g., "NK/T"),
# parameter values (N_PCS, SUBSET_RESOLUTION), and output paths
# must be adapted to your own dataset.
