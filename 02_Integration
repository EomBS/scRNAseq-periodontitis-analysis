# ============================================================
# 02_Integration.R
# Template for multi-sample integration and clustering (non-executable)
# ============================================================

library(Seurat)
library(harmony)

# ------------------------------------------------------------
# 1. Merge multiple Seurat objects (one per sample / donor)
# ------------------------------------------------------------
# Assume you have a list of Seurat objects: seu_list
# e.g., seu_list <- list(sample1, sample2, sample3, ...)

# Merge into a single Seurat object
# Adjust the variable "grouping" (e.g., orig.ident, sampleID) as needed.
merged <- merge(
  x = seu_list[[1]],
  y = seu_list[-1]
)

# Optionally, ensure grouping metadata exists (e.g., per sample)
# merged$group <- merged@meta.data$orig.ident  # or your own grouping column

# ------------------------------------------------------------
# 2. Split by sample/group and run SCTransform per subset
# ------------------------------------------------------------
obj.list <- SplitObject(merged, split.by = "orig.ident")  # replace with your column name

for (i in seq_along(obj.list)) {
  obj.list[[i]] <- SCTransform(obj.list[[i]])
}

# ------------------------------------------------------------
# 3. Select integration features and reassign variable features
# ------------------------------------------------------------
# Here we use 3,000 features as an example; adjust as needed.
var.features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 3000
)

VariableFeatures(merged) <- var.features

# ------------------------------------------------------------
# 4. PCA on merged object
# ------------------------------------------------------------
merged <- RunPCA(merged, verbose = FALSE)

# ------------------------------------------------------------
# 5. Harmony integration
# ------------------------------------------------------------
# Replace "orig.ident" with the metadata column indicating sample / batch ID.
# Replace 1:N_PCS with the number of PCs you wish to use.
merged <- RunHarmony(
  merged,
  assay.use    = "SCT",
  group.by.vars = "orig.ident",
  dims.use     = 1:N_PCS,
  max.iter.harmony = 50
)

# ------------------------------------------------------------
# 6. Graph construction and clustering
# ------------------------------------------------------------
merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims      = 1:N_PCS
)

# Optional: visualize elbow plot to choose PCs
# ElbowPlot(merged)

# Replace CLUSTER_RESOLUTION with your chosen resolution
merged <- FindClusters(
  object     = merged,
  resolution = CLUSTER_RESOLUTION
)

# ------------------------------------------------------------
# 7. UMAP on Harmony embeddings
# ------------------------------------------------------------
merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims      = 1:N_PCS
)

# Optional: save elbow plot / UMAP externally (example, commented)
# pdf("PATH_TO_OUTPUT/elbow_plot.pdf")
# ElbowPlot(merged)
# dev.off()

# ------------------------------------------------------------
# End of template script.
# All object names, metadata columns, and parameter values
# (N_PCS, CLUSTER_RESOLUTION, paths) must be customized
# for the user's own dataset.
