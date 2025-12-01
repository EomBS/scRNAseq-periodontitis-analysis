# ============================================================
# 01_Preprocessing.R
# Template script for scRNA-seq preprocessing (non-executable)
# ============================================================

# Load required libraries (replace with your environment)
library(Seurat)
library(DoubletFinder)

# ------------------------------------------------------------
# 1. Load raw count matrix
# ------------------------------------------------------------
# Replace "PATH_TO_DATA" with your dataset directory
counts <- Read10X(data.dir = "PATH_TO_DATA")

# Create Seurat object (set your own min.cells/min.features)
seu <- CreateSeuratObject(
  counts = counts,
  project = "PROJECT_NAME",
  min.cells = 3,
  min.features = 200
)

# ------------------------------------------------------------
# 2. Basic QC metrics
# ------------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Example QC filters (modify with your own thresholds)
seu <- subset(
  seu,
  subset = nFeature_RNA > MIN_FEATURE &
    nFeature_RNA < MAX_FEATURE &
    percent.mt < MAX_MT_PERCENT
)

# ------------------------------------------------------------
# 3. Normalization & dimension reduction
# ------------------------------------------------------------
seu <- SCTransform(seu)
seu <- RunPCA(seu)

# Optional exploratory embeddings
# seu <- RunTSNE(seu)
# seu <- RunUMAP(seu, dims = 1:N_PCS)

# ------------------------------------------------------------
# 4. Doublet identification (template only)
# ------------------------------------------------------------
# These steps require manual tuning; placeholders only.

# seu <- FindNeighbors(seu, dims = 1:N_PCS)
# seu <- FindClusters(seu, resolution = RESOLUTION)

# sweep.list <- paramSweep(seu, PCs = 1:N_PCS, sct = TRUE)
# sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
# bcmvn <- find.pK(sweep.stats)

# Estimate expected doublets (replace rate)
# homotypic.prop <- modelHomotypic(seu$seurat_clusters)
# nExp <- round(DOUBLET_RATE * nrow(seu@meta.data))
# nExp.adj <- round(nExp * (1 - homotypic.prop))

# doublet_results <- doubletFinder(
#    seu,
#    PCs = 1:N_PCS,
#    pN = PN_VALUE,
#    pK = PK_VALUE,
#    nExp = nExp.adj,
#    reuse.pANN = FALSE,
#    sct = TRUE
# )

# ------------------------------------------------------------
# 5. Subset singlets
# ------------------------------------------------------------
# seu <- subset(doublet_results, subset = DF.classifications == "Singlet")

# ------------------------------------------------------------
# 6. Additional QC checks (optional)
# ------------------------------------------------------------
# ElbowPlot(seu)
# DimPlot(seu)

# ------------------------------------------------------------
# 7. Marker identification (placeholder)
# ------------------------------------------------------------
# markers <- FindAllMarkers(
#   seu,
#   only.pos = TRUE,
#   min.pct = MIN_PCT,
#   logfc.threshold = LOGFC_THRESHOLD
# )

# write.csv(markers, file = "OUTPUT_MARKERS.csv")

# ------------------------------------------------------------
# End of template script
# Replace all placeholder values with your own data and parameters.
# This script is intentionally non-executable for public sharing.
