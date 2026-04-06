## Spatial_transcriptomics ##
## spatial data deconvolution using Robust Cell Type Decomposition (RCTD) and Spacexr Package ##

library(Seurat)
library(spacexr)
library(tidyverse)
library(patchwork)
set.seed(12345)

# load Visium HD data gene expression
colon_spatial <- Load10X_Spatial(data.dir = "spatial/Whole Transcriptome/",
                         filename = "filtered_feature_bc_matrix.h5",
                         assay = "Spatial",
                         slice = "slice1", filter.matrix = TRUE,
                         to.upper = FALSE, image = NULL)
colon_spatial <- NormalizeData(colon_spatial)
DefaultAssay(colon_spatial) <- "Spatial"
colon_spatial <- FindVariableFeatures(colon_spatial)
colon_spatial <- ScaleData(colon_spatial)
colon_spatial <- RunPCA(colon_spatial, assay = "Spatial", reduction.name = "pca.colon.Spatial", verbose = T)
colon_spatial <- FindNeighbors(colon_spatial, reduction = "pca.colon.Spatial", dims = 1:50)
colon_spatial <- FindClusters(colon_spatial, cluster.name = "seurat_cluster.Spatial")
colon_spatial <- RunUMAP(colon_spatial, reduction = "pca.colon.Spatial", reduction.name = "umap.colon.Spatial",
                  return.model = T, dims = 1:50, verbose = T)

## Create the RCTD query object using 'SpatialRNA' function
counts_hd <- colon_spatial[["Spatial"]]$counts
colon_cells_hd <- colnames(colon_spatial[["Spatial"]])
coords <- GetTissueCoordinates(colon_spatial)[colon_cells_hd, 1:2]
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))
## Load in a scRNA-Seq reference datasets
ref <- readRDS("srobject.RDS")
counts <- ref[["RNA"]]$counts
view(ref@meta.data)
cluster <- as.factor(ref@active.ident)
nUMI <- ref$nCount_RNA

# Create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

# Create RCTD object
RCTD <- create.RCTD(query, reference, max_cores = 8) # for parallel processing

# Run RCTD Algorithm
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

# Add results back to Seurat object
colon_spatial <- AddMetaData(colon_spatial, metadata = RCTD@results$results_df)