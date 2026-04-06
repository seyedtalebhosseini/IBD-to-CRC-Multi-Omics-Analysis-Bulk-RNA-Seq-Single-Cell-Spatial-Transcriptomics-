library(Seurat)
library(spacexr)
library(tidyverse)
library(patchwork)
library(hdf5r)
library(ggplot2)

set.seed(12345)

# CD "CEACAM5", "CENPE", "NUAK1", "CALCRL", "ADAMTS4", "CHI3L1", "SERPINE1", "ABL2", "MAPK10", "MMP3"
# UC "CEACAM5","MMP12","NAMPT","LGALS1","PDGFRA"
colorectalcapture <- readRDS("CRC/colorectalcapture_spatial_uc.RDS")
colonFPPE <- readRDS("CRC/colonFPPE_spatial_uc.RDS")
colorectalcancertranscriptome <- readRDS("CRC/colorectalcancertranscriptome_spatial_uc.RDS")

pdf("SpatialFeaturePlot_UC_colorectalcancertranscriptome.pdf", width = 10, height = 2.8)
SpatialFeaturePlot(colorectalcancertranscriptome, 
                        features = c("CEACAM5","MMP12","NAMPT","LGALS1","PDGFRA") , 
                        crop = TRUE, 
                        pt.size.factor = 2, 
                        ncol = 5, image.alpha = 0,
                   image.scale = "hires")
dev.off()
