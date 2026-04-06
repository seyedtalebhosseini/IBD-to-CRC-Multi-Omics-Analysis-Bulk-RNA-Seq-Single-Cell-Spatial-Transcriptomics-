### This script is analysis of enrichment (biological_process (BP),
### molecular_function (MF), cellular_component (CC) and kegg_pathway (KEGG)) 
### for markers of each cell type using clusterProfiler package (4.10.1) ###

## ClusterProfiler_package_Single_Cell ##
## Load packages ##
library(tidyverse)
library(tidytree)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(enrichplot)
library(gridExtra)

## Group_1:Crohns_Disease ##
Epithelial_cell <- readLines("E:/biodata/Epithelial cells.txt")
T_cell <- readLines("E:/biodata/T cells.txt")
NK_cell <- readLines("E:/biodata/NK cell.txt")
DC <- readLines("E:/biodata/DC.txt")
Macrophage <- readLines("E:/biodata/Macrophage.txt")
Astrocyte <- readLines("E:/biodata/Astrocyte.txt")
B_cell <- readLines("E:/biodata/B cell.txt")
Monocyte <- readLines("E:/biodata/Monocyte.txt")
Tissue_stem_cell <- readLines("E:/biodata/Tissue stem cells.txt")
Neutrophils <- readLines("E:/biodata/Neutrophils.txt")
Endothelial_cell <- readLines("E:/biodata/Endothelial cells.txt")

gene_list <- list(Epithelial_cell=Epithelial_cell, T_cell=T_cell, NK_cell=NK_cell,
                  DC=DC, Macrophage=Macrophage, Astrocyte=Astrocyte, B_cell=B_cell,
                  Monocyte=Monocyte, Tissue_stem_cell=Tissue_stem_cell, Neutrophils=Neutrophils,
                  Endothelial_cell=Endothelial_cell)


## Biological_Process ##
Enrich_BP_Normal_Crohns_Disease <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

## DotPlot_Biological_Process ##
pdf("E:/biodata/Enrich_BP_Normal_Crohns_Disease.pdf", width = 15, height = 15)
dotplot(Enrich_BP_Normal_Crohns_Disease, showCategory = 2) + 
  ggtitle("(A) GO_BP of CellTypes - Normal vs Crohns_Disease") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_BP_Normal_Crohns_Disease, "E:/biodata/Enrich_BP_Normal_Crohns_Disease.csv")

## Molecular_Function ##
Enrich_MF_Normal_Crohns_Disease <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  pAdjustMethod = "BH"
)

## DotPlot_Molecular_Function ##
pdf("E:/biodata/Enrich_MF_Normal_Crohns_Disease.pdf", width = 15, height = 10)
dotplot(Enrich_MF_Normal_Crohns_Disease, showCategory = 2) + 
  ggtitle("(B) GO_MF of CellTypes - Normal vs Crohns_Disease") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_MF_Normal_Crohns_Disease, "E:/biodata/Enrich_MF_Normal_Crohns_Disease.csv")

## Cellular_Component ##
Enrich_CC_Normal_Crohns_Disease <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH"
)

## DotPlot_Cellular_Component ##
pdf("E:/biodata/Enrich_CC_Normal_Crohns_Disease.pdf", width = 15, height = 10)
dotplot(Enrich_CC_Normal_Crohns_Disease, showCategory = 2) + 
  ggtitle("(C) GO_CC of CellTypes - Normal vs Crohns_Disease") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_CC_Normal_Crohns_Disease, "E:/biodata/Enrich_CC_Normal_Crohns_Disease.csv")

## KEGG_Pathway ##
gene_list_entrez <- lapply(gene_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
})

gene_list_entrez <- lapply(gene_list_entrez, function(df) df$ENTREZID)

Enrich_KEGG_Normal_Crohns_Disease <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichKEGG",
  pAdjustMethod = "BH",
  organism = "hsa"
)

## DotPlot_KEGG_Pathway ##
pdf("E:/biodata/Enrich_KEGG_Normal_Crohns_Disease.pdf", width = 15, height = 10)
dotplot(Enrich_KEGG_Normal_Crohns_Disease, showCategory = 2) + 
  ggtitle("(D) KEGG of CellTypes - Normal vs Crohns_Disease") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_KEGG_Normal_Crohns_Disease, "E:/biodata/Enrich_KEGG_Normal_Crohns_Disease.csv")

##########################################################################
## Group_2:Ulcerative_Colitis ##
B_cell <- readLines("E:/biodata/B cell.txt")
T_cell <- readLines("E:/biodata/T cell.txt")
Epithelial_cell <- readLines("E:/biodata/Epithelial cell.txt")
NK_cell <- readLines("E:/biodata/NK cell.txt")
Neutrophils <- readLines("E:/biodata/Neutrophils.txt")
Macrophage <- readLines("E:/biodata/Macrophage.txt")
Tissue_stem_cell <- readLines("E:/biodata/Tissue stem cell.txt")
Hepatocytes <- readLines("E:/biodata/Hepatocytes.txt")
Endothelial_cell <- readLines("E:/biodata/Endothelial cell.txt")


gene_list <- list(B_cell=B_cell, T_cell=T_cell, Epithelial_cell=Epithelial_cell,
                  NK_cell=NK_cell, Neutrophils=Neutrophils, Macrophage=Macrophage,
                  Tissue_stem_cell=Tissue_stem_cell, Hepatocytes=Hepatocytes, Endothelial_cell=Endothelial_cell)

## Biological_Process ##
Enrich_BP_Normal_Ulcerative_Colitis <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH"
)

## DotPlot_Biological_Process ##
pdf("E:/biodata/Enrich_BP_Normal_Ulcerative_Colitis.pdf", width = 15, height = 15)
dotplot(Enrich_BP_Normal_Ulcerative_Colitis, showCategory = 2) + 
  ggtitle("(E) GO_BP of CellTypes - Normal vs Ulcerative_Colitis") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_BP_Normal_Ulcerative_Colitis, "E:/biodata/Enrich_BP_Normal_Ulcerative_Colitis.csv")

## Molecular_Function ##
Enrich_MF_Normal_Ulcerative_Colitis <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  pAdjustMethod = "BH"
)
## DotPlot_Molecular_Function ##
pdf("E:/biodata/Enrich_MF_Normal_Ulcerative_Colitis.pdf", width = 15, height = 10)
dotplot(Enrich_MF_Normal_Ulcerative_Colitis, showCategory = 2) + 
  ggtitle("(F) GO_MF of CellTypes - Normal vs Ulcerative_Colitis") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_MF_Normal_Ulcerative_Colitis, "E:/biodata/Enrich_MF_Normal_Ulcerative_Colitis.csv")

## Cellular_Component ##
Enrich_CC_Normal_Ulcerative_Colitis <- compareCluster(
  geneCluster = gene_list,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  pAdjustMethod = "BH"
)

## DotPlot_Cellular_Component ##
pdf("E:/biodata/Enrich_CC_Normal_Ulcerative_Colitis.pdf", width = 15, height = 10)
dotplot(Enrich_CC_Normal_Ulcerative_Colitis, showCategory = 2) + 
  ggtitle("(G) GO_CC of CellTypes - Normal vs Ulcerative_Colitis") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_CC_Normal_Ulcerative_Colitis, "E:/biodata/Enrich_CC_Normal_Ulcerative_Colitis.csv")

## KEGG_Pathway ##
gene_list_entrez <- lapply(gene_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
})

gene_list_entrez <- lapply(gene_list_entrez, function(df) df$ENTREZID)


Enrich_KEGG_Normal_Ulcerative_Colitis <- compareCluster(
  geneCluster = gene_list_entrez,
  fun = "enrichKEGG",
  pAdjustMethod = "BH",
  organism = "hsa"
)

## DotPlot_KEGG_Pathway ##
pdf("E:/biodata/Enrich_KEGG_Normal_Ulcerative_Colitis.pdf", width = 15, height = 10)
dotplot(Enrich_KEGG_Normal_Ulcerative_Colitis, showCategory = 2) + 
  ggtitle("(H) KEGG of CellTypes - Normal vs Ulcerative_Colitis") +
  theme(plot.title = element_text(size = 25))
dev.off()
write.csv(Enrich_KEGG_Normal_Ulcerative_Colitis, "E:/biodata/Enrich_KEGG_Normal_Ulcerative_Colitis.csv")


#### FINISHED ###