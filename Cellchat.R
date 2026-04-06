library(CellChat)
library(Seurat)
library(reticulate)
# 1.Create CellChat Object
seurat_object <- readRDS("Object.RDS")
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
labels <- Idents(seurat_object)
meta <- data.frame(group = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
# 2.Load the Ligand-receptor interaction database
CellChatDB.human <- CellChatDB.human
CellChatDB.mouse <- CellChatDB.mouse
CellChatDB.zebrafish <- CellChatDB.zebrafish
## Show the ligand-receptor categories
showDatabaseCategory(CellChatDB.human)
showDatabaseCategory(CellChatDB.mouse)
showDatabaseCategory(CellChatDB.zebrafish)
# 3.Set all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB.human
# 4.Set a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
CellChatDB.use <- subsetDB(CellChatDB.human, search = "ECM-Receptor")
CellChatDB.use <- subsetDB(CellChatDB.human, search = "Cell-Cell Contact")
# 5.Update CellChatDB
interaction <- CellChatDB.human$interaction
complex <- CellChatDB.human$complex
cofactor <- CellChatDB.human$cofactor
geneInfo <- CellChatDB.human$geneInfo
write.csv(interaction, file = "cellchat/Ulcerative_Colitis/interaction.csv")
write.csv(complex, file = "cellchat/Ulcerative_Colitis/complex.csv")
write.csv(cofactor, file = "cellchat/Ulcerative_Colitis/cofactor.csv")
write.csv(geneInfo, file = "cellchat/Ulcerative_Colitis/geneInfo.csv")
# 6.Adding users curated ligand-receptor pairs
interaction <- read.csv(file = "cellchat/Ulcerative_Colitis/interaction.csv", row.names = 1)
complex <- read.csv(file = "cellchat/Ulcerative_Colitis/complex.csv", row.names = 1)
cofactor <- read.csv(file = "cellchat/Ulcerative_Colitis/cofactor.csv", row.names = 1)
geneInfo <- read.csv(file = "cellchat/Ulcerative_Colitis/geneInfo.csv", row.names = 1)
CellChatDB.human.updated <- list()
CellChatDB.human.updated$interaction <- interaction
CellChatDB.human.updated$complex <- complex
CellChatDB.human.updated$cofactor <- cofactor
CellChatDB.human.updated$geneInfo <- geneInfo
CellChatDB.use <- CellChatDB.human.updated
# 7.Add the Secreted Signaling database in the CellChat Object
cellchat@DB <- CellChatDB.use
# 8.Subset and pre-processing the expression data
## Subset the expression data to use less RAM
cellchat <- subsetData(cellchat)
# 9.Pre-processing the expression data
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
table(cellchat@idents)
## Optional:Project gene expression data onto protein-protein interaction (PPI)
cellchat <- projectData(cellchat, PPI.human)
# 10.Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# 11.Filter out the cell-cell communication if there are only few number of cells
## In certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
# 12.Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)
# 13.Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat@net$count
cellchat@net$weight
# 14.Visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")

dev.off()

# 15.Examine the signaling sent from each cell group
pdf("Ulcerative_Colitis/signal.pdf", width = 200, height = 150)
par(mfrow = c(4, 4), xpd = TRUE, mar = c(30, 30, 30, 30))
mat <- cellchat@net$weight
groupSize <- as.numeric(table(cellchat@idents))
scaled_mat <- mat / max(mat) * 0.3

for (i in 1:nrow(mat)) {
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]  
  
  
  netVisual_circle(
    mat2, 
    vertex.weight = groupSize, 
    weight.scale = F, 
    edge.weight.max = max(mat), 
    edge.width.max = 1, 
    vertex.label.cex = 15
  )
  
  
  mtext(rownames(mat)[i], side = 3, cex = 15, font = 40)
}
dev.off()

p <- netVisual_heatmap(cellchat, signaling = "Name of signaling", color.heatmap = "Reds")
par(mfrow = c(1,1), xpd = TRUE)
par(cex = 0.5)
netVisual_individual(cellchat, signaling = "COLLAGEN", layout = "circle")
netVisual_individual(cellchat, signaling = "LAMININ", layout = "circle")
netVisual_individual(cellchat, signaling = "APP", layout = "circle")
netVisual_individual(cellchat, signaling = "CypA", layout = "circle")
netVisual_individual(cellchat, signaling = "CD99", layout = "circle")
netVisual_individual(cellchat, signaling = "GALECTIN", layout = "circle")
netVisual_individual(cellchat, signaling = "VISFATIN", layout = "circle")
netVisual_individual(cellchat, signaling = "MIF", layout = "circle")
dev.off()