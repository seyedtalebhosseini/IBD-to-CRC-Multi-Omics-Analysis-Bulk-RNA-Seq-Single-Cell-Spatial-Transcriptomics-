### This script is to analysis of trajectory for each cell type using monocle package (2.30.0) ###

## library required packages ...
library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(patchwork)

srobject <- readRDS("object.RDS")
meta <- srobject@meta.data
meta$active.ident <- srobject@active.ident
exp <- srobject@assays$integrated@data
genes <- as.data.frame(rownames(exp))
colnames(genes) <- "gene_short_name"
rownames(genes) <- genes$gene_short_name
set.seed(6)
srobject <- as.data.frame(srobject@assays$integrated@data)
downsample_cells <- as.vector(sample_n(srobject, 2000, replace = FALSE))
downsample_cells <- sample(colnames(exp), 2000)
exp_data_down <- exp[,downsample_cells]
meta_data_down <- meta[downsample_cells,]

pd <- new("AnnotatedDataFrame", data = meta_data_down)
fd <- new("AnnotatedDataFrame", data = genes)
class(exp_data_down)
srobject_CDS <- newCellDataSet(exp_data_down, phenoData = pd,
                               featureData = fd, lowerDetectionLimit = 0.5,
                               expressionFamily = negbinomial.size())

srobject_CDS <- estimateSizeFactors(srobject_CDS)
srobject_CDS <- estimateDispersions(srobject_CDS)

srobject_CDS <- detectGenes(srobject_CDS, min_expr = 0.1)
print(head(fData(srobject_CDS)))

srobject_CDS <- reduceDimension(srobject_CDS, max_components = 2, num_dim = 5, 
                                reduction_method = "tSNE", verbose = TRUE)

srobject_CDS <- clusterCells(srobject_CDS, num_clusters = 4)
plot_cell_clusters(srobject_CDS)
pData(srobject_CDS)$my_colour <- pData(srobject_CDS)$Cluster == 1 | pData(srobject_CDS)$Cluster == 2 | pData(srobject_CDS)$Cluster == 3 | pData(srobject_CDS)$Cluster == 4
plot_cell_clusters(srobject_CDS, color_by = 'my_colour')
expressed_genes <- row.names(subset(fData(srobject_CDS), num_cells_expressed >- 100))
srobject_cds_subset <- srobject[expressed_genes, pData(srobject_CDS)$my_colour]
srobject_cds_subset <- detectGenes(srobject_CDS, min_expr = 0.1)
fData(srobject_cds_subset)$use_for_ordering <- fData(srobject_cds_subset)$num_cells_expressed > 0.05 * ncol(srobject_cds_subset)
srobject_cds_subset <- reduceDimension(srobject_cds_subset, max_components = 2,
                                       norm_method = 'log',
                                       num_dim = 10,
                                       reduction_method = 'tSNE', verbose = TRUE)

srobject_cds_subset <- clusterCells(srobject_cds_subset, verbose = FALSE)
plot_cell_clusters(srobject_cds_subset, color_by = "Cluster")
clustering_DEG_genes <- differentialGeneTest(srobject_cds_subset, fullModelFormulaStr = "~active.ident", cores = 8)
srobject_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:5000]
srobject_cds_subset <- setOrderingFilter(srobject_cds_subset, ordering_genes = srobject_ordering_genes)
srobject_cds_subset <- reduceDimension(srobject_cds_subset, max_components = 2, method = "DDRTree")
srobject_cds_subset <- orderCells(srobject_cds_subset)

srobject_cds_subset$orig.ident
srobject_cds_subset$active.ident
srobject_cds_subset$seurat_clusters

pdf("E:/biodata/Cell_density.pdf", width = 30, height = 15)
p1=plot_cell_trajectory(srobject_cds_subset, color_by = "seurat_clusters", cell_size = 3) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Clusters_Normal_Ulcerative_Colitis")
p2=plot_cell_trajectory(srobject_cds_subset, color_by = "active.ident", cell_size = 3) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Cell_Density_Distribution")

p3=plot_cell_trajectory(srobject_cds_subset, color_by = "Pseudotime", cell_size = 3) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")

p4=plot_cell_trajectory(srobject_cds_subset, color_by = "Class", cell_size = 3) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")

p5=plot_cell_trajectory(srobject_cds_subset, color_by = "State", cell_size = 3) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")

p6=plot_cell_trajectory(srobject_cds_subset, color_by = "State", cell_size = 3) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State") + facet_wrap(~State, nrow=1)
p1+p2+p3+p4+p5+p6
dev.off()


pdf("E:/biodata/pseudotime_trajectories_withD.pdf", width = 25, height = 15)
p1=plot_cell_trajectory(Tissue_Stem_Cells, markers= "LGALS1", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "(D) Tissue stem cells - Ulcerative_Colitis")

p2=plot_cell_trajectory(Endothelial_Cells, markers= "PTPRB", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Endothelial cells  - Ulcerative_Colitis")

p3=plot_cell_trajectory(Macrophage, markers= "MMP12", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Macrophage - Ulcerative_Colitis")


p4=plot_cell_trajectory(Macrophage, markers= "NAMPT", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Macrophage - Ulcerative_Colitis")

p5=plot_cell_trajectory(Tissue_Stem_Cells, markers= "PDGFRA", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Tissue stem cells - Ulcerative_Colitis")

p6=plot_cell_trajectory(Tissue_Stem_Cells, markers= "ABL2", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Tissue stem cells - Crohns_Disease")

p7=plot_cell_trajectory(Tissue_Stem_Cells, markers= "MAPK10", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Tissue stem cells - Crohns_Disease")

p8=plot_cell_trajectory(Endothelial_Cells, markers= "ADAMTS4", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Endothelial cells - Crohns_Disease")

p9=plot_cell_trajectory(Tissue_Stem_Cells, markers= "MMP3", cell_size = 3,
                        use_color_gradient = TRUE, show_tree = TRUE, 
                        show_backbone = FALSE, cell_link_size = 1, 
                        show_branch_points = TRUE, markers_linear = FALSE) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right")+ labs(title = "Tissue stem cells - Crohns_Disease")
p1+p2+p3+p4+p5
dev.off()
## plot pseudotime trajectories
pdf("E:/biodata/monocle_Normal_Crohns_Disease/pseudotime_trajectories_Pseudotime.pdf", width = 25, height = 15)
p1=plot_cell_trajectory(Endothelial_Cells, color_by =  "Pseudotime",markers="NUAK1", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Endothelial cells")

p2=plot_cell_trajectory(Epithelial_Cells, color_by =  "Pseudotime",markers="CENPE", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Epithelial cells")
p3=plot_cell_trajectory(Tissue_Stem_Cells, color_by =  "Pseudotime",markers="CHI3L1", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Tissue stem cells")
p4=plot_cell_trajectory(Endothelial_Cells, color_by =  "Pseudotime", markers="CALCRL", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Endothelial cells ")
p5=plot_cell_trajectory(Tissue_Stem_Cells, color_by =  "Pseudotime",markers="SERPINE1", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Tissue stem cells")
p6=plot_cell_trajectory(Tissue_Stem_Cells, color_by =  "Pseudotime",markers="ABL2", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Tissue stem cells")
p7=plot_cell_trajectory(Tissue_Stem_Cells, color_by =  "Pseudotime",markers="MAPK10", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Tissue stem cells")
p8=plot_cell_trajectory(Endothelial_Cells, color_by =  "Pseudotime",markers="ADAMTS4", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Endothelial cells ")
p9=plot_cell_trajectory(Tissue_Stem_Cells, color_by =  "Pseudotime",markers="MMP3", cell_size = 3) +
  theme(plot.title = element_text(size = 25),
        text = element_text(size = 50, face = "bold"),
        axis.title = element_text(size=20,face="bold"),
        axis.text.x = element_text(size=20, color="black"),
        axis.text.y = element_text(size=20, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"),
        legend.position = "right") + labs(title = "Tissue stem cells")
p1+p2+p3+p4+p5+p6+p7+p8+p9
dev.off()