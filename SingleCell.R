library(BiocManager)
library(Seurat)
library(dplyr)
library(sctransform)
library(limma)
library(cowplot)
library(ggplot2)

my_count=Read10X("barcodes-features-matrixes.tsv")
dim(my_count)
Normal_1=CreateSeuratObject(my_count, min.cells = 3, min.features = 200, project = "Normal_1")
dim(Normal_1)
#this step performed for each sample

#### STEP 2) MERGING
Normal=merge(Normal_1, Normal_2, ... 
            add.cell.ids = c("Normal_1","Normal_2", ...))
dim(Normal)
remove(my_count)
saveRDS(Normal, file = "Normal.RDS")

#### STEP 3) Filtering by number of features(genes) and percentage of mitochondoria
Normal[["MTpercent"]]=PercentageFeatureSet(Normal, pattern = "^MT-")
Normal$MTpercent

## Step 4) Subseting
pdf(file = "Feature_Normal.pdf", title = "Normal", width = 25)
VlnPlot(Normal, features = c("nFeature_RNA","nCount_RNA","MTpercent"), ncol = 3)
plot1=FeatureScatter(Normal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2=FeatureScatter(Normal, feature1 = "nCount_RNA", feature2 = "MTpercent")
plot1+plot2
Normal=subset(Normal, subset = nFeature_RNA>200 & nFeature_RNA<6000 & MTpercent <50)
VlnPlot(Normal, features = c("nFeature_RNA","nCount_RNA","MTpercent"), ncol = 3)
plot1=FeatureScatter(Normal, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2=FeatureScatter(Normal, feature1 = "nCount_RNA", feature2 = "MTpercent")
plot1+plot2
dev.off()

#Normal
nFeature_RNA<6000
#CD
nFeature_RNA<6000
#UC
nFeature_RNA<5000
## Step 5) Normalization
Normal=NormalizeData(Normal, normalization.method = "LogNormalize", scale.factor = 10000)
Normal[["RNA"]]

## Step 6) Highest Variable Features(Genes)
library(ggplot2)
Normal=FindVariableFeatures(Normal, selection.method = "vst", nfeatures = 2000)
topVariableFeatures=head(VariableFeatures(Normal),10)
plot1=VariableFeaturePlot(Normal)
plot2=LabelPoints(plot1,topVariableFeatures,size=4,xnudge = 0.1,ynudge = 0.1)
pdf(file = "FindVariableFeatures_Normal.pdf", width = 12)
plot2
dev.off()

#All steps (step 1-6) performed for all samples such as normal, CD and UC

## Step 7) Integration Between Samples (Normal and CD)
AnchorFeatures=SelectIntegrationFeatures(list(Normal, CD))
Anchors=FindIntegrationAnchors(list(Normal, CD), anchor.features = AnchorFeatures)
srobject_Normal_Crohns_Disease=IntegrateData(anchorset = Anchors)
saveRDS(Anchors, "anchors.RDS")
## Step 8) Centering and Scaling the Data (Mean=0, Variance=1)
srobject_Normal_Crohns_Disease=ScaleData(srobject_Normal_Crohns_Disease)

# Step 9) Dimension Reduction
srobject_Normal_Crohns_Disease=RunPCA(srobject_Normal_Crohns_Disease, features = VariableFeatures(srobject_Normal_Crohns_Disease))
# for see genes in PCA
print(srobject_Normal_Crohns_Disease[["pca"]],dims = 1:5, nfeatures = 10)
pdf(file = "VizDimLoadings.pdf", width = 12)
VizDimLoadings(srobject_Normal_Crohns_Disease, dims = 1:2, reduction = "pca")
dev.off()
pdf(file = "DimPlot.pdf", width = 12)
DimPlot(srobject_Normal_Crohns_Disease, reduction = "pca", pt.size = 1.5)
dev.off()
pdf(file = "DimHeatmap.pdf", width = 12)
DimHeatmap(srobject_Normal_Crohns_Disease, dims = 1, cells = 300)
dev.off()

## Step 10) Determine number of Dimension
srobject_Normal_Crohns_Disease=JackStraw(srobject_Normal_Crohns_Disease, num.replicate = 100)
srobject_Normal_Crohns_Disease=ScoreJackStraw(srobject_Normal_Crohns_Disease, dims = 1:20)
pdf(file = "JackStrawPlot.pdf", width = 12)
JackStrawPlot(srobject_Normal_Crohns_Disease, dims = 1:20)
dev.off()
pdf(file = "ElbowPlot.pdf", width = 12)
ElbowPlot(srobject_Normal_Crohns_Disease)
dev.off()

## Step 11) Clustering 
srobject_Normal_Crohns_Disease=FindNeighbors(srobject_Normal_Crohns_Disease, dims = 1:20)
srobject_Normal_Crohns_Disease=FindClusters(srobject_Normal_Crohns_Disease, resolution = 0.5)
Idents(srobject_Normal_Crohns_Disease)
length(levels(Idents(srobject_Normal_Crohns_Disease)))

## Step 12) UMAP
srobject_Normal_Crohns_Disease=RunUMAP(srobject_Normal_Crohns_Disease, dims = 1:20)
pdf(file = "UMAP.pdf", width = 70, height = 30)
DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", pt.size = 0.5, label = T, raster = FALSE, label.size = 15, split.by = "orig.ident", ncol = 6) + 
  theme(plot.title = element_text(size = 60),
        text = element_text(size = 130, face = "bold"),
        axis.title = element_text(size=25,face="bold"),
        axis.text.x = element_text(size=25, color="black"),
        axis.text.y = element_text(size=25, color="black"),
        legend.text=element_text(size=25, color="black"),
        legend.title=element_text(size=25, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"))
dev.off()

# UMAP_Split_By_Samples #
pdf(file = "UMAP_Split_By_Samples.pdf", width = 50, height = 20)
DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", pt.size = 0.5, label = T, raster = FALSE, label.size = 15, split.by = "orig.ident", ncol = 6) + 
  theme(plot.title = element_text(size = 40),
        text = element_text(size = 40, face = "bold"),
        axis.title = element_text(size=25,face="bold"),
        axis.text.x = element_text(size=25, color="black"),
        axis.text.y = element_text(size=25, color="black"),
        legend.text=element_text(size=25, color="black"),
        legend.title=element_text(size=25, color="black"),
        axis.line.x = element_line(size=2, color="black"),
        axis.line.y = element_line(size=2, color="black"))
dev.off()

## Step 13) Find Markers
markers_Normal_Crohns_Disease=FindAllMarkers(srobject_Normal_Crohns_Disease, min.pct = 0.25, logfc.threshold = 0.25)
dim(markers_Normal_Crohns_Disease)
max(markers_Normal_Crohns_Disease$avg_log2FC)
min(markers_Normal_Crohns_Disease$avg_log2FC)
View(markers)
View(markers_Normal_Crohns_Disease[which(markers_Normal_Crohns_Disease$avg_log2FC>=1),])
markers_Normal_Crohns_Disease=markers_Normal_Crohns_Disease[which(markers_Normal_Crohns_Disease$avg_log2FC>=1),]
write.csv(markers_Normal_Crohns_Disease,"markers_Normal_Crohns_Disease.csv")

## Step 14) Detection of Markers in Cell Types
# Cell Types using Packages (B Cell, T Cell, ... )
require(BiocManager)
install("celldex")
install("SingleR")
install("SingleCellExperiment")
require(celldex)
require(SingleR)
require(SingleCellExperiment)
myref=celldex::HumanPrimaryCellAtlasData(ensembl = F)
dim(myref)
View(myref)
myref$label.main[1:10]
View(as.data.frame(table(myref$label.main)))
View(as.data.frame(table(myref$label.fine)))
mylabels_Normal_Crohns_Disease=SingleR(as.SingleCellExperiment(srobject_Normal_Crohns_Disease), myref, myref$label.main)
length(mylabels_Normal_Crohns_Disease$labels)
mylabels_Normal_Crohns_Disease$labels[1000]
View(as.data.frame(table(mylabels_Normal_Crohns_Disease$labels)))
View(mylabels_Normal_Crohns_Disease$scores)
srobject_Normal_Crohns_Disease$CellType=mylabels_Normal_Crohns_Disease$labels

table(srobject_Normal_Crohns_Disease$seurat_clusters[which(mylabels_Normal_Crohns_Disease$labels=="Cell_Types_Names_Replaced")])


ClustersID=0:number of clusters
CellType=c("cell type 1", "cell type 2", "cell type 3", .......)

names(CellType)=levels(srobject_Normal_Crohns_Disease)
srobject_Normal_Crohns_Disease=RenameIdents(srobject_Normal_Crohns_Disease, CellType)
pdf(file = "UMAP_Name_celltypes.pdf", width = 22, height =13)
DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", label = T, pt.size = 0.5, label.size = 15, raster = FALSE)+ NoLegend() + 
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 70, face = "bold"),
        axis.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40, color="black"),
        axis.text.y = element_text(size=40, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), plot.margin = unit(c(1,1,1,1), "cm")) +
  coord_cartesian(xlim = c(-10.5,15))
dev.off()
# umap_name_split_by_samples
pdf(file = "UMAP_Name_celltypes_for_each_sample.pdf", width = 50, height =20)
DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", label = T, pt.size = 0.5, label.size = 10, raster = FALSE, split.by = "orig.ident", ncol = 6)+ NoLegend() + 
  theme(plot.title = element_text(size = 20),
        text = element_text(size = 40, face = "bold"),
        axis.title = element_text(size=25,face="bold"),
        axis.text.x = element_text(size=25, color="black"),
        axis.text.y = element_text(size=25, color="black"),
        legend.text=element_text(size=20, color="black"),
        legend.title=element_text(size=20, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"))
dev.off()

metadata <- srobject_Normal_Crohns_Disease@meta.data
head(metadata)
freq_table <- as.data.frame(table(metadata$orig.ident, srobject_Normal_Crohns_Disease@active.ident))
colnames(freq_table) <- c("Sample", "CellType", "Frequency")
head(freq_table)
library(ggplot2)
ggplot(freq_table, aes(x = Sample, y = Frequency, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    text = element_text(size = 15)
  ) +
  labs(
    title = "Cell Type Frequency per Sample",
    x = "Samples",
    y = "Frequency",
    fill = "Cell Type"
  )


# Step 16) Plot for abundance of cells in each group
srobject_Normal_Crohns_Disease$Class="Normal"
srobject_Normal_Crohns_Disease$Class[grep("Crohns_Disease", names(srobject_Normal_Crohns_Disease$orig.ident))]="Crohns_Disease"
table(srobject_Normal_Crohns_Disease$Class)
pdf(file = "UMAP_Class.pdf", width = 60, height =25)
p1=DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", group.by = "Class", label.size = 30, pt.size = 0.5, label = T, raster = FALSE) + 
  theme(plot.title = element_text(size = 60),
        text = element_text(size = 130, face = "bold"),
        axis.title = element_text(size=70,face="bold"),
        axis.text.x = element_text(size=70, color="black"),
        axis.text.y = element_text(size=70, color="black"),
        legend.text=element_text(size=50, color="black"),
        legend.title=element_text(size=50, color="black"),
        axis.line.x = element_line(size=3, color="black"),
        axis.line.y = element_line(size=3, color="black"))
p2=DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", label = T, pt.size = 0.5, label.size = 18, raster = FALSE) + NoLegend() + 
  theme(plot.title = element_text(size = 60),
        text = element_text(size = 130, face = "bold"),
        axis.title = element_text(size=70,face="bold"),
        axis.text.x = element_text(size=70, color="black"),
        axis.text.y = element_text(size=70, color="black"),
        legend.text=element_text(size=50, color="black"),
        legend.title=element_text(size=50, color="black"),
        axis.line.x = element_line(size=3, color="black"),
        axis.line.y = element_line(size=3, color="black"))
p1+p2
dev.off()
# Plot Based Groups
pdf(file = "UMAP_GROUP.pdf", width = 60, height =25)
DimPlot(srobject_Normal_Crohns_Disease, reduction = "umap", label = T, split.by = "Class", label.size = 20, pt.size = 0.5, raster = FALSE) + 
  theme(plot.title = element_text(size = 45),
        text = element_text(size = 130, face = "bold"),
        axis.title = element_text(size=70,face="bold"),
        axis.text.x = element_text(size=70, color="black"),
        axis.text.y = element_text(size=70, color="black"),
        legend.text=element_text(size=40, color="black"),
        legend.title=element_text(size=50, color="black"),
        axis.line.x = element_line(size=3, color="black"),
        axis.line.y = element_line(size=3, color="black"))
dev.off()
write.csv(freq_table, "freq_table.csv")