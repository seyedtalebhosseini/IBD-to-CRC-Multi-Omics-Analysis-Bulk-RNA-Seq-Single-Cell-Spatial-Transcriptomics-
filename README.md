# IBD to CRC Multi Omics Analysis (Bulk RNA Seq, Single Cell, Spatial Transcriptomics)
This repository contains analysis scripts used in the study:

Deciphering Immune and Cellular Reprogramming During the Progression from Inflammatory Bowel Disease to Colorectal Cancer Using Multi-Omics Single-Cell and Spatial Transcriptomics

Accepted for publication in Journal of Translational Medicine.

Citation details and DOI will be added once the article is published.

Project Overview:
This project investigates immune and cellular reprogramming during the progression from Inflammatory Bowel Disease (IBD) to Colorectal Cancer (CRC) using an integrative multi-omics framework combining:
Bulk RNA sequencing
Single-cell RNA sequencing
Spatial transcriptomics
TCGA transcriptomic datasets

The repository includes scripts used for data preprocessing, differential expression analysis, Single cell analysis, trajectory inference, spatial analysis, cell-cell communication analysis, and pathway enrichment.

Repository Structure:
Each directory corresponds to a specific analysis step used in the study.

AUC:
Scripts for calculating Area Under the Curve (AUC) scores to evaluate gene signatures and predictive performance.

Bulk-RNASeq Analysis in Ubuntu:
Pipeline and scripts for preprocessing and analyzing bulk RNA‑seq data in a Linux/Ubuntu environment.

CellChat:
Analysis of cell-cell communication networks using the CellChat framework to infer signaling interactions between cell populations.

Clusterprofiler:
Functional enrichment analyses including:

Gene Ontology (GO)
KEGG pathway enrichment
pathway visualization
performed using the clusterProfiler R package.

DESeq2 and SVA:
Differential gene expression analysis using DESeq2, with batch effect correction using SVA (Surrogate Variable Analysis).

Gene Expression Spatial:
Scripts for analyzing spatial gene expression patterns derived from spatial transcriptomics datasets.

GraphST For Spatial Distribution of cell types:
Implementation of GraphST for identifying and visualizing the spatial distribution of cell types within tissue sections.

Monocle Trajectories:
Trajectory inference analysis using Monocle2 to model dynamic cellular transitions and lineage relationships.

Moran’s I:
Spatial autocorrelation analysis using Moran’s I to identify spatially variable genes.

Pseudotime Trajectories Plots:
Visualization scripts for pseudotime trajectories, lineage progression, and gene expression dynamics.

RDSTOH5AD:
Utility scripts for converting RDS files to H5AD format, enabling compatibility between Seurat and Scanpy workflows.

SingleCell:
Core pipeline for single-cell RNA-seq analysis, including preprocessing, normalization, clustering, and cell type annotation.

Spatial Decomposition:
Methods for spatial deconvolution to estimate cell type composition across spatial transcriptomics spots.

TCGAbiolinks:
Scripts using TCGAbiolinks to download and process TCGA colorectal cancer datasets.

TCGA Gene Expression:
Analysis of TCGA gene expression profiles for validation and comparison with single‑cell and spatial transcriptomic findings.

Computational Environment
All analyses were performed using:

Ubuntu 22.04 LTS
R 4.3.2
Python

Software and Packages
Data Processing and Bioinformatics Tools:
SRAToolkits
FASTQC
TRIMMOMATIC
HISAT2
HTSeq

R Packages:
SVA
DESeq2
clusterProfiler
Seurat
Celldex
SingleCellExperiment
SingleR
Monocle2
CellChat
SeuratDisk
SeuratData
SpacexR
dplyr
spdep
TCGAbiolinks
pROC
survival
survminer

Python Packages:
Scanpy
Matplotlib
OS
Torch
Pandas
NumPy
Skmisc
Scikit-learn
OT
GraphST


Citation
If you use the scripts in this repository, please cite:

Deciphering Immune and Cellular Reprogramming During the Progression from Inflammatory Bowel Disease to Colorectal Cancer Using Multi‑Omics Single‑Cell and Spatial Transcriptomics

Journal of Translational Medicine (Accepted)

License
This repository is intended for academic and research purposes.

Contact
For questions or collaboration inquiries, please contact the repository author.
