library(pROC)
library(ggplot2)
library(dplyr)
library(readxl)
library(gridExtra)
library(grid)
file_path <- "BIO/Selected_Genes_TCGA_COAD_READ.xlsx"
coad_data <- read_excel(file_path, sheet = "COAD-NEW")
read_data <- read_excel(file_path, sheet = "READ-NEW")
crohn_genes <- c("CEACAM5", "CENPE", "NUAK1", "CALCRL", "ADAMTS4", 
                 "CHI3L1", "SERPINE1", "ABL2", "MAPK10", "MMP3")

ulcerative_genes <- c("CEACAM5", "MMP12", "NAMPT", "LGALS1", "PDGFRA")
gene_colors <- c(
  "#B40426", "#DD4F22", "#FCA50A", "#FDE725", "#C4E338",
  "#90D84C", "#63CB7E", "#78B7EB", "#6DA2E2", "#5278D3"
)
plot_roc_with_colored_cells <- function(df, genes, title){
  
  labels <- ifelse(grepl("Tumor", df[[1]], ignore.case = TRUE), 1, 0)
  genes <- genes[genes %in% colnames(df)]
  
  
  roc_list <- lapply(genes, function(g){ roc(labels, df[[g]]) })
  names(roc_list) <- genes
  
  
  auc_values <- sapply(roc_list, function(x) round(auc(x),3))
  auc_df <- data.frame(Gene = names(auc_values), AUC = auc_values)
  auc_df <- auc_df %>% arrange(desc(AUC))
  
  
  sorted_genes <- auc_df$Gene
  
  
  color_map <- gene_colors[1:length(sorted_genes)]
  names(color_map) <- sorted_genes
  
  
  
  roc_df <- bind_rows(lapply(sorted_genes, function(g){
    
    r <- roc_list[[g]]
    
    
    FPR_raw <- 1 - r$specificities
    TPR_raw <- r$sensitivities
    
    
    smooth_x <- seq(min(FPR_raw), max(FPR_raw), length.out = 300)
    smooth_y <- spline(FPR_raw, TPR_raw, xout = smooth_x)$y
    
    data.frame(FPR = smooth_x, TPR = smooth_y, Gene = g)
  }))
  
  
  ttheme_colored_cells <- ttheme_default(
    core = list(
      fg_params=list(cex=1, col="black", fontface="bold"),
      bg_params=list(
        fill=color_map,  
        col=NA
      )
    ),
    colhead = list(fg_params=list(cex=1.2, fontface="bold"))
  )
  
  auc_table <- tableGrob(auc_df, rows=NULL, theme=ttheme_colored_cells)
  
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR, color = Gene)) +
    geom_line(size=1.4) +
    geom_abline(intercept=0, slope=1, size=1.3, color="black") +
    scale_color_manual(values=color_map) +
    theme_minimal() +
    labs(title=title, x="100% - Specificity", y="Sensitivity %") +
    theme(
      plot.title = element_text(size=16, face="bold",color="black"),
      axis.title = element_text(size=16, face="bold",color="black"),
      axis.text = element_text(size=16, face="bold",color="black"),
      legend.position = "none"
    )
  
  arrangeGrob(p, auc_table, ncol=2, widths=c(3,1))
}

gA <- plot_roc_with_colored_cells(coad_data, crohn_genes, "Crohn's Disease Biomarkers (TCGA-COAD)")
gB <- plot_roc_with_colored_cells(read_data, crohn_genes, "Crohn's Disease Biomarkers (TCGA-READ)")
gC <- plot_roc_with_colored_cells(coad_data, ulcerative_genes, "Ulcerative Colitis Biomarkers (TCGA-COAD)")
gD <- plot_roc_with_colored_cells(read_data, ulcerative_genes, "Ulcerative Colitis Biomarkers (TCGA-READ)")
pdf("ROC_Figures_Crohn.pdf", width=13, height=4)
grid.arrange(gA, gB, nrow=1)
dev.off()
pdf("ROC_Figures_Ulcerative.pdf", width=13, height=4)
grid.arrange(gC, gD, nrow=1)
dev.off()