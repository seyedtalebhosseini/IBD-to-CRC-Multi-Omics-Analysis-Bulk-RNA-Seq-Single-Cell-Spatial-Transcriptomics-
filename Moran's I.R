library(Seurat)
library(dplyr)
library(spdep)
library(openxlsx)

csv_files <- c(
  "data/tissue-position/Human Colorectal Cancer - 11 mm Capture Area (FFPE).csv",
  "data/tissue-position/Human Colorectal Cancer - Whole Transcriptome Analysis.csv",
  "data/tissue-position/Human Intestine Cancer (FPPE).csv"
)

genes <- c("CEACAM5", "MMP12", "NAMPT", "LGALS1", "PDGFRA")

obj <- LoadH5Seurat("data/srobject_UC_H5Seurat/srobject_B_Cell.h5Seurat")

expr_mat <- GetAssayData(obj, assay = "RNA", slot = "data")
genes_available <- genes[genes %in% rownames(expr_mat)]

expr_target <- t(as.matrix(expr_mat[genes_available, , drop = FALSE]))

base_df <- data.frame(
  first_type = as.character(obj@active.ident),
  expr_target,
  stringsAsFactors = FALSE
)

N <- nrow(base_df)

for (csv_path in csv_files) {
  
  cat("Processing:", csv_path, "\n")
  
  
  excel_name <- paste0(tools::file_path_sans_ext(csv_path), "_MoranI.xlsx")
  
  
  csv_data <- read.csv(csv_path)
  M <- nrow(csv_data)
  
  
  x_all <- csv_data[[ncol(csv_data)-1]]
  y_all <- csv_data[[ncol(csv_data)]]
  
  
  if (M >= N) {
    x_coords <- x_all[1:N]
    y_coords <- y_all[1:N]
  } else {
    x_coords <- c(x_all, rep(NA, N - M))
    y_coords <- c(y_all, rep(NA, N - M))
  }
  
  df <- base_df
  df$x <- x_coords
  df$y <- y_coords
  
  
  df <- df[complete.cases(df[, c("x","y")]), ]
  
  cat(" → Final usable cells:", nrow(df), "\n")
  
  
  
  genes2 <- colnames(df)[!(colnames(df) %in% c("first_type","x","y"))]
  cell_types <- unique(df$first_type)
  
  Moran_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(genes2),
                         dimnames = list(cell_types, genes2))
  
  Pval_matrix <- Moran_matrix
  
  for (ct in cell_types) {
    for (gene in genes2) {
      sub_df <- df[df$first_type == ct, ]
      
      if (nrow(sub_df) < 4) next
      
      values <- sub_df[[gene]]
      coords <- sub_df[, c("x","y")]
      
      nb <- knn2nb(knearneigh(coords, k = 4))
      lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
      
      mi <- try(moran.test(values, lw, zero.policy = TRUE), silent = TRUE)
      
      if (!inherits(mi,"try-error")) {
        Moran_matrix[ct, gene] <- mi$estimate["Moran I statistic"]
        Pval_matrix[ct, gene]  <- mi$p.value
      }
    }
  }
  
  
  pvec <- as.vector(Pval_matrix)
  fdr_vec <- p.adjust(pvec, method = "BH")
  
  FDR_matrix <- matrix(
    fdr_vec,
    nrow = nrow(Pval_matrix),
    ncol = ncol(Pval_matrix),
    dimnames = dimnames(Pval_matrix)
  )
  
  
  wb <- createWorkbook()
  
  # Sheet1 - Moran’s I
  addWorksheet(wb, "Moran_I")
  writeData(wb, "Moran_I", data.frame(CellType = rownames(Moran_matrix), Moran_matrix))
  
  # Sheet2 - P-values
  addWorksheet(wb, "P_values")
  
  writeData(wb, "P_values", data.frame(CellType = rownames(Pval_matrix), Pval_matrix))
  
  # Sheet3 - FDR
  addWorksheet(wb, "FDR")
  writeData(wb, "FDR", data.frame(CellType = rownames(FDR_matrix), FDR_matrix))
  
  
  saveWorkbook(wb, excel_name, overwrite = TRUE)
  
  cat("Saved:", excel_name, "\n\n")
}

