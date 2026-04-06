library(tidyverse)
library(readxl)
library(patchwork)

file_path <- "CD.xlsx"
sheet_name <- "COAD"

data <- read_excel(file_path, sheet = sheet_name) %>%
  mutate(across(-Samples, ~ as.numeric(gsub("/", ".", .))))

genes <- colnames(data)[2:11]
colors <- c("Normal" = "#000000", "Tumor" = "#98B500")

plot_list <- list()

for (gene in genes) {
  df_gene <- data %>%
    select(Samples, !!sym(gene)) %>%
    rename(Expression = !!sym(gene)) %>%
    mutate(Samples = factor(Samples, levels = c("Normal", "Tumor")))
  
  pval <- t.test(Expression ~ Samples, data = df_gene)$p.value
  p_label <- paste0("p = ", format(pval, digits = 3, scientific = TRUE))
  
  p <- ggplot(df_gene, aes(x = Samples, y = Expression, fill = Samples)) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5, color = "black") +
    geom_jitter(aes(color = Samples), width = 0.1, size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_minimal(base_family = "sans") +
    labs(
      title = gene,
      x = NULL,
      y = "Expression level (TCGA-COAD)"
    ) +
    annotate(
      "text",
      x = 1.5,
      y = max(df_gene$Expression, na.rm = TRUE) * 1.07,
      label = p_label,
      fontface = "bold",
      size = 8
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 20),
      axis.text.y = element_text(face = "bold", size = 22, color = "black"),  
      axis.text.x = element_text(face = "bold", size = 21, color = "black")  
    )
  
  plot_list[[gene]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 10)

ggsave("Expression_COAD_top10_genes_row_pretty2.pdf",
       combined_plot,
       width = 3.2 * length(genes),
       height = 5,
       dpi = 1200)
