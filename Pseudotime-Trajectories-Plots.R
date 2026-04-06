library(Seurat)
library(monocle)
library(dplyr)
library(ggplot2)
library(patchwork)

#CD Objects
srobject_trajectori_CD <- readRDS("CD/srobject_trajectori_CD.RDS")
T_cell_trajectori_CD <- readRDS("CD/T_cell_trajectori_CD.RDS")
B_cell_trajectori_CD <- readRDS("CD/B_cell_trajectori_CD.RDS")
#UC Objects
srobject_trajectori_UC <- readRDS("UC/srobject_trajectori_UC.RDS")
T_cell_trajectori_UC <- readRDS("UC/T_cell_trajectori_UC.RDS")
B_cell_trajectori_UC <- readRDS("UC/B_cell_trajectori_UC.RDS")

trajectory_colors <- c(
  "#B40426",  
  "#DD4F22",  
  "#FCA50A",  
  "#FDE725",  
  "#C4E338",  
  "#90D84C",  
  "#63CB7E",  
  "#78B7EB",  
  "#6DA2E2",  
  "#5278D3",  
  "#3B4CC0"   
)

#CD
p1=plot_cell_trajectory(srobject_trajectori_CD, color_by = "active.ident", cell_size = 7) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 40),
        text = element_text(size = 40, face = "bold"),
        axis.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40, color="black"),
        axis.text.y = element_text(size=40, color="black"),
        legend.text=element_text(size=40, color="black"),
        legend.title=element_text(size=40, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Cell Density Distribution")

p2=plot_cell_trajectory(srobject_trajectori_CD, color_by = "Pseudotime", cell_size = 7) +
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 40),
        text = element_text(size = 40, face = "bold"),
        axis.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40, color="black"),
        axis.text.y = element_text(size=40, color="black"),
        legend.text=element_text(size=40, color="black"),
        legend.title=element_text(size=40, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")

p3=plot_cell_trajectory(srobject_trajectori_CD, color_by = "Class", cell_size = 7) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 40),
        text = element_text(size = 40, face = "bold"),
        axis.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40, color="black"),
        axis.text.y = element_text(size=40, color="black"),
        legend.text=element_text(size=40, color="black"),
        legend.title=element_text(size=40, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")

p4=plot_cell_trajectory(srobject_trajectori_CD, color_by = "State", cell_size = 7) + 
  theme(legend.position = "right") +
  theme(plot.title = element_text(size = 40),
        text = element_text(size = 40, face = "bold"),
        axis.title = element_text(size=40,face="bold"),
        axis.text.x = element_text(size=40, color="black"),
        axis.text.y = element_text(size=40, color="black"),
        legend.text=element_text(size=40, color="black"),
        legend.title=element_text(size=40, color="black"),
        axis.line.x = element_line(size=1, color="black"),
        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")

p5=plot_cell_trajectory(T_cell_trajectori_CD, color_by = "active.ident", cell_size = 7) + 
    theme(legend.position = "right") +
      theme(plot.title = element_text(size = 40),
            text = element_text(size = 40, face = "bold"),
            axis.title = element_text(size=40,face="bold"),
            axis.text.x = element_text(size=40, color="black"),
            axis.text.y = element_text(size=40, color="black"),
            legend.text=element_text(size=40, color="black"),
            legend.title=element_text(size=40, color="black"),
            axis.line.x = element_line(size=1, color="black"),
            axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "T Cell - Subclusters")
    
p6=plot_cell_trajectory(T_cell_trajectori_CD, color_by = "Pseudotime", cell_size = 7) + 
      theme(legend.position = "right") +
        theme(plot.title = element_text(size = 40),
              text = element_text(size = 40, face = "bold"),
              axis.title = element_text(size=40,face="bold"),
              axis.text.x = element_text(size=40, color="black"),
              axis.text.y = element_text(size=40, color="black"),
              legend.text=element_text(size=40, color="black"),
              legend.title=element_text(size=40, color="black"),
              axis.line.x = element_line(size=1, color="black"),
              axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")
      
p7=plot_cell_trajectory(T_cell_trajectori_CD, color_by = "Class", cell_size = 7) + 
        theme(legend.position = "right") +
          theme(plot.title = element_text(size = 40),
                text = element_text(size = 40, face = "bold"),
                axis.title = element_text(size=40,face="bold"),
                axis.text.x = element_text(size=40, color="black"),
                axis.text.y = element_text(size=40, color="black"),
                legend.text=element_text(size=40, color="black"),
                legend.title=element_text(size=40, color="black"),
                axis.line.x = element_line(size=1, color="black"),
                axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")
        
p8=plot_cell_trajectory(T_cell_trajectori_CD, color_by = "State", cell_size = 7) + 
          theme(legend.position = "right") +
            theme(plot.title = element_text(size = 40),
                  text = element_text(size = 40, face = "bold"),
                  axis.title = element_text(size=40,face="bold"),
                  axis.text.x = element_text(size=40, color="black"),
                  axis.text.y = element_text(size=40, color="black"),
                  legend.text=element_text(size=40, color="black"),
                  legend.title=element_text(size=40, color="black"),
                  axis.line.x = element_line(size=1, color="black"),
                  axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")

p9=plot_cell_trajectory(B_cell_trajectori_CD, color_by = "active.ident", cell_size = 7) + 
            theme(legend.position = "right") +
              theme(plot.title = element_text(size = 40),
                    text = element_text(size = 40, face = "bold"),
                    axis.title = element_text(size=40,face="bold"),
                    axis.text.x = element_text(size=40, color="black"),
                    axis.text.y = element_text(size=40, color="black"),
                    legend.text=element_text(size=40, color="black"),
                    legend.title=element_text(size=40, color="black"),
                    axis.line.x = element_line(size=1, color="black"),
                    axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "B Cell - Subclusters")
            
p10=plot_cell_trajectory(B_cell_trajectori_CD, color_by = "Pseudotime", cell_size = 7) + 
              theme(legend.position = "right") +
                theme(plot.title = element_text(size = 40),
                      text = element_text(size = 40, face = "bold"),
                      axis.title = element_text(size=40,face="bold"),
                      axis.text.x = element_text(size=40, color="black"),
                      axis.text.y = element_text(size=40, color="black"),
                      legend.text=element_text(size=40, color="black"),
                      legend.title=element_text(size=40, color="black"),
                      axis.line.x = element_line(size=1, color="black"),
                      axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")
              
p11=plot_cell_trajectory(B_cell_trajectori_CD, color_by = "Class", cell_size = 7) + 
                theme(legend.position = "right") +
                  theme(plot.title = element_text(size = 40),
                        text = element_text(size = 40, face = "bold"),
                        axis.title = element_text(size=40,face="bold"),
                        axis.text.x = element_text(size=40, color="black"),
                        axis.text.y = element_text(size=40, color="black"),
                        legend.text=element_text(size=40, color="black"),
                        legend.title=element_text(size=40, color="black"),
                        axis.line.x = element_line(size=1, color="black"),
                        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")
                
p12=plot_cell_trajectory(B_cell_trajectori_CD, color_by = "State", cell_size = 7) + 
                  theme(legend.position = "right") +
                    theme(plot.title = element_text(size = 40),
                          text = element_text(size = 40, face = "bold"),
                          axis.title = element_text(size=40,face="bold"),
                          axis.text.x = element_text(size=40, color="black"),
                          axis.text.y = element_text(size=40, color="black"),
                          legend.text=element_text(size=40, color="black"),
                          legend.title=element_text(size=40, color="black"),
                          axis.line.x = element_line(size=1, color="black"),
                          axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")
                  
# UC
p13=plot_cell_trajectory(srobject_trajectori_UC, color_by = "active.ident", cell_size = 7) + 
            theme(legend.position = "right") +
              theme(plot.title = element_text(size = 40),
                    text = element_text(size = 40, face = "bold"),
                    axis.title = element_text(size=40,face="bold"),
                    axis.text.x = element_text(size=40, color="black"),
                    axis.text.y = element_text(size=40, color="black"),
                    legend.text=element_text(size=40, color="black"),
                    legend.title=element_text(size=40, color="black"),
                    axis.line.x = element_line(size=1, color="black"),
                    axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Cell Density Distribution")
            
p14=plot_cell_trajectory(srobject_trajectori_UC, color_by = "Pseudotime", cell_size = 7) + 
              theme(legend.position = "right") +
                theme(plot.title = element_text(size = 40),
                      text = element_text(size = 40, face = "bold"),
                      axis.title = element_text(size=40,face="bold"),
                      axis.text.x = element_text(size=40, color="black"),
                      axis.text.y = element_text(size=40, color="black"),
                      legend.text=element_text(size=40, color="black"),
                      legend.title=element_text(size=40, color="black"),
                      axis.line.x = element_line(size=1, color="black"),
                      axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")
              
p15=plot_cell_trajectory(srobject_trajectori_UC, color_by = "Class", cell_size = 7) + 
                theme(legend.position = "right") +
                  theme(plot.title = element_text(size = 40),
                        text = element_text(size = 40, face = "bold"),
                        axis.title = element_text(size=40,face="bold"),
                        axis.text.x = element_text(size=40, color="black"),
                        axis.text.y = element_text(size=40, color="black"),
                        legend.text=element_text(size=40, color="black"),
                        legend.title=element_text(size=40, color="black"),
                        axis.line.x = element_line(size=1, color="black"),
                        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")
                
p16=plot_cell_trajectory(srobject_trajectori_UC, color_by = "State", cell_size = 7) + 
                  theme(legend.position = "right") +
                    theme(plot.title = element_text(size = 40),
                          text = element_text(size = 40, face = "bold"),
                          axis.title = element_text(size=40,face="bold"),
                          axis.text.x = element_text(size=40, color="black"),
                          axis.text.y = element_text(size=40, color="black"),
                          legend.text=element_text(size=40, color="black"),
                          legend.title=element_text(size=40, color="black"),
                          axis.line.x = element_line(size=1, color="black"),
                          axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")
                  
p17=plot_cell_trajectory(T_cell_trajectori_UC, color_by = "active.ident", cell_size = 7) + 
                    theme(legend.position = "right") +
                      theme(plot.title = element_text(size = 40),
                            text = element_text(size = 40, face = "bold"),
                            axis.title = element_text(size=40,face="bold"),
                            axis.text.x = element_text(size=40, color="black"),
                            axis.text.y = element_text(size=40, color="black"),
                            legend.text=element_text(size=40, color="black"),
                            legend.title=element_text(size=40, color="black"),
                            axis.line.x = element_line(size=1, color="black"),
                            axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "T Cell - Subclusters")
                    
p18=plot_cell_trajectory(T_cell_trajectori_UC, color_by = "Pseudotime", cell_size = 7) + 
                      theme(legend.position = "right") +
                        theme(plot.title = element_text(size = 40),
                              text = element_text(size = 40, face = "bold"),
                              axis.title = element_text(size=40,face="bold"),
                              axis.text.x = element_text(size=40, color="black"),
                              axis.text.y = element_text(size=40, color="black"),
                              legend.text=element_text(size=40, color="black"),
                              legend.title=element_text(size=40, color="black"),
                              axis.line.x = element_line(size=1, color="black"),
                              axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")
                      
p19=plot_cell_trajectory(T_cell_trajectori_UC, color_by = "Class", cell_size = 7) + 
                        theme(legend.position = "right") +
                          theme(plot.title = element_text(size = 40),
                                text = element_text(size = 40, face = "bold"),
                                axis.title = element_text(size=40,face="bold"),
                                axis.text.x = element_text(size=40, color="black"),
                                axis.text.y = element_text(size=40, color="black"),
                                legend.text=element_text(size=40, color="black"),
                                legend.title=element_text(size=40, color="black"),
                                axis.line.x = element_line(size=1, color="black"),
                                axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")
                        
p20=plot_cell_trajectory(T_cell_trajectori_UC, color_by = "State", cell_size = 7) + 
                          theme(legend.position = "right") +
                            theme(plot.title = element_text(size = 40),
                                  text = element_text(size = 40, face = "bold"),
                                  axis.title = element_text(size=40,face="bold"),
                                  axis.text.x = element_text(size=40, color="black"),
                                  axis.text.y = element_text(size=40, color="black"),
                                  legend.text=element_text(size=40, color="black"),
                                  legend.title=element_text(size=40, color="black"),
                                  axis.line.x = element_line(size=1, color="black"),
                                  axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")
                          
p21=plot_cell_trajectory(B_cell_trajectori_UC, color_by = "active.ident", cell_size = 7) + 
                            theme(legend.position = "right") +
                              theme(plot.title = element_text(size = 40),
                                    text = element_text(size = 40, face = "bold"),
                                    axis.title = element_text(size=40,face="bold"),
                                    axis.text.x = element_text(size=40, color="black"),
                                    axis.text.y = element_text(size=40, color="black"),
                                    legend.text=element_text(size=40, color="black"),
                                    legend.title=element_text(size=40, color="black"),
                                    axis.line.x = element_line(size=1, color="black"),
                                    axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "B Cell - Subclusters")
                            
p22=plot_cell_trajectory(B_cell_trajectori_UC, color_by = "Pseudotime", cell_size = 7) + 
                              theme(legend.position = "right") +
                                theme(plot.title = element_text(size = 40),
                                      text = element_text(size = 40, face = "bold"),
                                      axis.title = element_text(size=40,face="bold"),
                                      axis.text.x = element_text(size=40, color="black"),
                                      axis.text.y = element_text(size=40, color="black"),
                                      legend.text=element_text(size=40, color="black"),
                                      legend.title=element_text(size=40, color="black"),
                                      axis.line.x = element_line(size=1, color="black"),
                                      axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Pseudotime")
                              
p23=plot_cell_trajectory(B_cell_trajectori_UC, color_by = "Class", cell_size = 7) + 
                                theme(legend.position = "right") +
                                  theme(plot.title = element_text(size = 40),
                                        text = element_text(size = 40, face = "bold"),
                                        axis.title = element_text(size=40,face="bold"),
                                        axis.text.x = element_text(size=40, color="black"),
                                        axis.text.y = element_text(size=40, color="black"),
                                        legend.text=element_text(size=40, color="black"),
                                        legend.title=element_text(size=40, color="black"),
                                        axis.line.x = element_line(size=1, color="black"),
                                        axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "Class")
                                
p24=plot_cell_trajectory(B_cell_trajectori_UC, color_by = "State", cell_size = 7) + 
                                  theme(legend.position = "right") +
                                    theme(plot.title = element_text(size = 40),
                                          text = element_text(size = 40, face = "bold"),
                                          axis.title = element_text(size=40,face="bold"),
                                          axis.text.x = element_text(size=40, color="black"),
                                          axis.text.y = element_text(size=40, color="black"),
                                          legend.text=element_text(size=40, color="black"),
                                          legend.title=element_text(size=40, color="black"),
                                          axis.line.x = element_line(size=1, color="black"),
                                          axis.line.y = element_line(size=1, color="black"), legend.position = "right") + labs(title = "State")



all_plots <- list(
  p1, p2, p3, p4,
  p5, p6, p7, p8,
  p9, p10, p11, p12,
  p13, p14, p15, p16,
  p17, p18, p19, p20,
  p21, p22, p23, p24
)




label_positions <- c(1, 5, 9, 13, 17, 21)


tags <- paste0("(", LETTERS[1:length(label_positions)],")")   # (A) ... (F)


all_plots_tagged <- all_plots

for (i in seq_along(label_positions)) {
  
  idx <- label_positions[i]   
  tag <- tags[i]              
  
  all_plots_tagged[[idx]] <- all_plots[[idx]] +
    labs(tag = tag) +
    theme(
      plot.tag = element_text(
        size = 60,
        face = "bold",
        hjust = 0,
        vjust = 1
      ),
      plot.tag.position = c(0.005, 0.98)
    )
}


pdf("Figure4.pdf", width = 50, height = 45)
wrap_plots(all_plots_tagged, ncol = 4)
dev.off()

