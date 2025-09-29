library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(ggplot2)

# Plotting only the common introns that meet the jrsum=10 cutoff (for all introns)
# in all 3 time points

setwd("./")
wd = getwd()

outdir = wd

# Color pallete
my.colors = read.delim("./ENCODE_colors_16CellLines.txt", header = T, comment.char = "")
head(my.colors)

# Common introns across the 3 time points
common_introns = "24272"

bruseq.common = read.delim("intersection_16MRCL_",common_introns,"introns_jrsum10_MR_ALL_ISI_bruseq.bed",
                           sep = "\t", header = T, na.strings = c(""," ","Inf"),check.names = F)
bruchase2h.common = read.delim("intersection_16MRCL_",common_introns,"introns_jrsum10_MR_ALL_ISI_bruchase2h.bed",
                           sep = "\t", header = T, na.strings = c(""," ","Inf"),check.names = F)
bruchase6h.common = read.delim("intersection_16MRCL_",common_introns,"introns_jrsum10_MR_ALL_ISI_bruchase6h.bed",
                               sep = "\t", header = T, na.strings = c(""," ","Inf"),check.names = F)

df_list = list(bruseq.common,bruchase2h.common,bruchase6h.common)

# Combine and melt the 3 df's
merged.df = df_list %>% reduce(full_join, by = "gene_ensg_featureID")
merged.df.melt = reshape2::melt(merged.df, value.name = "SI")

# Create new variables columns
merged.df.melt$cell_line = gsub("_SI_[026]h","",merged.df.melt$variable)
merged.df.melt$cell_line = gsub("-","",merged.df.melt$cell_line)
merged.df.melt = merged.df.melt %>%
  mutate(time_point = case_when(
    str_detect(variable, "_0h") ~ "0h",
    str_detect(variable, "_2h") ~ "2h",
    str_detect(variable, "_6h") ~ "6h"
  )) %>%
  as.data.frame()

# Plotting each time point separately
time_points = c("0h","2h","6h")
# my_time_point = "0h"

for (my_time_point in time_points) {

  merged.df.melt.tp = merged.df.melt %>%
    filter(time_point == my_time_point) %>%
    as.data.frame()

  if (my_time_point=="0h") {
    plot.title = paste0("n = ",common_introns," introns at Bruseq 0h")
  } else if(my_time_point=="2h"){
    plot.title = paste0("n = ",common_introns," introns at BruChaseseq 2h")
  } else{
    plot.title = paste0("n = ",common_introns," introns at BruChaseseq 6h")
  }

  CL_order = c("IMR90","OCILY7","GM12878","K562","MCF7","A673","Calu3","HepG2",
               "Caco2","PC3","HUVEC","HMEC","MCF10A","HCT116","PC9","panc1")

  merged.df.melt.tp$cell_line = factor(merged.df.melt.tp$cell_line, levels = CL_order)
  my.colors.reordered = my.colors[match(CL_order,my.colors$cell_line),]

  # Plot each TP separately as boxplots
  myplot = ggplot(merged.df.melt.tp, aes(x = cell_line, y = SI, fill = cell_line)) +
    geom_boxplot(outlier.size=0.5, outlier.alpha = 0.1) +
    scale_fill_manual(values = my.colors.reordered$color) +
    labs(x = "Cell Lines", y = "Splicing Index",title = plot.title) +
    theme_bw() +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 8,face = "bold"),
          axis.title = element_text(size = 12,face = "bold"))
  myplot
  plot.filename = paste0(outdir,"intersection_16MRCL_",common_introns,"introns_jrsum10_MR_ALL_ISI_boxplots_REORDERED_",my_time_point)
  ggsave(filename = paste0(plot.filename,".jpg"),plot = myplot, dpi = 300)
  ggsave(filename = paste0(plot.filename,".pdf"),plot = myplot, dpi = 300)


  # Plot each TP separately as violin plots
  myplot = ggplot(merged.df.melt.tp, aes(x = cell_line, y = SI, fill = cell_line)) +
    # geom_boxplot(outlier.size=0.5, outlier.alpha = 0.1) +
    geom_violin(trim = T) +
    scale_fill_manual(values = my.colors.reordered$color) +
    labs(x = "Cell Lines", y = "Splicing Index",title = plot.title) +
    theme_bw() +
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 8,face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12,face = "bold"))
  myplot= myplot + geom_boxplot(width=0.1, lwd=0.2,outlier.size = 0.1, outlier.shape = 4,fill="white")
  myplot
  plot.filename = paste0(outdir,"intersection_16MRCL_",common_introns,"introns_jrsum10_MR_ALL_ISI_violin&boxplots_REORDERED_",my_time_point)
  ggsave(filename = paste0(plot.filename,".jpg"),plot = myplot, dpi = 300, width = 6, height = 4, units = "in")
  ggsave(filename = paste0(plot.filename,".pdf"),plot = myplot, dpi = 300, width = 6, height = 4, units = "in")

}





