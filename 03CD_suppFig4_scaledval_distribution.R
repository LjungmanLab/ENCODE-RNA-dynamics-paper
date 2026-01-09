#############################################################################################################

## This code was used to generate Main Figures 3C and 3D, as well as Supplemental Figure 4.              ## 
## Input file is Supplemental Table S7                                                                    ##                                               

#############################################################################################################

library(tidyverse)

# read the input file
df.main<-read.delim("Supplemental_Table_S7.xlsx - NormLFCscaled.tsv", sep = "\t", header = T,
               na.strings = c(""," ","NA"), check.names = F)

# Subset the 2h values
df.stab2h <- df.main %>% mutate(gene_ID = paste(gene, `gene name`, sep="/")) %>%
  select(gene_ID, matches("_2v0_ES.LFC.scaled"))

# subset the 6h values
df.stab6h <- df.main %>% mutate(gene_ID = paste(gene, `gene name`, sep="/")) %>%
  select(gene_ID, matches("_6v2_ES.LFC.scaled"))


# check for any differences in the two gene lists
df.stab2h$gene_ID[!(df.stab2h$gene_ID %in% df.stab6h$gene_ID)]


# Convert from wide to long table format
df2h <- df.stab2h %>% pivot_longer (!gene_ID, names_to = 'col_names', values_to = 'value') %>% 
  separate('col_names', c('cellline', 'comparison', 'type'),  sep='_', remove=F)


df6h <- df.stab6h %>% pivot_longer (!gene_ID, names_to = 'col_names', values_to = 'value') %>% 
  separate('col_names', c('cellline', 'comparison', 'type'),  sep='_', remove=F)


# generate individual plots

list_cl = c('A673','Caco2','Calu3', 'GM12878', 'HCT116', 'HepG2', 'HUVEC', 'HMEC', 'IMR90', 'K562', 'MCF7', 'MCF10A', 'OCILY7', 'panc1', 'PC3', 'PC9')
i = 1

for ( i in 1:length(list_cl)){
  cell_line = list_cl[[i]]
  
  cat (cell_line,"\n")
  
  d1 <- df2h[df2h$cellline == cell_line,]
  
  p<-ggplot(d1, aes(x=value)) +
    geom_histogram(fill="#345d9e", color="#345d9e",binwidth=0.01, alpha=1) +
    ggtitle(cell_line) +
    theme_bw(base_size=16)+
    scale_x_continuous(breaks=seq(0, 1, 0.2))+
    xlab("scaled LFC 2V0")+
    theme(axis.text.x=element_text(face="bold")) +
    theme(axis.text.y=element_text(face="bold"))

  
  ggsave(file=paste0("ENCODE16CL_2v0_scaledLFCdist_", cell_line , ".png"), 
         plot = p , 
         width = 6, 
         height = 5, 
         dpi = 300)
  
  ggsave(file=paste0("ENCODE16CL_2v0_scaledLFCdist_", cell_line , ".pdf"), 
         plot = p , 
         width = 6, 
         height = 5, 
         dpi = 300)
  
  #####6v2
  
  d2 <- df6h[df6h$cellline == cell_line,]
  
  
  p1<-ggplot(d2, aes(x=value)) +
    geom_histogram(fill="#b05c38", color="#b05c38", binwidth=0.01, alpha = 1) +
    ggtitle(cell_line) +
    theme_bw(base_size=16)+
    scale_x_continuous(breaks=seq(0, 1, 0.2))+
    xlab("scaled LFC 6v2")+
    theme(axis.text.x=element_text(face="bold")) +
    theme(axis.text.y=element_text(face="bold"))
  
  
  ggsave(file=paste0("ENCODE16CL_6v2_scaledLFCdist_", cell_line , ".png"), 
         plot = p1 , 
         width = 6, 
         height = 5, 
         dpi = 300)
  
  ggsave(file=paste0("ENCODE16CL_6v2_scaledLFCdist_", cell_line , ".pdf"), 
         plot = p1 , 
         width = 6, 
         height = 5, 
         dpi = 300)
  
}


## generate a plot containing all the cell lines 

# 2v0 scaled stability values

p2<-ggplot(df2h, aes(x=value)) +
  geom_histogram(fill="#345d9e", color="#345d9e",binwidth=0.01, alpha=1) +
  theme_bw(base_size=24)+
  scale_x_continuous(breaks=seq(0, 1, 0.2))+
  scale_y_continuous(breaks=seq(0, 200, 50))+
  xlab("\nscaled LFC 2V0")+
  theme(axis.text.x=element_text(face="bold", size=14)) +
  theme(axis.text.y=element_text(face="bold"))+
  facet_wrap(~cellline, nrow=2)+
  theme(strip.background = element_rect(fill = NA, color = "#FFFFFF"),
    strip.text = element_text(size=22, face="bold"),
    panel.spacing = unit(-.01,"cm"))


ggsave("ENCODE16CL_2v0_scaledLFCdist_allcelllines.png", 
       plot = p2 , 
       width = 25, 
       height = 8, 
       dpi = 300)


ggsave("ENCODE16CL_2v0_scaledLFCdist_allcelllines.pdf", 
       plot = p2 , 
       width = 25, 
       height = 8, 
       dpi = 300)




# 6v2 scaled stability values

p3<-ggplot(df6h, aes(x=value)) +
  geom_histogram(fill="#b05c38", color="#b05c38",binwidth=0.01, alpha=1) +
  theme_bw(base_size=24)+
  scale_x_continuous(breaks=seq(0, 1, 0.2))+
  scale_y_continuous(breaks=seq(0, 280, 50))+
  xlab("\nscaled LFC 6v2")+
  theme(axis.text.x=element_text(face="bold", size=14)) +
  theme(axis.text.y=element_text(face="bold"))+
  facet_wrap(~cellline, nrow=2)+
  theme(strip.background = element_rect(fill = NA, color = "#FFFFFF"),
    strip.text = element_text(size=22, face="bold"),
    panel.spacing = unit(-.01,"cm"))


ggsave("ENCODE16CL_6v2_scaledLFCdist_allcelllines.png", 
       plot = p3 , 
       width = 25, 
       height = 8, 
       dpi = 300)


ggsave("ENCODE16CL_6v2_scaledLFCdist_allcelllines.pdf", 
       plot = p3 , 
       width = 25, 
       height = 8, 
       dpi = 300)



