##########################################################################################
##  Retrieve stability values from RDS object and plot correlation heatmaps.            ##
##  Used in generating Figure 2A and B in Bedi et al 2024.                              ##
##  Written by Brian Magnuson (bmagnuso@umich.edu).                                     ##
##  Correspondences to Mats Ljungman (ljungman@umich.edu).                              ##
##########################################################################################

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

# load exonic stability
ex.20 <- readRDS("ENCODE4.16CellLines.Stability.2v0.ExonSense.rds")
ex.60 <- readRDS("ENCODE4.16CellLines.Stability.6v0.ExonSense.rds")
ex.62 <- readRDS("ENCODE4.16CellLines.Stability.6v2.ExonSense.rds")

# load genic stability
gn.20 <- readRDS("ENCODE4.16CellLines.Stability.2v0.GeneSense.rds")
gn.60 <- readRDS("ENCODE4.16CellLines.Stability.6v0.GeneSense.rds")
gn.62 <- readRDS("ENCODE4.16CellLines.Stability.6v2.GeneSense.rds")

# load sample metadata
metadata <- read.delim("ENCODE4.16CellLines.BruChase.Sample.Metadata.txt",sep="\t",header=TRUE)
# load gene annotation
annot <- read.delim("hg38.gencode_29.ensembl_gene.transcript_biotype_name_entrez.gencode.txt",sep="\t",header=TRUE)


### process exonic stabilities

# extract LFC values from exonic stability objects
cellline_lfcs <- lapply(unique(metadata$cell_line),function(cellline){
  lfc <- data.frame(
    gene=rownames(ex.20[[cellline]]$res),
    lfc_20=ex.20[[cellline]]$res$log2FoldChange,
    lfc_62=ex.62[[cellline]]$res$log2FoldChange,
    lfc_60=ex.60[[cellline]]$res$log2FoldChange
  )
  lfc.1 <- lfc[apply(lfc,1,function(x)!any(is.na(x[-1]))),]
  lfc.1
})
names(cellline_lfcs) <- unique(metadata$cell_line)
cellline_lfcs.df <- bind_rows(cellline_lfcs,.id="cellline")

# filter for values present in all 16 cell lines and pivot wider (one row per gene)
cellline_lfcs.df.1 <- cellline_lfcs.df %>% group_by(gene) %>% mutate(n=n()) %>% ungroup %>% filter(n==16) %>% pivot_wider(id_cols = gene,names_from=cellline,values_from=c(lfc_20,lfc_62,lfc_60))

# check the dataframe dimensions
dim(cellline_lfcs.gn.df.1)

# extract the comparison ("20", "60", "62") and cell line names
cellline_lfcs.df.1.metadata <- data.frame(
  comparison=case_when(grepl("_20_",colnames(cellline_lfcs.df.1)[-1]) ~ "20",grepl("_60_",colnames(cellline_lfcs.df.1)[-1]) ~ "60",grepl("_62_",colnames(cellline_lfcs.df.1)[-1]) ~ "62",TRUE ~ "none"),
  cellline=gsub("^lfc_[026]+_(.*)$","\\1",colnames(cellline_lfcs.df.1)[-1])
)
rownames(cellline_lfcs.df.1.metadata) <- colnames(cellline_lfcs.df.1)[-1]

cellline_lfcs.df.1.celllines <- cellline_lfcs.df.1.metadata[ ,"cellline",drop=FALSE]
cellline_lfcs.df.1b <- cellline_lfcs.df.1
colnames(cellline_lfcs.df.1b) <- gsub("lfc_([026]+)_(.*)$","\\2_\\1",colnames(cellline_lfcs.df.1))
rownames(cellline_lfcs.df.1.celllines) <- colnames(cellline_lfcs.df.1b)[-1]

# prepare plot colors
paletteLength <- 250
myColorFunc <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))
myBreaks <- c(
  seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
  seq(1/paletteLength, 1, length.out=floor(paletteLength/2))
)
rgb2hex <- function(r, g, b) {rgb(r, g, b, maxColorValue = 255)}
# load ENCODE portal color scheme
enc_colors <- read.delim("ENCODE4.16CellLines.Colors.FromPortal.txt")
tmp.colors <- data.frame(color=apply(enc_colors[enc_colors$style=="normal",c("r","g","b")],1,function(x){ rgb2hex(x["r"],x["g"],x["b"]) }))
rownames(tmp.colors) <- enc_colors$cell_line[enc_colors$style=="normal"]
tmp.colors.2 <- tmp.colors[,"color"]
names(tmp.colors.2) <- rownames(tmp.colors)
cellline_lfcs.df.1.colors <- list(cellline = tmp.colors.2)

## 2h v 0h heatmap

# generate correlation matrix
mat.20 <- as.matrix(cor(cellline_lfcs.df.1b[,grepl("_20$",colnames(cellline_lfcs.df.1b))]))
rownames(mat.20) <- gsub("_[026]+$","",rownames(mat.20))
colnames(mat.20) <- gsub("_[026]+$","",rownames(mat.20))

# get annotation
annot.20 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_20",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.20) <- annot.20$cell_line

# select colors
colors.20 <- list(cell_line=tmp.colors.2)

# plot without legend (default clustering)
hm1 <- pheatmap(mat = mat.20,
                scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
                # annotation_col=cellline_lfcs.df.1.celllines,annotation_row=cellline_lfcs.df.1.celllines,
                # annotation_colors = cellline_lfcs.df.1.colors,
                annotation_col=annot.20,
                annotation_row=annot.20,
                annotation_colors = colors.20,
                legend=FALSE,annotation_legend=FALSE)
hm1

## 6h v 2h heatmap

# generate correlation matrix
mat.62 <- as.matrix(cor(cellline_lfcs.df.1b[,grepl("_62$",colnames(cellline_lfcs.df.1b))]))
rownames(mat.62) <- gsub("_[026]+$","",rownames(mat.62))
colnames(mat.62) <- gsub("_[026]+$","",rownames(mat.62))

# get annotation
annot.62 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_62",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.62) <- annot.62$cell_line

# select colors
colors.62 <- list(cell_line=tmp.colors.2)

# plot without legend (default clustering)
hm2 <- pheatmap(mat = mat.62,
                scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
                # annotation_col=cellline_lfcs.df.1.celllines,annotation_row=cellline_lfcs.df.1.celllines,
                # annotation_colors = cellline_lfcs.df.1.colors,
                annotation_col=annot.62,
                annotation_row=annot.62,
                annotation_colors = colors.62,
                legend=FALSE,annotation_legend=FALSE)
hm2

# save heatmaps together as a png and a pdf
ggsave("ENCODE4.16celllines.stability_correlations.exonic.2v0.6v2.202211114.png",plot_grid(plotlist = list(hm1$gtable,hm2$gtable),nrow=1),width=10,height=5)
ggsave("ENCODE4.16celllines.stability_correlations.exonic.2v0.6v2.202211114.pdf",plot_grid(plotlist = list(hm1$gtable,hm2$gtable),nrow=1),width=10,height=5)

### process genic stabilities
               
cellline_lfcs.gn <- lapply(unique(metadata$cell_line),function(cellline){
  lfc <- data.frame(
    gene=rownames(gn.20[[cellline]]$res),
    lfc_20=gn.20[[cellline]]$res$log2FoldChange,
    lfc_62=gn.62[[cellline]]$res$log2FoldChange,
    lfc_60=gn.60[[cellline]]$res$log2FoldChange
  )
  lfc.1 <- lfc[apply(lfc,1,function(x)!any(is.na(x[-1]))),]
  lfc.1
})
names(cellline_lfcs.gn) <- unique(metadata$cell_line)
cellline_lfcs.gn.df <- bind_rows(cellline_lfcs.gn,.id="cellline")

cellline_lfcs.gn.df.1 <- cellline_lfcs.gn.df %>% group_by(gene) %>% mutate(n=n()) %>% ungroup %>% filter(n==16) %>% pivot_wider(id_cols = gene,names_from=cellline,values_from=c(lfc_20,lfc_62,lfc_60))
cellline_lfcs.gn.df.1b <- cellline_lfcs.gn.df.1
colnames(cellline_lfcs.gn.df.1b) <- gsub("lfc_([026]+)_(.*)$","\\2_\\1",colnames(cellline_lfcs.gn.df.1))

cellline_lfcs.df.1b[cellline_lfcs.df.1b$gene=="ENSG00000227232.5/WASH7P",] %>% as.data.frame
cellline_lfcs.gn.df.1b[cellline_lfcs.gn.df.1b$gene=="ENSG00000227232.5/WASH7P",] %>% as.data.frame

mat.gn.20 <- as.matrix(cor(cellline_lfcs.gn.df.1b[,grepl("_20$",colnames(cellline_lfcs.gn.df.1b))]))
rownames(mat.gn.20) <- gsub("_[026]+$","",rownames(mat.gn.20))
colnames(mat.gn.20) <- gsub("_[026]+$","",rownames(mat.gn.20))

annot.gn.20 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_20",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.gn.20) <- annot.gn.20$cell_line
colors.gn.20 <- list(cell_line=tmp.colors.2)

# plot without legend
hm1.gn <- pheatmap(mat = mat.gn.20,
                scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
                # annotation_col=cellline_lfcs.df.1.celllines,annotation_row=cellline_lfcs.df.1.celllines,
                # annotation_colors = cellline_lfcs.df.1.colors,
                annotation_col=annot.gn.20,
                annotation_row=annot.gn.20,
                annotation_colors = colors.gn.20,
                legend=FALSE,annotation_legend=FALSE)
hm1.gn


mat.gn.62 <- as.matrix(cor(cellline_lfcs.gn.df.1b[,grepl("_62$",colnames(cellline_lfcs.gn.df.1b))]))
rownames(mat.gn.62) <- gsub("_[026]+$","",rownames(mat.gn.62))
colnames(mat.gn.62) <- gsub("_[026]+$","",rownames(mat.gn.62))

annot.gn.62 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_62",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.gn.62) <- annot.gn.62$cell_line
colors.gn.62 <- list(cell_line=tmp.colors.2)

# plot without legend
hm2.gn <- pheatmap(mat = mat.gn.62,
                scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
                # annotation_col=cellline_lfcs.df.1.celllines,annotation_row=cellline_lfcs.df.1.celllines,
                # annotation_colors = cellline_lfcs.df.1.colors,
                annotation_col=annot.gn.62,
                annotation_row=annot.gn.62,
                annotation_colors = colors.gn.62,
                legend=FALSE,annotation_legend=FALSE)
hm2.gn

ggsave("ENCODE4.16celllines.stability_correlations.genic.2v0.6v2.202211114.png",plot_grid(plotlist = list(hm1.gn$gtable,hm2.gn$gtable),nrow=1),width=10,height=5)
ggsave("ENCODE4.16celllines.stability_correlations.genic.2v0.6v2.202211114.pdf",plot_grid(plotlist = list(hm1.gn$gtable,hm2.gn$gtable),nrow=1),width=10,height=5)

### build legend for heatmaps

# plot with legend
hm.forlegends <- 
  pheatmap(mat = mat.20,
           scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
           annotation_col=annot.20,
           annotation_row=annot.20,
           annotation_colors = colors.20
  )

# extract legends (legend and annotation legend)
leg1 <- gtable::gtable_filter(x = hm.forlegends$gtable, pattern="^legend$")
leg2 <- gtable::gtable_filter(x = hm.forlegends$gtable, pattern="^annotation_legend$")

# save legends as png and pdf
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend1.20221114.png",leg1,width=1,height=5)
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend2.20221114.png",leg2,width=2,height=10)
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend1.20221114.pdf",leg1,width=1,height=5)
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend2.20221114.pdf",leg2,width=2,height=10)

### end

