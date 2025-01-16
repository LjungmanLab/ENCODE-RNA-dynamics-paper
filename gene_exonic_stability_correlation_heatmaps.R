library(tidyverse)
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(cowplot)

#install.packages("cowplot")

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))


ex.20 <- readRDS("ENCODE4.16CellLines.Stability.2v0.ExonSense.rds")
ex.60 <- readRDS("ENCODE4.16CellLines.Stability.6v0.ExonSense.rds")
ex.62 <- readRDS("ENCODE4.16CellLines.Stability.6v2.ExonSense.rds")

gn.20 <- readRDS("ENCODE4.16CellLines.Stability.2v0.GeneSense.rds")
gn.60 <- readRDS("ENCODE4.16CellLines.Stability.6v0.GeneSense.rds")
gn.62 <- readRDS("ENCODE4.16CellLines.Stability.6v2.GeneSense.rds")


names(ex.60)
names(ex.60$A673)
head(ex.60$A673$res)

# es.stab.60 <- read.delim("ENCODE4.16CellLines.Stability.6v0.ExonSense.Annotated.txt",header=TRUE,sep="\t")
# es.stab.62 <- read.delim("ENCODE4.16CellLines.Stability.6v2.ExonSense.Annotated.txt",header=TRUE,sep="\t")
# es.stab.20 <- read.delim("ENCODE4.16CellLines.Stability.2v0.ExonSense.Annotated.txt",header=TRUE,sep="\t")
# gs.stab.60 <- read.delim("ENCODE4.16CellLines.Stability.6v0.GeneSense.Annotated.txt",header=TRUE,sep="\t")
# gs.stab.62 <- read.delim("ENCODE4.16CellLines.Stability.6v2.GeneSense.Annotated.txt",header=TRUE,sep="\t")
# gs.stab.20 <- read.delim("ENCODE4.16CellLines.Stability.2v0.GeneSense.Annotated.txt",header=TRUE,sep="\t")

metadata <- read.delim("ENCODE4.16CellLines.BruChase.Sample.Metadata.txt",sep="\t",header=TRUE)
head(metadata)
annot <- read.delim("hg38.gencode_29.ensembl_gene.transcript_biotype_name_entrez.gencode.txt",sep="\t",header=TRUE)
head(annot)

head(es.stab.60)
colnames(ex.20$A673$res)
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
head(cellline_lfcs.df)

cellline_lfcs.df.1 <- cellline_lfcs.df %>% group_by(gene) %>% mutate(n=n()) %>% ungroup %>% filter(n==16) %>% pivot_wider(id_cols = gene,names_from=cellline,values_from=c(lfc_20,lfc_62,lfc_60))

dim(cellline_lfcs.gn.df.1)
head(cellline_lfcs.df.1)


cellline_lfcs.df.1.metadata <- data.frame(
  comparison=case_when(grepl("_20_",colnames(cellline_lfcs.df.1)[-1]) ~ "20",grepl("_60_",colnames(cellline_lfcs.df.1)[-1]) ~ "60",grepl("_62_",colnames(cellline_lfcs.df.1)[-1]) ~ "62",TRUE ~ "none"),
  cellline=gsub("^lfc_[026]+_(.*)$","\\1",colnames(cellline_lfcs.df.1)[-1])
)
rownames(cellline_lfcs.df.1.metadata) <- colnames(cellline_lfcs.df.1)[-1]
cellline_lfcs.df.1.metadata

cellline_lfcs.df.1.celllines <- cellline_lfcs.df.1.metadata[ ,"cellline",drop=FALSE]
cellline_lfcs.df.1.celllines

cellline_lfcs.df.1b <- cellline_lfcs.df.1

colnames(cellline_lfcs.df.1b) <- gsub("lfc_([026]+)_(.*)$","\\2_\\1",colnames(cellline_lfcs.df.1))
colnames(cellline_lfcs.df.1b)         

rownames(cellline_lfcs.df.1.celllines) <- colnames(cellline_lfcs.df.1b)[-1]

paletteLength <- 250
myColorFunc <- colorRampPalette(brewer.pal(n = 7, name ="RdBu"))
myBreaks <- c(
  seq(-1, 0, length.out=ceiling(paletteLength/2) + 1), 
  seq(1/paletteLength, 1, length.out=floor(paletteLength/2))
)


rgb2hex <- function(r, g, b) {rgb(r, g, b, maxColorValue = 255)}

enc_colors <- read.delim("ENCODE4.16CellLines.Colors.FromPortal.txt")
enc_colors

tail(cellline_lfcs.df.1.celllines)

tmp.colors <- data.frame(color=apply(enc_colors[enc_colors$style=="normal",c("r","g","b")],1,function(x){ rgb2hex(x["r"],x["g"],x["b"]) }))
rownames(tmp.colors) <- enc_colors$cell_line[enc_colors$style=="normal"]
tmp.colors
tmp.colors.2 <- tmp.colors[,"color"]
names(tmp.colors.2) <- rownames(tmp.colors)

cellline_lfcs.df.1.colors <- list(cellline = tmp.colors.2)
head(cellline_lfcs.df.1.colors)


cellline_lfcs.df.1b[,grepl("_20$",colnames(cellline_lfcs.df.1b))] %>% head

mat.20 <- as.matrix(cor(cellline_lfcs.df.1b[,grepl("_20$",colnames(cellline_lfcs.df.1b))]))
rownames(mat.20) <- gsub("_[026]+$","",rownames(mat.20))
colnames(mat.20) <- gsub("_[026]+$","",rownames(mat.20))
head(mat.20)
annot.20 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_20",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.20) <- annot.20$cell_line
annot.20
colors.20 <- list(cell_line=tmp.colors.2)
colors.20

hm1 <- pheatmap(mat = mat.20,
                scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
                # annotation_col=cellline_lfcs.df.1.celllines,annotation_row=cellline_lfcs.df.1.celllines,
                # annotation_colors = cellline_lfcs.df.1.colors,
                annotation_col=annot.20,
                annotation_row=annot.20,
                annotation_colors = colors.20,
                legend=FALSE,annotation_legend=FALSE)
hm1

mat.62 <- as.matrix(cor(cellline_lfcs.df.1b[,grepl("_62$",colnames(cellline_lfcs.df.1b))]))
rownames(mat.62) <- gsub("_[026]+$","",rownames(mat.62))
colnames(mat.62) <- gsub("_[026]+$","",rownames(mat.62))
head(mat.62)
annot.62 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_62",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.62) <- annot.62$cell_line
colors.62 <- list(cell_line=tmp.colors.2)

hm2 <- pheatmap(mat = mat.62,
                scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
                # annotation_col=cellline_lfcs.df.1.celllines,annotation_row=cellline_lfcs.df.1.celllines,
                # annotation_colors = cellline_lfcs.df.1.colors,
                annotation_col=annot.62,
                annotation_row=annot.62,
                annotation_colors = colors.62,
                legend=FALSE,annotation_legend=FALSE)
hm2

ggsave("ENCODE4.16celllines.stability_correlations.exonic.2v0.6v2.202211114.png",plot_grid(plotlist = list(hm1$gtable,hm2$gtable),nrow=1),width=10,height=5)
ggsave("ENCODE4.16celllines.stability_correlations.exonic.2v0.6v2.202211114.pdf",plot_grid(plotlist = list(hm1$gtable,hm2$gtable),nrow=1),width=10,height=5)

###

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
head(cellline_lfcs.gn.df)
head(cellline_lfcs.df)

cellline_lfcs.gn.df.1 <- cellline_lfcs.gn.df %>% group_by(gene) %>% mutate(n=n()) %>% ungroup %>% filter(n==16) %>% pivot_wider(id_cols = gene,names_from=cellline,values_from=c(lfc_20,lfc_62,lfc_60))
cellline_lfcs.gn.df.1b <- cellline_lfcs.gn.df.1
colnames(cellline_lfcs.gn.df.1b) <- gsub("lfc_([026]+)_(.*)$","\\2_\\1",colnames(cellline_lfcs.gn.df.1))
head(cellline_lfcs.gn.df.1b %>% as.data.frame)
head(cellline_lfcs.df.1b %>% as.data.frame)

cellline_lfcs.df.1b[cellline_lfcs.df.1b$gene=="ENSG00000227232.5/WASH7P",] %>% as.data.frame
cellline_lfcs.gn.df.1b[cellline_lfcs.gn.df.1b$gene=="ENSG00000227232.5/WASH7P",] %>% as.data.frame

mat.gn.20 <- as.matrix(cor(cellline_lfcs.gn.df.1b[,grepl("_20$",colnames(cellline_lfcs.gn.df.1b))]))
rownames(mat.gn.20) <- gsub("_[026]+$","",rownames(mat.gn.20))
colnames(mat.gn.20) <- gsub("_[026]+$","",rownames(mat.gn.20))
head(mat.gn.20)
head(mat.20)

annot.gn.20 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_20",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.gn.20) <- annot.gn.20$cell_line
annot.gn.20
colors.gn.20 <- list(cell_line=tmp.colors.2)
colors.gn.20

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
head(mat.gn.62)
head(mat.62)

annot.gn.62 <- data.frame(cell_line=cellline_lfcs.df.1.celllines$cellline[grepl("_62",rownames(cellline_lfcs.df.1.celllines))])
rownames(annot.gn.62) <- annot.gn.62$cell_line
colors.gn.62 <- list(cell_line=tmp.colors.2)

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


###

hm.forlegends <- 
  pheatmap(mat = mat.20,
           scale = "none",cluster_rows = TRUE,cluster_cols=TRUE,color=myColorFunc(paletteLength),breaks=myBreaks,
           annotation_col=annot.20,
           annotation_row=annot.20,
           annotation_colors = colors.20
  )
#class(hm.forlegends)
#hm.forlegends$gtable

leg1 <- gtable::gtable_filter(x = hm.forlegends$gtable, pattern="^legend$")
leg2 <- gtable::gtable_filter(x = hm.forlegends$gtable, pattern="^annotation_legend$")

ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend1.20221114.png",leg1,width=1,height=5)
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend2.20221114.png",leg2,width=2,height=10)
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend1.20221114.pdf",leg1,width=1,height=5)
ggsave(filename = "ENCODE4.16celllines.stability_correlations.legend2.20221114.pdf",leg2,width=2,height=10)
# 
# names(cellline_lfcs)
# hms <- lapply(cellline_lfcs,function(x){
#   p <- pheatmap(mat=as.matrix(cor(x[,-1])),color=red(n=250),breaks = seq(-1,1,length.out = 251), cluster_rows = FALSE,cluster_cols = FALSE)
# })
# 
# plot_grid(plotlist = lapply(hms,function(x)x$gtable),nrow=4)
# 
# 



