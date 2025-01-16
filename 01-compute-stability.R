options(stringsAsFactors = FALSE)

library(tidyverse)
library(cowplot)
library(hexbin)
library(ggdendro)
library(ggrepel)
library(DESeq2)
library(pheatmap)
library(openxlsx)
library(scales)

getwd()

wd = getwd()

outdir = paste0(wd,"/outputs")
dir.create(outdir,showWarnings = F, recursive = T)

tmp_dir = paste0(wd,"/tmp")
dir.create(tmp_dir,showWarnings = F,recursive = T)

#################################################
# Part A: Shortlist genes based on expression

# Biotype info
biotype.df = read.delim("input_files/hg38.gencode_29.ensembl_gene.transcript_biotype_name_entrez.gencode.txt",
                        sep = "\t", header = T, na.strings = c(""," ","NA"),check.names = F)
biotype.df.2 = biotype.df %>%
  mutate(ensembl_gene = paste(ensembl_gene,gene_name,sep = "/")) %>%
  select(ensembl_gene,gene_type) %>%
  as.data.frame()

# Find out all genes at 0h and 2h whose mean exon sense RPKM (mean across 2 replicates) is >= 0.5 across all cell lines
es.df = read.delim("input_files/ENCODE4.16CellLines.ReplicatesMeanExonSenseRPKM.bed", sep = "\t", header = T, na.strings = c(""," ","NA"),
                   check.names = F)
head(es.df)
colnames(es.df)
fixed.gene.cols = colnames(es.df)[1:7]

es.df.2 = es.df %>%
  select(all_of(fixed.gene.cols),grep(paste(c("_ES_0h","_ES_2h"),collapse = "|"),colnames(.),value = T)) %>%
  filter(if_all(ends_with("_0h"), ~ . > 0.5)) %>%
  filter(if_all(ends_with("_2h"), ~ . > 0.5)) %>%
  as.data.frame()
es.df.2 = left_join(es.df.2,biotype.df.2, by = c("gene"="ensembl_gene")) %>%
  relocate(gene_type, .after = strand) %>%
  select(-featuretype) %>%
  filter(gene_type == "protein_coding") %>%
  as.data.frame()

# Remove the following gene because of duplicated symbol
# ENSG00000187522.15/HSPA14
# ENSG00000284024.2/HSPA14
del.genes = c("ENSG00000187522.15/HSPA14","ENSG00000284024.2/HSPA14")
es.df.2 = es.df.2[! (es.df.2$gene %in% del.genes),]
colnames(es.df.2) = gsub("Panc1","panc1",colnames(es.df.2))

#################################################
# Part B: Compute stabilities

metadata <- read.delim("input_files/ENCODE4.16CellLines.BruChase.Sample.Metadata.txt",sep="\t",header=TRUE)
head(metadata)

counts <- read.delim("input_files/ENCODE4.16CellLines.BruChase.FeatureTypes.Counts.bed",sep="\t",header = T)
head(counts)
unique(counts$featuretype)
# filter for exon_sense and remove featuretype column
counts <- counts[counts$featuretype=="exon_sense",-which(colnames(counts)=="featuretype")]
rownames(counts) <- counts$gene
identical(colnames(counts)[-c(1:6)],metadata$sample_id)
# TRUE

# massive set of chase comparisons, loop through all cell lines and compare 2 chase time points (2 v 0, and 6 v 2)

mysamples.s20 <- metadata$sample_id[metadata$chase_time %in% c(0,2)]
mysamples.s62 <- metadata$sample_id[metadata$chase_time %in% c(2,6)]

# 2h vs 0h Log2FC
cl.s20.ex.results <- lapply(unique(metadata$cell_line),function(cl){
  print(cl)
  mysamples <- metadata$sample_id[metadata$cell_line==cl & metadata$sample_id %in% mysamples.s20]
  coldata <- data.frame(
    chase_time = factor(metadata$chase_time[metadata$sample_id %in% mysamples])
  )
  rownames(coldata) <- mysamples
  cts <- round(counts[,mysamples])
  dds <- DESeqDataSetFromMatrix(cts,coldata,design = ~chase_time)
  dds <- DESeq(dds)
  rn <- resultsNames(dds)
  res <- results(dds,name=rn[2]) #
  shrunk <- lfcShrink(dds,coef = 2,type = "apeglm")
  return(list(dds=dds,res=res,shrunk=shrunk))
})
names(cl.s20.ex.results) <- unique(metadata$cell_line)
s20.LFC <- sapply(cl.s20.ex.results,function(cl){
  cl$res$log2FoldChange
})
s20.LFC.df <- as.data.frame(s20.LFC,row.names = rownames(counts))
colnames(s20.LFC.df) <- paste0(colnames(s20.LFC.df),"_2vs0_log2FC")


# 6h vs 2h Log2FC
cl.s62.ex.results <- lapply(unique(metadata$cell_line),function(cl){
  print(cl)
  mysamples <- metadata$sample_id[metadata$cell_line==cl & metadata$sample_id %in% mysamples.s62]
  coldata <- data.frame(
    chase_time = factor(metadata$chase_time[metadata$sample_id %in% mysamples])
  )
  rownames(coldata) <- mysamples
  cts <- round(counts[,mysamples])
  dds <- DESeqDataSetFromMatrix(cts,coldata,design = ~chase_time)
  dds <- DESeq(dds)
  rn <- resultsNames(dds)
  res <- results(dds,name=rn[2]) #
  shrunk <- lfcShrink(dds,coef = 2,type = "apeglm")
  return(list(dds=dds,res=res,shrunk=shrunk))
})
names(cl.s62.ex.results) <- unique(metadata$cell_line)
s62.LFC <- sapply(cl.s62.ex.results,function(cl){
  cl$res$log2FoldChange
})
s62.LFC.df <- as.data.frame(s62.LFC,row.names = rownames(counts))
colnames(s62.LFC.df) <- paste0(colnames(s62.LFC.df),"_6vs2_log2FC")

# Merge the 2 LFC df's into 1

stab.LFC.df = merge(s20.LFC.df,s62.LFC.df,by = 0, all = TRUE)
colnames(stab.LFC.df)[1] = "genes"
stab.sorted.cols = sort(colnames(stab.LFC.df)[2:ncol(stab.LFC.df)])
stab.LFC.df = stab.LFC.df[,c("genes",stab.sorted.cols)]

# Remove rows where all log2FC values are 0 or NA
stab.LFC.df.2 = stab.LFC.df[rowSums(is.na(stab.LFC.df[2:ncol(stab.LFC.df)])) != ncol(stab.LFC.df[2:ncol(stab.LFC.df)]), ]

# Shortlist only expressed genes from Part A
stab.LFC.df.3 = stab.LFC.df.2 %>%
  filter(genes %in% es.df.2$gene)
head(stab.LFC.df.3)
colnames(stab.LFC.df.3)[2:ncol(stab.LFC.df.3)] = gsub("_log2FC","_ES.LFC",colnames(stab.LFC.df.3)[2:ncol(stab.LFC.df.3)])
stab.fixed.cols = colnames(stab.LFC.df.3)[1]


# Rescale the values between 0 and 1
stab.df.all.pc.scaled = data.frame(matrix(nrow = nrow(stab.LFC.df.3), ncol = ncol(stab.LFC.df.3)))
str(stab.df.all.pc.scaled)

stab.df.all.pc.scaled[,c(1)] = stab.LFC.df.3[,c(1)]
colnames(stab.df.all.pc.scaled) = colnames(stab.LFC.df.3)
colnames(stab.df.all.pc.scaled)[2:ncol(stab.df.all.pc.scaled)] = paste0(colnames(stab.df.all.pc.scaled)[2:ncol(stab.df.all.pc.scaled)],".scaled")
for (i in 2:ncol(stab.df.all.pc.scaled)) {

  stab.df.all.pc.scaled[,i] = rescale(stab.LFC.df.3[,i], to = c(0,1))

}

# Merge bed file columns
stab.df.all.pc.scaled = left_join(stab.df.all.pc.scaled,es.df.2[,1:6],join_by("genes" == "gene")) %>%
  relocate(c(chr,start,end), .before = "genes") %>%
  relocate(c(score,strand),.after = "genes")

write.xlsx(stab.df.all.pc.scaled,paste0(outdir,"/scaled_stabilities.xlsx"),overwrite = T)












