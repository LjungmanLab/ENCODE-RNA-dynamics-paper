# I first calculated SI of all introns in all known transcripts in gc29.
# Then, based on a strict split (SS) reads filter of 3, I shortlisted all transcripts
# (for genes with more than 1 isoforms) in which ALL introns met the SS filter in the bruseq 0h sample.
# I then shortlisted only these introns from the chase 2h and 6h counterpart samples (without any SS filter)
# I was able to obtain 1076 genes, with 2549 transcript isoforms.
# Across the 16 cell lines, I obtained 20796 common introns.
# The common introns SI files are located in C:/bedik/ENCODE_analysis/PEsplicing

# Now, I want to take these 2549 transcript isoforms and look at their IsoPct changs across the 3 timepoints,
# using RSEM results from re_alignment_exonsONLY attenpt

# For each cell line, I find out what percentage of genes have isoforms that follow certain patterns,
# such as being dominant in all three time points or loose dominance over time, etc etc

library(dplyr)
library(ggplot2)
library(ggside)
library(tidyr)
library(tidyverse)
library(broom)
library(quantmod)
library(reshape2)
library(WebGestaltR)
library(stringr)
library(matrixStats)

write.old.files = TRUE

mytp = c("026h")

rsem_subfolder = mytp
# wd = paste0("C:/bedik/ENCODE_analysis/rsem/")

wd = getwd()

# si_dir = "C:/bedik/ENCODE_analysis/PEsplicing/"
si_dir = wd

# outdir = paste0(wd,mytp,"/SIbased_EO/") # EO is for exonONLY
outdir = wd

my_metadata1 = read.delim("./cell_lines.txt",sep = "\t",stringsAsFactors = F, header = F)
head(my_metadata1)
# i=1

my_matrix = matrix(ncol = 5, nrow = 16)

for (i in 1:nrow(my_metadata1)) {

  cell_line = my_metadata1[i, "V1"]

  new_cell_line = ifelse(cell_line=="Panc-1","panc1",cell_line)

  bruseq_SI = read.delim(paste0(si_dir,"/",new_cell_line,"_ss3_MR_ALL_TFgt1_ISI.bed"),
                         sep = "\t",header = F, stringsAsFactors = F)
  bruchase2h_SI = read.delim(paste0(si_dir,"/",new_cell_line,"_ss3_MR_ALL_TFgt1_ISI.bed"),
                             sep = "\t",header = F, stringsAsFactors = F)
  bruchase6h_SI = read.delim(paste0(si_dir,"/",new_cell_line,"_ss3_MR_ALL_TFgt1_ISI.bed"),
                             sep = "\t",header = F, stringsAsFactors = F)

  si_merge = Reduce(function(x, y) merge(x, y, by = colnames(x)[1:16]), list(bruseq_SI,bruchase2h_SI, bruchase6h_SI))
  head(si_merge)
  si_merge = si_merge[order(si_merge$V1,si_merge$V13,si_merge$V4, si_merge$V2),]

  si_tx_list = unique(sort(si_merge$V4))


  bruseq_IsoPct = read.delim(paste0(wd,"/",cell_line,"_0h_exonsONLY_meanIsoPct.txt"),
                             sep = "\t",header = F, stringsAsFactors = F)
  bruchase2h_IsoPct = read.delim(paste0(wd,"/",cell_line,"_2h_exonsONLY_meanIsoPct.txt"),
                                 sep = "\t",header = F, stringsAsFactors = F)
  bruchase6h_IsoPct = read.delim(paste0(wd,"/",cell_line,"_6h_exonsONLY_meanIsoPct.txt"),
                                 sep = "\t",header = F, stringsAsFactors = F)

  isopct_merge = Reduce(function(x, y) merge(x, y, by = colnames(x)[1]), list(bruseq_IsoPct,bruchase2h_IsoPct, bruchase6h_IsoPct))
  head(isopct_merge)
  str(isopct_merge)

  # Shortlist isoforms from SI tx list
  isopct_merge_sub1 = isopct_merge[isopct_merge$V1 %in% si_tx_list,]

  # Combine to get gene id's
  head(si_merge)
  tail(si_merge)
  ensg_enst = unique(si_merge[,c(13,4,14)])
  ensg_enst_2 = ensg_enst %>%
    separate(V14,c("col1","col2","col3"), sep = "@@") %>%
    separate(col2,c("GL","ENSG","symbol","IL"), sep = ":") %>%
    select(symbol,V13,V4) %>%
    unique() %>%
    as.data.frame()

  isopct_merge_sub2 = merge(ensg_enst_2,isopct_merge_sub1,by.x = "V4", by.y="V1")
  isopct_merge_sub2 = isopct_merge_sub2[,c(2,3,1,4:ncol(isopct_merge_sub2))]

  colnames(isopct_merge_sub2) = c("symbol","gene", "transcript", "IsoPct_0h", "IsoPct_2h", "IsoPct_6h")
  isopct_merge_sub2 = isopct_merge_sub2[order(isopct_merge_sub2$gene,isopct_merge_sub2$transcript),]
  str(isopct_merge_sub2)

  cols.num = c("IsoPct_0h", "IsoPct_2h", "IsoPct_6h")
  isopct_merge_sub2[cols.num] = sapply(isopct_merge_sub2[cols.num], as.numeric)

  # Add row wise variance column
  isopct_merge_sub2 = isopct_merge_sub2 %>%
    mutate(row_wise_var = rowVars(as.matrix(isopct_merge_sub2[,grep("IsoPct_",colnames(isopct_merge_sub2),value = T)]))) %>%
    as.data.frame()


  if(write.old.files){

    write.table(isopct_merge_sub2,paste0(outdir,cell_line,"_SIbased_exonsONLY_IsoPcts.txt"), sep = "\t", col.names = T,
                row.names = F, quote = F)
  }

  isopct_merge_sub3 = isopct_merge_sub2[order(isopct_merge_sub2$IsoPct_6h,-isopct_merge_sub2$IsoPct_0h),]


  isopct_merge_sub4 = isopct_merge_sub3 %>%
    select(., transcript, IsoPct_0h,IsoPct_2h, IsoPct_6h)


  isopct_merge_sub5 = transform(isopct_merge_sub2, SD_06=apply(isopct_merge_sub2[,c(4,6)],1, sd, na.rm = TRUE))
  isopct_merge_sub5$mean_IsoPct = rowMeans(isopct_merge_sub5[,c(4:6)])

  isopct_merge_sub5$trend = ifelse(isopct_merge_sub5$IsoPct_0h>isopct_merge_sub5$IsoPct_2h &
                                     isopct_merge_sub5$IsoPct_2h>isopct_merge_sub5$IsoPct_6h,"down-down",
                                   ifelse(isopct_merge_sub5$IsoPct_0h<isopct_merge_sub5$IsoPct_2h &
                                            isopct_merge_sub5$IsoPct_2h<isopct_merge_sub5$IsoPct_6h,"up-up",
                                          ifelse(isopct_merge_sub5$IsoPct_0h>isopct_merge_sub5$IsoPct_2h &
                                                   isopct_merge_sub5$IsoPct_2h<isopct_merge_sub5$IsoPct_6h,"down-up",
                                                 ifelse(isopct_merge_sub5$IsoPct_0h<isopct_merge_sub5$IsoPct_2h &
                                                          isopct_merge_sub5$IsoPct_2h>isopct_merge_sub5$IsoPct_6h,"up-down","nochange"))))

  allgenes = length(unique(sort(isopct_merge_sub5$symbol)))
  # 2899 A673

  # # Write background gene list
  # write.table(as.data.frame(unique(sort(isopct_merge_sub5$symbol))),paste0(outdir,"/",cell_line,"_exonsONLY_background.txt"),
  #             col.names = F, row.names = F, quote = F)


  pct0hgt50_1_all = isopct_merge_sub5[isopct_merge_sub5$IsoPct_0h>50,]
  length(unique(sort(pct0hgt50_1_all$symbol))) # 1433 A673
  # d026
  pct0hgt50_2_pct2hgt50_pct6hgt50 = pct0hgt50_1_all[pct0hgt50_1_all$IsoPct_2h>50 & pct0hgt50_1_all$IsoPct_6h>50,]
  d026 = length(unique(sort(pct0hgt50_2_pct2hgt50_pct6hgt50$symbol))) # 922/1433 A673
  # d02r6
  pct0hgt50_3_pct2hgt50_pct6hlte50 = pct0hgt50_1_all[pct0hgt50_1_all$IsoPct_2h>50 & pct0hgt50_1_all$IsoPct_6h<=50,]
  d02r6 = length(unique(sort(pct0hgt50_3_pct2hgt50_pct6hlte50$symbol))) # 98/1433 A673
  # d0r26
  pct0hgt50_4_pct2hlte50_pct6hlte50 = pct0hgt50_1_all[pct0hgt50_1_all$IsoPct_2h<=50 & pct0hgt50_1_all$IsoPct_6h<=50,]
  d0r26 = length(unique(sort(pct0hgt50_4_pct2hlte50_pct6hlte50$symbol))) # 348/1433 A673
  # d06r2
  pct0hgt50_5_pct2hlte50_pct6hgt50 = pct0hgt50_1_all[pct0hgt50_1_all$IsoPct_2h<=50 & pct0hgt50_1_all$IsoPct_6h>50,]
  d06r2 = length(unique(sort(pct0hgt50_5_pct2hlte50_pct6hgt50$symbol))) # 65/1433 A673


  pct0hlte50_1_all = isopct_merge_sub5[isopct_merge_sub5$IsoPct_0h<=50,]
  length(unique(sort(pct0hlte50_1_all$symbol))) # 2899 A673
  # r026
  pct0hlte50_2_pct2hlte50_pct6hlte50 = pct0hlte50_1_all[pct0hlte50_1_all$IsoPct_2h<=50 & pct0hlte50_1_all$IsoPct_6h<=50,]
  r026 = length(unique(sort(pct0hlte50_2_pct2hlte50_pct6hlte50$symbol))) # 2760/2899
  # r02d6
  pct0hlte50_3_pct2hlte50_pct6hgt50 = pct0hlte50_1_all[pct0hlte50_1_all$IsoPct_2h<=50 & pct0hlte50_1_all$IsoPct_6h>50,]
  r02d6 = length(unique(sort(pct0hlte50_3_pct2hlte50_pct6hgt50$symbol))) # 228/2899
  # r0d26
  pct0hlte50_4_pct2hgt50_pct6hgt50 = pct0hlte50_1_all[pct0hlte50_1_all$IsoPct_2h>50 & pct0hlte50_1_all$IsoPct_6h>50,]
  r0d26 = length(unique(sort(pct0hlte50_4_pct2hgt50_pct6hgt50$symbol))) # 396/2899
  # r06d2
  pct0hlte50_5_pct2hgt50_pct6hlte50 = pct0hlte50_1_all[pct0hlte50_1_all$IsoPct_2h>50 & pct0hlte50_1_all$IsoPct_6h<=50,]
  r06d2 = length(unique(sort(pct0hlte50_5_pct2hgt50_pct6hlte50$symbol))) # 83/2899


  # d026 d02r6 and d0r26 OF MOST INTEREST
  # r0d26 and r02d6 OF MOST INTEREST

  # % always recessive
  pAR = round(100*r026/allgenes, digits = 0)

  # % always dominant
  pAD = round(100*d026/allgenes, digits = 0)
  # # Write always dominant gene list
  # write.table(as.data.frame(unique(sort(pct0hgt50_2_pct2hgt50_pct6hgt50$symbol))),paste0(outdir,"/",cell_line,"_exonsONLY_AD.txt"),
  #             col.names = F, row.names = F, quote = F)

  # % that loose dominance over time
  pLD = round(100*(d02r6+d0r26)/allgenes, digits = 0)
  # Write loose dominance over time gene list
  pLDgl = c(unique(sort(pct0hgt50_3_pct2hgt50_pct6hlte50$symbol)),unique(sort(pct0hgt50_4_pct2hlte50_pct6hlte50$symbol)))
  # write.table(as.data.frame(unique(sort(pLDgl))),paste0(outdir,"/",cell_line,"_exonsONLY_LD.txt"),
  #             col.names = F, row.names = F, quote = F)

  # % that gain dominance over time
  pGD = round(100*(r0d26+r02d6)/allgenes, digits = 0)
  # Write gain dominance over time gene list
  pGDgl = c(unique(sort(pct0hlte50_4_pct2hgt50_pct6hgt50$symbol)),unique(sort(pct0hlte50_3_pct2hlte50_pct6hgt50$symbol)))
  # write.table(as.data.frame(unique(sort(pGDgl))),paste0(outdir,"/",cell_line,"_exonsONLY_GD.txt"),
  #             col.names = F, row.names = F, quote = F)

  my_matrix[i,] = c(cell_line,pAR,pAD,pLD,pGD)


}

pct_df = as.data.frame(my_matrix)
str(pct_df)
colnames(pct_df) = c("cell_line", "AR", "AD", "LD","GD")
cols.num = colnames(pct_df)[2:ncol(pct_df)]
pct_df[cols.num] = sapply(pct_df[cols.num], as.numeric)

plot.data.filename = paste0(outdir,"/A0_gene_isoform_percentages")

if(write.old.files){
  write.table(pct_df, paste0(plot.data.filename,".txt"), sep = "\t", col.names = T, row.names = F,
              quote = F)
}

pct_dfm = melt(pct_df[,c(1,3:ncol(pct_df))],id.vars = "cell_line")


myplot = ggplot(pct_dfm, aes(x = cell_line, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  labs(x = "Cell Lines", y = "% genes") +
  guides(fill=guide_legend(title="Gene isoform type")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 9, angle = 45))
myplot

if(write.old.files){
  ggsave(myplot, filename = paste0(plot.data.filename,".jpg"),dpi = 600)
  ggsave(myplot, filename = paste0(plot.data.filename,".pdf"),dpi = 600)
}



