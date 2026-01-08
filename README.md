# ENCODE-RNA-dynamics-paper

This repo contains data generation and analysis code used in this paper:

Isoform and pathway-specific regulation of post-transcriptional RNA processing in human cells.
Karan Bedi, Brian Magnuson, Ishwarya Venkata Narayanan, Ariel McShane, Mario Ashaka, Michelle T. Paulsen, Thomas E. Wilson, Mats Ljungman
bioRxiv 2024.06.12.598705; doi: https://doi.org/10.1101/2024.06.12.598705 

----------------------------------------------------------------------------------------

1. compute-stability.R

Calculates the relative stability of exons during the 0 to 2 hour and 2 to 6 hour chase periods. RNA stability is represented as Log2 fold change values, which are further scaled from 0 to 1 for each cell line. 

Input files: general_input_files.zip
Supplemental table S13 lists the ENCODE accession numbers for the genic features quantifications files that contain the exon_sense read counts used in this analysis.

2.  

Determines the relative stability of intronic RNA, enhancer RNA (eRNA), PROMPTs and readthrough (RT) RNA, using exonic RNA as the reference. This analysis is related to Main Figure 1b. 

3. gene_exonic_stability_correlation_heatmaps.R


