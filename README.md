# ENCODE-RNA-dynamics-paper

This repo contains data generation and analysis code used in this paper:

Isoform and pathway-specific regulation of post-transcriptional RNA processing in human cells.
Karan Bedi, Brian Magnuson, Ishwarya Venkata Narayanan, Ariel McShane, Mario Ashaka, Michelle T. Paulsen, Thomas E. Wilson, Mats Ljungman
bioRxiv 2024.06.12.598705; doi: https://doi.org/10.1101/2024.06.12.598705 

----------------------------------------------------------------------------------------

1. **compute-stability.R**

    Calculates the relative stability of exonic RNA during the 0 to 2 hour and 2 to 6 hour chase periods. RNA stability is represented as Log2 fold change values, which are further scaled from 0 to 1 for each cell line. 
- Input files: general_input_files.zip
- Supplemental Table S13 provides the ENCODE accession numbers for the genic features quantifications files on the ENCODE portal (https://www.encodeproject.org/), containing the exon_sense read counts used in this analysis.

2. **<>Fig 1B**

    Determines the relative stability of intronic RNA, enhancer RNA (eRNA), PROMPTs and readthrough (RT) RNA, using exonic RNA as the reference. This analysis is related to Main Figure 1B.
-  Input file: Supplemental Table S2, sheet "".

3. **01D_bru026_SI_plots.R**

    Generates plots of splicing indices for common introns across 16 cell lines at 0, 2 and 6 hours, as displayed in Main Figure 1D.
- Input files: input_files_Fig01D.zip

4. **02AB_gene_exonic_stability_correlation_heatmaps.R**

    Computes the correlation coefficients among the 16 cell lines using exonic relative stability values (Log2 fold change). This analysis corresponds to Main Figures 2A and 2B.

5. **<> dist code**

<>

6. **04F_026h_SIbased_exonsONLY_IsoPct_trends.R**

    Measures isoform percentages and their prevalence over time. Findings are shown in Main Figure 4F. 
- Input files: input_files_Fig04F_1.zip, input_files_Fig04F_2.zip, input_files_Fig04F_3.zip, input_files_Fig04F_4.zip, input_files_Fig04F_5.zip

7. **<> Fig 6**

<>

8. **suppFig_12and13.R**

    Defines thresholds for classifying introns as "Stable" or "Retained", and generates the plots shown in Supplemental Figures 12 and 13.
- Input files: Supplemental Table S8, sheet "6vs0".

