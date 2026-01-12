# ENCODE-RNA-dynamics-paper

This repository contains data generation and analysis code used in this paper.

----------------------------------------------------------------------------------------

1. **compute-stability.R**

    Calculates the relative stability of exonic RNA during the 0 to 2 hour and 2 to 6 hour chase periods. RNA stability is represented as Log2 fold change values, which are further scaled from 0 to 1 for each cell line after applying expression filters.
   - Input files: general_input_files.zip
   - Supplemental Table S13 provides the ENCODE accession numbers for the genic features quantifications files on the ENCODE portal (https://www.encodeproject.org/) which contain the exon_sense read counts used in this analysis.

2. **01B_rna_categories_stability.R**

    Determines the relative stability of intronic RNA, enhancer RNA (eRNA), PROMPTs and readthrough (RT) RNA, using exonic RNA as the reference. This analysis is related to Main Figure 1B.
   - Input file: Supplemental Table S2

4. **01D_bru026_SI_plots.R**

    Generates plots of splicing indices for common introns across 16 cell lines at 0, 2 and 6 hours, as displayed in Main Figure 1D.
   - Input files: input_files_Fig01D.zip

4. **02AB_gene_exonic_stability_correlation_heatmaps.R**

    Computes correlation coefficients among the 16 cell lines using exonic relative stability values (Log2 fold change). This analysis corresponds to Main Figures 2A and 2B.
   - Supplemental Table S5 lists the log2 fold change values for the two chase periods used in this analysis.

6. **03CD_suppFig4_scaledval_distribution.R**

    Plots the distribution of scaled RNA stability values for each cell line. The resulting plots can be found in Main Figures 3C and 3D, as well as Supplemental Figure 5.
   - Input file: Supplemental Table S7

7. **04F_026h_SIbased_exonsONLY_IsoPct_trends.R**

    Measures isoform percentages and their prevalence over time. Findings are shown in Main Figure 4F.
   - Input files: input_files_Fig04F_1.zip, input_files_Fig04F_2.zip, input_files_Fig04F_3.zip, input_files_Fig04F_4.zip, input_files_Fig04F_5.zip

8. **06BD_fraction_splicing_category**

    Calculates the fractions of four splicing pattern categories per sample, then generates the stacked bar plot for each of the three timepoints. The 0h plot is shown in Main Figure 6B. This script also creates Main Figure 6D, which shows the mean and standard deviation of the fractions for each pattern at the three time points.
   - Input files: Supplemental Tables S10, S11 and S12

8. **suppFig_12and13.R**

    Defines thresholds for classifying introns as "Stable" or "Retained", and generates the plots shown in Supplemental Figures 12 and 13.
   - Input files: Supplemental Table S8, sheet "6vs0".

