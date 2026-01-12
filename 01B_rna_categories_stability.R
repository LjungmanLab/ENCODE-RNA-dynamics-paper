# load libraries
library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)


# load counts file 

all_feats <- read_excel("all_features_counts_84samples.xlsx")

norm_feats_counts <- all_feats

#### Exon ####

# normalizing to nMapUniq
norm_feats_counts$exon_nmap_norm_cts <- (all_feats$exon_total_counts/(all_feats$nMapUniq/1e6))


#### Intron ####

# normalizing to nMapUniq
norm_feats_counts$intron_nmap_norm_cts <- (all_feats$intron_total_counts/(all_feats$nMapUniq/1e6))

#### PROMPTS ####

# normalizing to nMapUniq
norm_feats_counts$prompts_nmap_norm_cts <- (all_feats$prompts_total_counts/(all_feats$nMapUniq/1e6))

#### eRNA ####

# normalizing to nMapUniq
norm_feats_counts$erna_nmap_norm_cts <- (all_feats$erna_total_counts/(all_feats$nMapUniq/1e6))

#### RT ####

# normalizing to nMapUniq
norm_feats_counts$readthrough_nmap_norm_cts <- (all_feats$readthrough_total_counts/(all_feats$nMapUniq/1e6))

# group by assay and get cumulative normalized counts 
all_feats_norm_grouped <- norm_feats_counts %>% 
  group_by(assay) %>% 
  summarise(exon_sum_nmap_norm_cts = sum(exon_nmap_norm_cts), intron_sum_nmap_norm_cts = sum(intron_nmap_norm_cts),
            prompts_sum_nmap_norm_cts = sum(prompts_nmap_norm_cts), erna_sum_nmap_norm_cts = sum(erna_nmap_norm_cts),
            readthrough_sum_nmap_norm_cts = sum(readthrough_nmap_norm_cts))


#normalize to 0hr counts

# store names of columns in a variable
cols <- colnames(all_feats_norm_grouped[-1])


# store 0h value of each column to divide each timepoint by
zero_vals <- all_feats_norm_grouped[1, cols]

# divide each timepoiont for each feature by the 0h normalized sum then make a new column for those values
all_feats_norm_grouped <- all_feats_norm_grouped %>%
  mutate(across(all_of(cols), ~ .x / zero_vals[[cur_column()]],
                .names = "{.col}_0hr"))

# store names of the new columns in a variable 
norm_cols_nmap <- grep("nmap_norm_cts_0hr$", names(all_feats_norm_grouped), value = TRUE)

# divide each features 0h normal counts by its corresponding timepoint's exon_sum_nmap_norm_cts_0h
all_feats_norm_grouped <- all_feats_norm_grouped %>%
  mutate(across(all_of(norm_cols_nmap),
                ~ .x / exon_sum_nmap_norm_cts_0hr,
                .names = "{.col}_exon"))

###### steps for making 1b plot

# long formed table

long <- all_feats_norm_grouped %>% 
  select(c(1,12:16)) %>% 
  pivot_longer(cols = c(2:6),
               names_to = "features",
               values_to = "exon_norm" )

#take out 0hr
only_2_6 <- long[long$assay != "Bru-seq",]

# PLOT 

only_2_6_plot <- ggplot() +
  geom_col(data = only_2_6, aes(x = features, y = exon_norm, fill = assay), position = "dodge") +
  labs(title = "RNA Stability compared to Exons",
       x = "RNA Categories", y = "Relative Stability") +
  scale_x_discrete(limits = c("exon_sum_nmap_norm_cts_0hr_exon", "intron_sum_nmap_norm_cts_0hr_exon", 
                              "prompts_sum_nmap_norm_cts_0hr_exon", "erna_sum_nmap_norm_cts_0hr_exon",
                              "readthrough_sum_nmap_norm_cts_0hr_exon"), 
                   labels = c("exon_sum_nmap_norm_cts_0hr_exon" = "Exon", "intron_sum_nmap_norm_cts_0hr_exon" = "Intron", 
                              "prompts_sum_nmap_norm_cts_0hr_exon" = "PROMPTs", "erna_sum_nmap_norm_cts_0hr_exon" = "eRNA",
                              "readthrough_sum_nmap_norm_cts_0hr_exon" = "Readthrough")) +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = c("BruChase-seq_2h" = "#4F94CD", "BruChase-seq_6h" = "#EE0000"), 
                    labels = c("BruChase-seq_2h"  = "2h", "BruChase-seq_6h" = "6h"))+
  theme_bw()+
  theme(legend.position = c(.93, .83),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(.5, "cm"), 
        legend.background = element_rect(fill = "transparent", colour = "transparent"))

ggsave(
  filename = "only_2_6_exon_norm_plot.png", 
  plot = only_2_6_plot,
  width = 6, 
  height = 2, 
  units = "in", 
  dpi = 600
)
