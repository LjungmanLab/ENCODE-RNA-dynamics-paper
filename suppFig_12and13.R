#############################################################################################################

## This code was used to generate Supplemental Figures 12 and 13, and Supplemental Table S9              ## 
## Input file is Supplemental Table S8 Relative intron stability and SI values.xlsx - 6vs0.tsv             ##                                               

#############################################################################################################

library(tidyverse)

# Read the input file
d0<-read.delim("Table S8 Relative intron stability and SI values.xlsx - 6vs0.tsv", 
               h= TRUE, stringsAsFactors = FALSE, sep = "\t", na.strings = c(""," ","NA"), check.names = F)

# Create a working dataset filtered by isoform count and intron status
d1<-d0 %>% filter(isoform_count == 1 & is_TRUE_intron == "Yes")


# Subset the required metadata and 6h Vs 0h relative fold change columns
dfc<- d1 %>% mutate(ID = paste(gene,intron_position, sep="_"), .before = 1) %>% select(1,12:27)

# Convert the table into a long format and generate quartiles for the FC values for each cell line
dfc1 <- dfc %>% pivot_longer(cols = -ID, names_to = "sample_FC", values_to = "6vs0_FC") %>%
  mutate(cell_line = substr(sample_FC, 1, nchar(sample_FC) - 8)) %>% group_by(cell_line) %>%
  mutate(quartile_FC = ntile(`6vs0_FC`,4)) %>% ungroup()

# Subset the required metadata and 6h splicing index columns
dsi<- d1 %>% mutate(ID = paste(gene,intron_position, sep="_"), .before = 1) %>% select(1,34:49)

# Convert the table into a long format and generate quartiles for the SI values for each cell line
dsi1<- dsi %>% pivot_longer(cols = -ID, names_to = "sample_SI", values_to = "6h_SI") %>%
  mutate(cell_line = substr(sample_SI, 1, nchar(sample_SI) - 6)) %>% group_by(cell_line) %>%
  mutate(quartile_SI = ntile(`6h_SI`,4)) %>% ungroup()

# Merge the long format tables
# Perform statistics on the values
# Assign thresholds to classify an intron as STABLE or RETAINED
df<-merge(dfc1,dsi1,by=c("ID","cell_line")) %>% group_by(cell_line) %>%
  mutate(mean_6vs0_FC = mean(`6vs0_FC`, na.rm = TRUE),
    median_6vs0_FC = median(`6vs0_FC`, na.rm = TRUE),
    sd_6vs0_FC   = sd(`6vs0_FC`, na.rm = TRUE),
    z_6vs0_FC    = (`6vs0_FC` - mean_6vs0_FC) / sd_6vs0_FC,
    mean_6h_SI   = mean(`6h_SI`, na.rm = TRUE),
    median_6h_SI = median(`6h_SI`, na.rm = TRUE),
    sd_6h_SI     = sd(`6h_SI`, na.rm = TRUE),
    z_6h_SI      = (`6h_SI` - mean_6h_SI) / sd_6h_SI) %>%
  ungroup() %>% drop_na() %>%
  mutate(IS_intron_retained = (`z_6vs0_FC` > 3 & `z_6h_SI` < -3)) %>%
  mutate(IS_intron_stable = (`z_6vs0_FC` > 3 & `6h_SI` > mean_6h_SI))

# Count how many cell lines have common introns
df.summary.retained <- df %>% group_by(ID) %>%
  summarize(n_true = sum(IS_intron_retained, na.rm = TRUE)) %>%
  count(n_true) %>% filter(n_true!= 0)

df.summary.stable <- df %>% group_by(ID) %>%
  summarize(n_true = sum(IS_intron_stable, na.rm = TRUE)) %>%
  count(n_true) %>% filter(n_true!= 0)


# Compile the final supplementary table with additional metadata columns
d2<-d1 %>% select(1:10,77:ncol(d1)) %>% mutate(ID = paste(gene,intron_position, sep="_"))

df.final<-merge(d2,df, by=c("ID")) %>% select(2:17,19,20,22:last_col())

write.table(df.final,paste0("Derived_from_SuppTable_S8"),
            sep = "\t", col.names = T, row.names = F, quote = F)



## Generate Plots 

# Scatter plot - Retained introns

p1<- ggplot(df, aes(x = `6vs0_FC`, y = `6h_SI`)) +
  geom_point(aes(color = IS_intron_retained), size=1) +
  scale_color_manual(values = c('FALSE' = "#B3B3B3",
                                'TRUE'  = "#36648B")) +
  facet_wrap(~cell_line) +
  labs(x = "6vs0_FC", y = "6h_SI", color = "Retained introns") + 
  theme_minimal() +
  guides(color="none")+
  ggtitle("Retained Introns")+
  theme(axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold"),
    strip.text = element_text(face = "bold"))+
  geom_hline(yintercept = 0, linetype = "solid", color = "#000000", linewidth=0.1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "#000000", linewidth=0.1)


ggsave("retained_intron_dist_1.png", 
       plot = p1 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)

ggsave("retained_intron_dist_1.pdf", 
       plot = p1 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)


# Scatter plot - Stable introns

p2<- ggplot(df, aes(x = `6vs0_FC`, y = `6h_SI`)) +
  geom_point(aes(color = IS_intron_stable), size=1) +
  scale_color_manual(values = c('FALSE' = "#B3B3B3",
                                'TRUE'  = "#8B3A3A")) +
  facet_wrap(~cell_line) +
  labs(x = "6vs0_FC", y = "6h_SI", color = "Stable introns") +
  theme_minimal()+
  guides(color="none")+
  ggtitle("Stable Introns")+
  theme(axis.title = element_text(face = "bold"),
    axis.text  = element_text(face = "bold"),
    strip.text = element_text(face = "bold"))+
  geom_hline(yintercept = 0, linetype = "solid", color = "#000000", linewidth=0.1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "#000000", linewidth=0.1)


ggsave("stable_intron_dist_1.png", 
       plot = p2 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)

ggsave("stable_intron_dist_1.pdf", 
       plot = p2 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)


# Bar plot - Retained introns

p3<- ggplot(df.summary.retained, aes(x = factor(n_true), y = n)) +
  geom_col(fill = "#36648B") +
  geom_text(aes(label = n), vjust = -0.5, fontface = "bold") +
  labs(x = "Number of cell lines with matching retained introns",
    y = "Number of Introns",
    title = "Distribution of retained introns across cell lines") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"))+
  geom_hline(yintercept = 0, linetype = "solid", color = "#000000", linewidth=0.1)


ggsave("retained_intron_count_1.png", 
       plot = p3 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)

ggsave("retained_intron_count_1.pdf", 
       plot = p3 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)


# Bar plot - Stable introns

p4<- ggplot(df.summary.stable, aes(x = factor(n_true), y = n)) +
  geom_col(fill = "#8B3A3A") +
  geom_text(aes(label = n), vjust = -0.5, fontface = "bold") +
  labs(x = "Number of cell lines with matching stable introns",
    y = "Number of Introns",
    title = "Distribution of stable introns across cell lines") +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold"))+
  geom_hline(yintercept = 0, linetype = "solid", color = "#000000", linewidth=0.1)


ggsave("stable_intron_count_1.png", 
       plot = p4 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)

ggsave("stable_intron_count_1.pdf", 
       plot = p4 , 
       width = 8, 
       height = 6,
       bg="#FFFFFF",
       dpi = 300)





