library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)

# load input tables
table_0h <- read_tsv("Supplemental_Table_S10_0h.tsv")
table_2h <- read_tsv("Supplemental_Table_S11_2h.tsv")
table_6h <- read_tsv("Supplemental_Table_S12_6h.tsv")


#### 0hr 

# filter for the four columns we are interested in 
# find total of counts for each column
# make table into long form to make later steps easier
# split 'Splicing Pattern" column to get cell line and pattern information separately 
table_0_sum <- table_0h %>% 
  select(matches("US_SPL_DS_UNSPL|US_UNSPL_DS_SPL|USDS_SPL|USDS_UNSPL")) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>% 
  pivot_longer(cols = 1:64, names_to = "Splicing Pattern", values_to = "Sum") %>% 
  mutate(cellline = str_extract(`Splicing Pattern`, "^[^_]+")) %>% 
  mutate(pattern = sub("^[^_]+_", "", `Splicing Pattern`))

# group by cell line and find total per cell line
total_0h <- table_0_sum %>% 
  group_by(cellline) %>% 
  summarise(total = sum(Sum, na.rm = TRUE))

# divide each total per pattern by the total per corresponding cell line
table_0_sum$fraction<- apply(table_0_sum, 1, function(x){
  
  cl <- x["cellline"]
  sum <- as.numeric(x["Sum"])
  total <- as.numeric(total_0h$total[total_0h == cl])
  fraction <- sum/total
  
  return(fraction)
})

# plot
fraction_0_plot <- ggplot()+
  geom_col(data = table_0_sum, aes(x = cellline, y = fraction, fill = pattern), position = position_stack(reverse = TRUE), width = .7) +
  labs(title = "Fraction Spliced 0hr", x = element_blank(), y = "splicing category fraction") +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = c("US_SPL_DS_UNSPL" = "#CD0000", "US_UNSPL_DS_SPL" = "#EEC900", "USDS_SPL" = "#66CD00","USDS_UNSPL" = "#4F94CD" ),
                    labels = c("US_SPL_DS_UNSPL" = "S-U", "US_UNSPL_DS_SPL" = "U-S", "USDS_SPL" = "S-S","USDS_UNSPL" = "U-U"))+
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(
  filename = "fig6_b_0h_plot.png", 
  plot = fraction_0_plot,
  width = 6, 
  height = 3, 
  units = "in", 
  dpi = 600
)

### 2hr

# filter for the four columns we are interested in 
# find total of counts for each column
# make table into long form to make later steps easier
# split 'Splicing Pattern" column to get cell line and pattern information separately 
table_2_sum <- table_2h %>% 
  select(matches("US_SPL_DS_UNSPL|US_UNSPL_DS_SPL|USDS_SPL|USDS_UNSPL")) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>% 
  pivot_longer(cols = 1:64, names_to = "Splicing Pattern", values_to = "Sum") %>% 
  mutate(cellline = str_extract(`Splicing Pattern`, "^[^_]+")) %>% 
  mutate(pattern = sub("^[^_]+_", "", `Splicing Pattern`))

# group by cell line and find total per cell line
total_2h <- table_2_sum %>% 
  group_by(cellline) %>% 
  summarise(total = sum(Sum, na.rm = TRUE))

# divide each total per pattern by the total per corresponding cell line
table_2_sum$fraction<- apply(table_2_sum, 1, function(x){
  
  cl <- x["cellline"]
  sum <- as.numeric(x["Sum"])
  total <- as.numeric(total_2h$total[total_2h == cl])
  fraction <- sum/total
  
  return(fraction)
})

# plot
fraction_2_plot <- ggplot()+
  geom_col(data = table_2_sum, aes(x = cellline, y = fraction, fill = pattern), position = position_stack(reverse = TRUE)) +
  labs(title = "Fraction Spliced 2hr", x = element_blank(), y = "splicing category fraction") +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = c("US_SPL_DS_UNSPL" = "#CD0000", "US_UNSPL_DS_SPL" = "#EEC900", "USDS_SPL" = "#66CD00","USDS_UNSPL" = "#4F94CD" ),
                    labels = c("US_SPL_DS_UNSPL" = "S-U", "US_UNSPL_DS_SPL" = "U-S", "USDS_SPL" = "S-S","USDS_UNSPL" = "U-U"))+
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(
  filename = "fig6_b_2h_plot.png", 
  plot = fraction_2_plot,
  width = 6, 
  height = 3, 
  units = "in", 
  dpi = 600
)

####### 6hr

# filter for the four columns we are interested in 
# find total of counts for each column
# make table into long form to make later steps easier
# split 'Splicing Pattern" column to get cell line and pattern information separately 
table_6_sum <- table_6h %>% 
  select(matches("US_SPL_DS_UNSPL|US_UNSPL_DS_SPL|USDS_SPL|USDS_UNSPL")) %>% 
  summarise(across(everything(), sum, na.rm = TRUE)) %>% 
  pivot_longer(cols = 1:64, names_to = "Splicing Pattern", values_to = "Sum") %>% 
  mutate(cellline = str_extract(`Splicing Pattern`, "^[^_]+")) %>% 
  mutate(pattern = sub("^[^_]+_", "", `Splicing Pattern`))

# group by cell line and find total per cell line
total_6h <- table_6_sum %>% 
  group_by(cellline) %>% 
  summarise(total = sum(Sum, na.rm = TRUE))

# divide each total per pattern by the total per corresponding cell line
table_6_sum$fraction<- apply(table_6_sum, 1, function(x){
  
  cl <- x["cellline"]
  sum <- as.numeric(x["Sum"])
  total <- as.numeric(total_6h$total[total_6h == cl])
  fraction <- sum/total
  
  return(fraction)
})

# plot
fraction_6_plot <- ggplot()+
  geom_col(data = table_6_sum, aes(x = cellline, y = fraction, fill = pattern), position = position_stack(reverse = TRUE)) +
  labs(title = "Fraction Spliced 6hr", x = element_blank(), y = "splicing category fraction") +
  scale_y_continuous(expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = c("US_SPL_DS_UNSPL" = "#CD0000", "US_UNSPL_DS_SPL" = "#EEC900", "USDS_SPL" = "#66CD00","USDS_UNSPL" = "#4F94CD" ),
                    labels = c("US_SPL_DS_UNSPL" = "S-U", "US_UNSPL_DS_SPL" = "U-S", "USDS_SPL" = "S-S","USDS_UNSPL" = "U-U"))+
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(
  filename = "fig6_b_6h_plot.png", 
  plot = fraction_6_plot,
  width = 6, 
  height = 3, 
  units = "in", 
  dpi = 600
)

####### fig 6d

# group by pattern 
grouped_0 <- table_0_sum %>% 
  group_by(pattern) %>% 
  summarise(avg = mean(fraction), SD = sd(fraction))

# add timepoint information 
grouped_0$time <- "0h"

# group by pattern 
grouped_2 <- table_2_sum %>% 
  group_by(pattern) %>% 
  summarise(avg = mean(fraction), SD = sd(fraction))

# add timepoint information 
grouped_2$time <- "2h"

# group by pattern 
grouped_6 <- table_6_sum %>% 
  group_by(pattern) %>% 
  summarise(avg = mean(fraction), SD = sd(fraction))

# add timepoint information 
grouped_6$time <- "6h"

# make one data table that includes all three timepoints
all <- bind_rows(grouped_0, grouped_2, grouped_6)

# plot
plot_6D <- ggplot(all, aes(x = pattern, y = avg, fill = time)) +
  geom_col(position = position_dodge(), width = .7) +
  labs(title = "Fraction Spliced", x = element_blank(), y = "splicing category fraction") +
  geom_errorbar(aes( ymin = avg - SD, ymax = avg + SD), 
                width = 0.2,
                position = position_dodge(width = 0.9)) +
  scale_x_discrete(labels = c("US_SPL_DS_UNSPL" = "S-U", "US_UNSPL_DS_SPL" = "U-S", "USDS_SPL" = "S-S","USDS_UNSPL" = "U-U")) +
  scale_y_continuous(expand = expansion(mult = 0)) +
  scale_fill_manual(values = c("0h" = "limegreen", "2h" = "steelblue2", "6h" = "orangered" ))+
  theme_bw() +
  theme(legend.position = c(.15, .83),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent", colour = "transparent"))

ggsave(
  filename = "fig6_D_plot.png", 
  plot = plot_6D,
  width = 6, 
  height = 3, 
  units = "in", 
  dpi = 600
)
