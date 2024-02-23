library("tidyverse")
library("RColorBrewer")
library("dplyr")
library("ggpubr")
fill_colors <- c("steelblue","gray83", "skyblue3","lightcyan3", "gray60", "lightsteelblue2", "gray40", "skyblue4", "lightskyblue2", "gray83")

####import####
setwd("/Users/mangofrog7/Documents/Graduate_School/Projects/bioinformatics/TC_final/K9_peric_enr")

# Get the list of files in the directory
files <- list.files(pattern = "\\.tab$")

# Iterate through each file
for (file in files) {
  # Extract relevant information using a regular expression
  sat <- gsub(".*K9peaks_(.*)_peakmatrix.*", "\\1", file)
  
  # Read the file into a data frame
  df <- read_delim(file, delim = '\t', skip = 2)
  
  # Name the data frame based on the specified format (excluding "AVG" prefix)
  time <- as.numeric(str_extract(file, "\\d+\\.?\\d*"))
  df_name <- paste(sat, time, "K9", sep = "_")
  assign(df_name, df)
}

periC_only_2hpf <- peric_only_2hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2hpf")
periC_only_2.5hpf <- peric_only_2.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2.5hpf")
periC_only_3hpf <- peric_only_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
periC_only_3.5hpf <- peric_only_3.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3.5hpf")
periC_only_4hpf <- peric_only_4hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4hpf")
periC_only_4.5hpf <- peric_only_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

piRNA_only_2hpf <- piRNA_only_2hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2hpf")
piRNA_only_2.5hpf <- piRNA_only_2.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2.5hpf")
piRNA_only_3hpf <- piRNA_only_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
piRNA_only_3.5hpf <- piRNA_only_3.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3.5hpf")
piRNA_only_4hpf <- piRNA_only_4hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4hpf")
piRNA_only_4.5hpf <- piRNA_only_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

piRNA_N_peric_2hpf <- piRNA_N_peric_2hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2hpf")
piRNA_N_peric_2.5hpf <- piRNA_N_peric_2.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2.5hpf")
piRNA_N_peric_3hpf <- piRNA_N_peric_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
piRNA_N_peric_3.5hpf <- piRNA_N_peric_3.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3.5hpf")
piRNA_N_peric_4hpf <- piRNA_N_peric_4hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4hpf")
piRNA_N_peric_4.5hpf <- piRNA_N_peric_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

noperic_nopiRNA_2hpf <- noperic_nopiRNA_2hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2hpf")
noperic_nopiRNA_2.5hpf <- noperic_nopiRNA_2.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2.5hpf")
noperic_nopiRNA_3hpf <- noperic_nopiRNA_3hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
noperic_nopiRNA_3.5hpf <- noperic_nopiRNA_3.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3.5hpf")
noperic_nopiRNA_4hpf <- noperic_nopiRNA_4hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4hpf")
noperic_nopiRNA_4.5hpf <- noperic_nopiRNA_4.5hpf_9_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

#####
periC_only <- bind_rows(periC_only_2hpf, periC_only_2.5hpf, periC_only_3hpf, periC_only_3.5hpf, periC_only_4hpf, periC_only_4.5hpf)
periC_only_enr <- periC_only %>% select(time_point, total) %>% 
  mutate(time_point = factor(time_point, levels = c("2hpf", "2.5hpf", "3hpf", "3.5hpf", "4hpf", "4.5hpf"))) %>%
  mutate(cat = "Peric.")

piRNA_only <- bind_rows(piRNA_only_2hpf,piRNA_only_2.5hpf, piRNA_only_3hpf, piRNA_only_3.5hpf, piRNA_only_4hpf, piRNA_only_4.5hpf)
piRNA_only_enr <- piRNA_only %>% select(time_point, total) %>% 
  mutate(time_point = factor(time_point, levels = c("2hpf", "2.5hpf", "3hpf", "3.5hpf", "4hpf", "4.5hpf"))) %>%
  mutate(cat = "piRNA")

piRNA_N_peric <- bind_rows(piRNA_N_peric_2hpf,piRNA_N_peric_2.5hpf, piRNA_N_peric_3hpf, piRNA_N_peric_3.5hpf, piRNA_N_peric_4hpf, piRNA_N_peric_4.5hpf)
piRNA_N_peric_enr <- piRNA_N_peric %>% select(time_point, total) %>% 
  mutate(time_point = factor(time_point, levels = c("2hpf", "2.5hpf", "3hpf", "3.5hpf", "4hpf", "4.5hpf"))) %>%
  mutate(cat = "Peric. & piRNA")

noperic_nopiRNA <- bind_rows(noperic_nopiRNA_2hpf,noperic_nopiRNA_2.5hpf, noperic_nopiRNA_3hpf, noperic_nopiRNA_3.5hpf, noperic_nopiRNA_4hpf, noperic_nopiRNA_4.5hpf)
noperic_nopiRNA_enr <- noperic_nopiRNA %>% select(time_point, total) %>% 
  mutate(time_point = factor(time_point, levels = c("2hpf", "2.5hpf", "3hpf", "3.5hpf", "4hpf", "4.5hpf"))) %>%
  mutate(cat = "Neither")

total <- bind_rows(periC_only_enr, piRNA_only_enr, piRNA_N_peric_enr, noperic_nopiRNA_enr)

t.test <- compare_means(total ~ cat, p.adjust.method = "bonferroni", data=total, method = "t.test", group.by = "time_point")

total_early <- total %>% filter(time_point %in% c("2hpf", "2.5hpf", "3hpf")) %>%
  mutate(cat = factor(cat, levels = c("Neither", "Peric.", "piRNA", "Peric. & piRNA")))
  
t.test_early <- t.test %>% filter(time_point %in% c("2hpf", "2.5hpf", "3hpf")) %>% 
  mutate(y.position = c(330, 360, 320, 350, 340, 370, 
                        330, 360, 320, 350, 340, 370, 
                        330, 360, 320, 350, 340, 370)) %>%
  mutate(p.signif = ifelse(p < 0.05, "*", "ns"))

early_peak_groups_no_out <- ggplot(total_early, aes(x = fct_rev(cat), y = total)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cat), notch = TRUE) +
  labs(y = "H3K9me3 Enrichment") +
  theme_minimal() +  facet_wrap(~time_point, nrow = 3, strip.position = "right") +
  scale_fill_manual(values = fill_colors) + ylim(0,371) +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 8)) + coord_flip() +
  theme(axis.title.y = element_blank()) + 
  stat_pvalue_manual(t.test_early, label = "p.signif", size = 2, tip.length = 2) 
ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig5/early_peak_groups_noOut_new.pdf", plot = early_peak_groups_no_out, width = 4.5, height = 3.15, units = "in")

####going to add in the postEGA peaks to look at last three timepoints####
postEGA_2hpf <- postEGApeaks_comp_2hpf_peakmatrix.tab_2_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2hpf")
postEGA_2.5hpf <- postEGApeaks_comp_2.5hpf_peakmatrix.tab_2.5_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "2.5hpf")
postEGA_3hpf <- postEGApeaks_comp_3hpf_peakmatrix.tab_3_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3hpf")
postEGA_3.5hpf <- postEGApeaks_comp_3.5hpf_peakmatrix.tab_3.5_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "3.5hpf")
postEGA_4hpf <- postEGApeaks_comp_4hpf_peakmatrix.tab_4_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4hpf")
postEGA_4.5hpf <- postEGApeaks_comp_4.5hpf_peakmatrix.tab_4.5_K9 %>% rowwise() %>% mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(time_point = "4.5hpf")

postEGA_enr <- postEGA_4.5hpf %>% select(time_point, total) %>%  mutate(cat = "postEGA")

total2 <- bind_rows(postEGA_enr, total)
total2_late <- total2 %>% filter(time_point %in% c("3.5hpf", "4hpf", "4.5hpf")) %>%
  mutate(cat = factor(cat, levels = c("Neither", "Peric.", "piRNA", "Peric. & piRNA", "postEGA")))

t.test2 <- compare_means(total ~ cat, p.adjust.method = "bonferroni", data=total2, method = "t.test", group.by = "time_point")
t.test_4.5h <- t.test2 %>% filter(time_point %in% c("4.5hpf")) %>% 
  mutate(y.position = c(30000, 29000, 28000, 31000, 26000, 33000, 25000, 32000, 27000, 34000)) %>%
  mutate(p.signif = ifelse(p < 0.05, "*", "ns"))

late_peak_groups <- ggplot(filter(total2_late, time_point == "4.5hpf"), aes(x = fct_rev(cat), y = total)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cat), notch = TRUE) +
  labs(y = "H3K9me3 Enrichment") + ylim(0,34050) +
  theme_minimal() + facet_wrap(~time_point, strip.position = "right") +
  scale_fill_manual(values = fill_colors) +
  theme(legend.position = "none") + 
  theme(text = element_text(size = 8)) +
  theme(axis.title.y = element_blank()) + coord_flip() +
  stat_pvalue_manual(t.test_4.5h, label = "p.signif", size = 2, tip.length = 0)

ggsave("/Users/mangofrog7/Documents/Graduate_School/Papers/K9_Kinetics/Fig5/late_peak_groups_new.pdf", plot = late_peak_groups, width = 4.5, height = 1.5, units = "in")

