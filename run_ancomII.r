#module load r/intel/4.0.4
#R
install.packages("nlme")
#78
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("compositions")
install.packages("readr")

#setwd('C:/Gresham/Project_Gravimondo/Project_Impact_2/PeerJ_Submission/ANCOM-master/ANCOM-master')
setwd('/scratch/ps163/Dr_Carolina/Project_Impact/ANCOM/')
library(nlme)

library(tidyverse)
library(ggplot2)
library(compositions)
source("scripts/ancom_v2.1.R")

library(readr)

otu_data = read_tsv("Impact/feature-table.biom.txt")
otu_id = otu_data$`#OTU_ID`
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

meta_data = read_tsv("Impact/MappingFile_Impact.csv")
meta_data = meta_data %>% rename(Sample.ID = 'sample-id')

source("scripts/ancom_v2.1.R")

# Step 1: Data preprocessing

feature_table = otu_data; sample_var = "Sample.ID"; group_var = NULL
out_cut = 0.25; zero_cut = 0.90; lib_cut = 100; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, meta_data, sample_var, group_var, 
                                   out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM

main_var = "Site"; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = NULL
t_start = Sys.time()
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, 
            alpha, adj_formula, rand_formula)
t_end = Sys.time()
t_run = t_end - t_start # around 30s

write_csv(res$out, "outputs/res_impact_3.csv")

# Step 3: Volcano Plot

# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")

fig = res$fig +  
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") + 
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE)
fig  
