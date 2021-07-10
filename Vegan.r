#Test correlations otus X environmental variables
#### Only using Vegan:
library(vegan)
info_data <- read.csv("C:/Gresham/Project_Gravimondo/Project_Impact_2/Supplemental/unassigned_csv/level-7_TSS.csv", row.names=1)

info_data <- read.csv("C:/Gresham/Project_Gravimondo/Project_Impact_2/PeerJ_Submission/Otu_frequency_TSS.csv", row.names=1)
Environmental_variables <- read.delim("C:/Gresham/Project_Gravimondo/Project_Impact_2/metadata/env_variables_expanded.txt", row.names = 1)

Otu_veg <- t(info_data) # Samples become rows and OTUs are columns
Otus.pca <- rda(Otu_veg)
OTUs_dist <- vegdist(Otu_veg) # Creates distance matrix of the OTUs
#Plot PCA
PCA_plot<-ordiplot(Otus.pca)
PCA_plot
# Plotting environmental variables
Envdt <- Environmental_variables
#Envdt <- Environmental_variables[,2:20]
Tax_env <- envfit(OTUs_dist, Envdt, permu=999)
Tax_env
plot(Tax_env, p.max = 0.1, col = "blue")
plot(Tax_env, p.max = 0.05, col = "red")

# add labels to the samples in the ordiplot
orditorp(PCA_plot, "site", pch="+", pcol="grey")


# ----
## Bray-Curtis distances between samples
taxa_counts <- read.delim("C:/Gresham/Project_Gravimondo/Project_Impact_2/qiime_results/asv2taxa_feature_taxonomy.tab", sep='\t')
my_varespec = taxa_counts
rownames(my_varespec) <- my_varespec$taxa
my_varespec <- my_varespec[-c(1)]

my_varespec <- as.data.frame(t(my_varespec))
my_dis <- vegdist(my_varespec)

## groups 3 and 3
my_groups <- factor(c(rep(1,3), rep(2,3)), labels = c("Pristine","Impacted"))

## Calculate multivariate dispersions
my_mod <- betadisper(my_dis, my_groups)

## Perform test
anova(my_mod)