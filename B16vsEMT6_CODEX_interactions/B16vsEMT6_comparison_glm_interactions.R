################################################################################################################################
################################################################################################################
########################################################################################
######## Comparison of the two CODEX runs by Nicole for the Ad-IL12 story. The idea is to see if IL12 has different effects in the two tumor models. Ideally, it has much more effects
######## in the EMT6 model, where the effect on tumor growth and survival are much bigger than in B16. The basis of this analysis is the glm negative binomial model run in separate R files
######## here I just load these results and perform some plots
######## Marcel Trefny 11.5.21
################################################################################################
########################################################################################################
########################################################################################################################

#import

library(viridis)
library(devtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(circlize)
library(ComplexHeatmap)
library(tidyr)
library(edgeR)
library(multcomp)
library(gridExtra)
library(RColorBrewer)
library(forcats)
library(readxl)
library(pscl)
set.seed(123) #define randomness generation 
#initial definitions

myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

#make custom functions
#custom function to write tab deliminated text with the correct settings
write.table.tab <- function(x, file = ""){
  write.table(x,file,sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
}

### Import Data files
path <- "/Volumes/CIMM$/CIMM/Projects_Human_Team/Marcel_Nicole/AdV5_IL12_Project/B16vsEMT6_CODEX_interactions_Marcel/" #change to file path
setwd(path)

files <- list.files(pattern = "glm_neg_binomial")
files

#clean and combine the two lists of statistics
interaction_statistics <- sapply(files, function(x) read.table(x,  sep = "\t", header = TRUE), simplify = FALSE, USE.NAMES = TRUE)
interaction_statistics

names(interaction_statistics) <- gsub(".txt", "", gsub("glm_neg_binomial_interactions_supercluster_", "",names(interaction_statistics) ))
interaction_statistics <-  sapply(names(interaction_statistics), function(x) {interaction_statistics[[x]]$tumortype <- x; interaction_statistics[[x]] }, simplify = FALSE, USE.NAMES = TRUE  )

interaction_statistics_combined <- inner_join(interaction_statistics[[1]], interaction_statistics[[2]], by = c("comparison", "condition"), suffix = paste0("_",names(interaction_statistics)) )
#also calculate log fold changes from estimate (is in exponential space)
interaction_statistics_combined$lfc_B16 <- log2(exp(interaction_statistics_combined$estimate_B16))
interaction_statistics_combined$lfc_EMT6 <- log2(exp(interaction_statistics_combined$estimate_EMT6))

if (!dir.exists("./lists")){dir.create("./lists/")}
write.table( interaction_statistics_combined, file = "./lists/interaction_statistics_combined_EMT6_B16.txt", sep = "\t", quote = FALSE, row.names = FALSE)



#define a function to plot z-values on axis and log fold changes and p values as other elements
plot_subdata_z<- function(term = ""){
subdata <- interaction_statistics_combined %>% filter(grepl(term, comparison))
subdata
p <- ggplot(subdata, aes( x= zvalue_B16,y=zvalue_EMT6, label = comparison,size = abs(lfc_EMT6), color =  -log10(pvalue_EMT6) ))  +
  ggtitle(paste0(" Interaction Analysis --  ", term)) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("B16 z-value AdV5-IL12 vs untreated") +
  ylab("EMT6 z-value AdV5-IL12 vs untreated") +
  geom_point()+
  geom_label_repel(data = subset(subdata, pvalue_EMT6 < 0.05 ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "EMT6 -log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log2 fold change\nAdV5-IL12 vs untreated in EMT6") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p
}


plot_subdata_OR <- function(term = ""){
  subdata <- interaction_statistics_combined %>% filter(grepl(term, comparison))
  subdata
  p <- ggplot(subdata, aes( x= log10(OR_B16),y=log10(OR_EMT6), label = comparison,size = abs(lfc_EMT6), color =  -log10(pvalue_EMT6) ))  +
    ggtitle(paste0(" Interaction Analysis --  ", term)) +
    geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
    geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
    geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
    xlab("B16 log10 Odds Ratio AdV5-IL12 vs untreated") +
    ylab("EMT6 log10 Odds Ratio AdV5-IL12 vs untreated") +
    geom_point()+
    geom_label_repel(data = subset(subdata, pvalue_EMT6 < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
    scale_color_gradient( low="blue",high="red", space ="Lab" , name = "EMT6 -log10 p-value\nAdV5-IL12 vs untreated") +
    #scale_color_manual(values=c("#57BFFA", "#DE2902"))
    scale_size(name="abs log2 fold change\nAdV5-IL12 vs untreated in EMT6") +
    theme_bw() + 
    theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
  p
}

pdf(file ="./plots/NegBinomial_Models_Comparisons_EMT6_vs_B16.pdf", width = 8, height = 6)
#select different subsets to plot
plot_subdata_OR("") #all comparisons
plot_subdata_OR("CD8")
plot_subdata_OR("NK")
plot_subdata_OR("DC")
plot_subdata_OR("CD8|NK")

#select different subsets to plot
plot_subdata_z("") #all comparisons
plot_subdata_z("CD8")
plot_subdata_z("NK")
plot_subdata_z("DC")
plot_subdata_z("CD8|NK") ## This plot was used as Figure 3F showing the differneces in cell interactions of CD8 and NK cells with other cell types
dev.off()
