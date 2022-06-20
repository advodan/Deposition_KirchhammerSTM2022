# this is the analysis of the B16 Codex Muliple Imaging Interaction Counts from the Kirchhammer et al STM 2022 Manuscript
# Performed by Marcel P. Trefny in 2021
# forms parts of the analysis shown in Figure 1F

# ### Install necessary packages
# #only need to do this once, thus commented out
# install.packages("circlize")
# install.packages("ggplot2")
# install.packages("devtools")
# install.packages("multcomp")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("flowCore")
# BiocManager::install("ggcyto")
# BiocManager::install("edgeR")
# BiocManager::install("ggrepel")
# BiocManager::install("gridExtra")
# BiocManager::install("dplyr")
# BiocManager::install("purrr")
# BiocManager::install("plyr")
# BiocManager::install("tidyr")
# BiocManager::install("reshape2")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("forcats")
# BiocManager::install("tidyverse")
# install.packages("viridis")
#install.packages("readxl")

# Load Packages
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


myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

#make custom functions
#custom function to write tab deliminated text with the correct settings
write.table.tab <- function(x, file = ""){
  write.table(x,file,sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
}

### Import Data files

path <- "/Volumes/CIMM$/CIMM/Projects_Human_Team/Marcel_Nicole/AdV5_IL12_Project/B16_CODEX_interactions_Marcel/InteractionCount/" #set to path of files


setwd(path)
interaction_count_files <- list.files(pattern = "CytoMAP")
interaction_count_files

### extract interaction counts in list format
interactions_counts <- sapply(interaction_count_files, function(x) data.matrix(read_excel(x, col_names = FALSE)), simplify = FALSE, USE.NAMES = TRUE)
names(interactions_counts) <- gsub("_ClusteringA.xls" , "", names(interactions_counts))
names(interactions_counts) <- gsub("CytoMAP_Sample_" , "", names(interactions_counts))
head(interactions_counts)

headers <- read_excel("Headings.xls") #import information on rows and columns
headers
column_names <- gsub("DistTo_All_", "", colnames(headers)) #change to match row names
column_names <- column_names[2:length(column_names)] #get rid off the row names column
row_names <-  headers$...1  #rownames are the first column
row_names <- gsub(" ", "_" , row_names) #row names contained spaces -> change to underline

interactions_counts <- sapply(interactions_counts, function(x) {colnames(x) <-column_names ; x}, simplify = FALSE, USE.NAMES = TRUE)
interactions_counts <- sapply(interactions_counts, function(x) {rownames(x) <-row_names ; x}, simplify = FALSE, USE.NAMES = TRUE)

#columns and rows are not in order, reorder 
interactions_counts <- sapply(interactions_counts, function(x) { x <- x[,row_names]}, simplify = FALSE, USE.NAMES = TRUE)
head(interactions_counts)

# are column names and row names identical? should be all true if correct
sapply(interactions_counts, function(x) {summary(colnames(x) == rownames(x))}, simplify = TRUE) #all true

#now clean up data to have all clusters in all samples
clusters = unique(colnames(interactions_counts[[1]])) #names of clusters

#this is here just in case not all the files have all clusters. Will generate missing clusters and add 0 counts to this one
for (i in 1:length(interactions_counts)){ #for each sample
  for (j in 1:length(clusters)){ #for each cluster
    if (!clusters[j] %in% colnames(interactions_counts[[i]])){ #check if that cluster is in the dataset, if not
      a <- colnames(interactions_counts[[i]])
      b <- rownames(interactions_counts[[i]])
      interactions_counts[[i]] <- cbind(interactions_counts[[i]],  rep(0,nrow(interactions_counts[[i]])))  # add this cluster with all zeros to both columns
      interactions_counts[[i]] <- rbind(interactions_counts[[i]], rep(0,ncol(interactions_counts[[i]]))) # and rows
      colnames(interactions_counts[[i]]) <- c(a, clusters[j]) # clean up column names (did not find a better way... is a bit ugly)
      rownames(interactions_counts[[i]]) <- c(b, clusters[j])
    }
  }
}
#now there should be all true and equal number of values
sapply(interactions_counts, function(x) {summary(colnames(x) == rownames(x))}, simplify = TRUE) #68 true all

#sort all columns and rows so that each dataframe looks the same, is important for duplication removal during listing of comparisons
interactions_counts <- lapply(interactions_counts, function(df){
  df[order(rownames(df)),order(colnames(df))]
})
head(interactions_counts)


#define which clusters to exclude
#beware they need to be excluded again for the the superclusters
exclusion_clusters <- c("") #nothing excluded
# exclusion_clusters <- c("L")
#exclude these clusters
interactions_counts <- sapply(interactions_counts, function(x) {a <- colnames(x) == exclusion_clusters; b <- rownames(x) == exclusion_clusters; x[!b,!a] } , simplify = FALSE) #68 true all

interactions_counts

###############################################################################################
#### adpot nomenclature from the EMT6 Codex analysis. Therefore rename and summarize clusters to match the clusters in that experiment
########################################################################################################
supercluster_mapping <- read.table("../supercluster_annotation.txt", header = TRUE , sep = "\t")
supercluster_mapping
#also exclude the clusters from above here:
supercluster_mapping <- supercluster_mapping[!(supercluster_mapping$cluster == exclusion_clusters),]
#which are all unique superclusters
supercluster_mapping_unique <- unique(supercluster_mapping$supercluster)


interactions_counts_superclusters <- list()
#make the interaction_counts_supercluster from cluster data by summing the individual cluster
for (h in 1:length(interactions_counts)){  #for each sample
  #setup the new dataframe as an empty = 0 dataframe with the right columns and rows
  interactions_counts_superclusters[[h]] <- matrix(0, nrow = length(supercluster_mapping_unique), ncol=length(supercluster_mapping_unique))
  colnames(interactions_counts_superclusters[[h]]) <- supercluster_mapping_unique
  rownames(interactions_counts_superclusters[[h]]) <- supercluster_mapping_unique
  
  for (i in 1:length(supercluster_mapping_unique)){ # for each supercluster, find the columns and rows of the clusters
    for (j in 1:length(supercluster_mapping_unique)){
      #now find submatrix
      x <- supercluster_mapping_unique[i]
      y <- supercluster_mapping_unique[j]
      clusters_rows <- rownames(interactions_counts[[h]]) %in% supercluster_mapping$cluster[supercluster_mapping$supercluster == x] #which rows are needed
      clusters_cols <- colnames(interactions_counts[[h]]) %in% supercluster_mapping$cluster[supercluster_mapping$supercluster == y] # which columns are needed
      interactions_counts_superclusters[[h]][as.character(x),as.character(y)] <- sum(matrix(interactions_counts[[h]][clusters_rows, clusters_cols ] )) #make the sum of this submatrix of rows and columns
    }
  }
  
}
names(interactions_counts_superclusters) <- names(interactions_counts)
interactions_counts_superclusters

####################
#Overwrite the clusters according to the supercluster mapping to reach the same clusters as in EMT6 data
###############
interactions_counts <- interactions_counts_superclusters


#load the mapping of which condition is which sample
condition_mapping <- read.table("../Overview_MouseID_Region.csv", sep = ",", header = TRUE)
condition_mapping

#exclude regions
regions_keep <- condition_mapping %>% filter(Exclude == 0) %>% dplyr::select(Region) #define which regions to exclude
regions_keep <- regions_keep$Region
interactions_counts <- interactions_counts[ grepl(paste(regions_keep, collapse="|"), names(interactions_counts))]
head(interactions_counts)

#add multiple regions from one sample to each other
names(interactions_counts) %in% condition_mapping$Region
interactions_counts_added <- list()
mice <- unique(condition_mapping$MouseNr)
for(i in 1:length(mice)){
  mouse <- mice[i]
  regions <- condition_mapping %>% filter(MouseNr == mouse) %>% dplyr::select(Region)
  regions <- regions$Region
  interactions_counts_added[[as.character(mouse)]] <- Reduce("+", interactions_counts[names(interactions_counts) %in% regions])
  print(mouse)
}

head(interactions_counts_added)
length(names(interactions_counts_added))
#there are only superclusters this time, thus directly assign them
interactions_counts_superclusters <- interactions_counts_added #all names are correctly matching one region (all true)


#save files now with all columns and rows even if empty
if (!dir.exists("../interactions_counts_superclusters")){dir.create("../interactions_counts_superclusters/")}
sapply(names(interactions_counts_superclusters),function(x) write.table(interactions_counts_superclusters[[x]], file = paste0("../interactions_counts_superclusters/", gsub(".csv", "", x) , "_filled.txt"),sep = "\t", quote = FALSE))


########################################################################################################
#perform odds_ratio analysis

generate_log_odds_ratios <- function(interactions_counts_df){ #provide single interactions_counts/_superclusters data.frames and return their odds_ratios
  clusters <- unique(colnames(interactions_counts_df))
  x= 10^-4 #prior count to deal with 0 -> problem with logs, chosen to match Maths software data
  
  #reformat data into three columns using the melt function
  interaction_counts_melted <- setNames(melt(interactions_counts_df+x), c("cell1", "cell2", "interactions")) 
  
  #for each cluster count the total interactions for that cluster with all other clusters
  interactions_counts_sums <- c() #make emtpy vector to store the data
  for (i in 1:length(clusters)){
    interactions_counts_sums[clusters[i]] <- sum(subset(interaction_counts_melted, cell1 == clusters[i], select=interactions))  #important not to count them twice (thus only cell1)
  }
  total_interactions <- sum(interactions_counts_sums) #calculate total interactions of all cells with one another
  interactions_counts_sums <- data.frame(interactions_counts_sums, "interactions_counts_sums_frequency" = interactions_counts_sums/total_interactions) #calculate frequency of interations based on total interations
  
  expected_frequencies = interactions_counts_sums$interactions_counts_sums_frequency %o% interactions_counts_sums$interactions_counts_sums_frequency # from the vector above create its outer matrix product  = matrix with clusters and n columns and n rows, this now has
  expected_counts = expected_frequencies * total_interactions
  observed_frequencies = (interactions_counts_df+x)/total_interactions  #calculate observed frequencies
  odds_ratios <- observed_frequencies/expected_frequencies #odds ratio is just observed / expected
  colnames(odds_ratios) <- colnames(interactions_counts_df) #match colnames
  rownames(odds_ratios) <- rownames(interactions_counts_df) #match rownames
  return(log10(odds_ratios)) #return the new table of log_odds_ratios as the output of this function
}

log_odds_ratios_superclusters <- sapply(interactions_counts_superclusters, function(x) generate_log_odds_ratios(x), simplify = FALSE) #for each interactions_counts_superclusters file generate its odds ratios


generate_expected_counts <- function(interactions_counts_df){ #provide single interactions_counts/_superclusters data.frames and return their odds_ratios
  clusters <- unique(colnames(interactions_counts_df))
  x= 10^-4 #prior count to deal with 0 -> problem with logs, chosen to match Maths software data
  
  #reformat data into three columns using the melt function
  interaction_counts_melted <- setNames(melt(interactions_counts_df+x), c("cell1", "cell2", "interactions")) 
  
  #for each cluster count the total interactions for that cluster with all other clusters
  interactions_counts_sums <- c() #make emtpy vector to store the data
  for (i in 1:length(clusters)){
    interactions_counts_sums[clusters[i]] <- sum(subset(interaction_counts_melted, cell1 == clusters[i], select=interactions))  #important not to count them twice (thus only cell1)
  }
  total_interactions <- sum(interactions_counts_sums) #calculate total interactions of all cells with one another
  interactions_counts_sums <- data.frame(interactions_counts_sums, "interactions_counts_sums_frequency" = interactions_counts_sums/total_interactions) #calculate frequency of interations based on total interations
  
  expected_frequencies = interactions_counts_sums$interactions_counts_sums_frequency %o% interactions_counts_sums$interactions_counts_sums_frequency # from the vector above create its outer matrix product  = matrix with clusters and n columns and n rows, this now has
  expected_counts = expected_frequencies * total_interactions
  colnames(expected_counts) <- colnames(interactions_counts_df) #match colnames
  rownames(expected_counts) <- rownames(interactions_counts_df) #match rownames
  return(expected_counts) #return the new table of expected_counts
}

expected_counts_superclusters <- sapply(interactions_counts_superclusters, function(x) generate_expected_counts(x), simplify = FALSE) #for each interactions_counts_superclusters file generate its odds ratios


###############################################################################################################################################
#now also generate just the normalized interaction counts, without generating the odds ratios

generate_normalized_interactions <- function(interactions_counts_df){ #provide single interactions_counts/_superclusters data.frames and return their odds_ratios
  clusters <- unique(colnames(interactions_counts_df))
  x= 10^-4 #prior count to deal with 0 -> problem with logs, chosen to match Maths software data
  
  #reformat data into three columns using the melt function
  interaction_counts_melted <- setNames(melt(interactions_counts_df+x), c("cell1", "cell2", "interactions")) 
  
  #for each cluster count the total interactions for that cluster with all other clusters
  interactions_counts_sums <- c() #make emtpy vector to store the data
  for (i in 1:length(clusters)){
    interactions_counts_sums[clusters[i]] <- sum(subset(interaction_counts_melted, cell1 == clusters[i], select=interactions))  #important not to count them twice (thus only cell1)
  }
  total_interactions <- sum(interactions_counts_sums) #calculate total interactions of all cells with one another
  interactions_counts_sums <- data.frame(interactions_counts_sums, "interactions_counts_sums_frequency" = interactions_counts_sums/total_interactions) #calculate frequency of interations based on total interations
  
  expected_frequencies = interactions_counts_sums$interactions_counts_sums_frequency %o% interactions_counts_sums$interactions_counts_sums_frequency # from the vector above create its outer matrix product  = matrix with clusters and n columns and n rows, this now has
  observed_frequencies = (interactions_counts_df+x)/total_interactions  #calculate observed frequencies
  return(observed_frequencies) #return the new table of log_odds_ratios as the output of this function
}


norm_interactions_superclusters <- sapply(interactions_counts_superclusters, function(x) generate_normalized_interactions(x), simplify = FALSE) #for each interactions_counts_superclusters file generate its odds ratios




#create empty dataframe and then melt all the three comparisons of superclusters together 
superclusters_listed <- data.frame(matrix(0, nrow=0, ncol=4))

#same with conditions spread out as columns and not melted into one dataframe with other variables
a <- setNames( melt(log_odds_ratios_superclusters[[1]]) , c("cell1", "cell2", "log_odds_ratio"))
superclusters_logodds <- data.frame("comparison" = paste0(a$cell1,"_vs_", a$cell2))
superclusters_interactions <- data.frame("comparison" = paste0(a$cell1,"_vs_", a$cell2))
superclusters_expected_counts <- data.frame("comparison" = paste0(a$cell1,"_vs_", a$cell2))


for (i in 1:length(log_odds_ratios_superclusters)) {# for each condition
  a <- setNames( melt(log_odds_ratios_superclusters[[i]]) , c("cell1", "cell2", "log_odds_ratio"))
  b <- setNames( melt(norm_interactions_superclusters[[i]]) , c("cell1", "cell2", "norm_interactions_counts"))
  c <- setNames(melt(interactions_counts_superclusters[[i]]), c("cell1", "cell2", "interactions_counts")  )
  d <- setNames(melt(expected_counts_superclusters[[i]]), c("cell1", "cell2", "expected_counts")  )
  
  #this generates nxn tables of the described variables, here we want all conditions
  superclusters_logodds[,i+1] <- a$log_odds_ratio 
  superclusters_interactions[,i+1] <- c$interactions_counts
  superclusters_expected_counts[,i+1] <- d$expected_counts
  
  #now for the list with only need each condition once e.g. N_vs_CD8 is the same as CD8_vs_N now -> delete those duplicated rows
  for (j in 1:nrow(a)){ #sort the cells by place of the string in the alphabet
    k <- as.character(a$cell1[j])
    l <- as.character(a$cell2[j])
    if (k <= l) {
      cell1 <- k
      cell2 <- l
    }    else
    {
      cell2 <- k
      cell1 <- l
    }
    a$comparison[j] <- paste0(cell1,"_vs_", cell2)
  }
  keep <- !duplicated(a$comparison)
  a <- a[keep,] #now the conditions with the different directions are deleted
  b <- b[keep,]
  c <- c[keep,]
  d <- d[keep,]
  condition <- condition_mapping[grep(names(log_odds_ratios_superclusters)[i], condition_mapping$MouseNr)[1], "Treatment" ]
  #this generates a list of all comparisons for all conditions, this is important for the steps below. 
  superclusters_listed <- rbind(superclusters_listed,  data.frame( a$comparison, a$log_odds_ratio,b$norm_interactions_counts, c$interactions_counts, d$expected_counts, condition,names(log_odds_ratios_superclusters)[i]  ))
  
}
colnames(superclusters_listed) <- c(  "comparison", "log_odds_ratio","norm_interactions","interactions_counts", "expected_counts", "condition", "mouse")
colnames(superclusters_logodds) <- c(  "comparison",names(log_odds_ratios_superclusters) )
colnames(superclusters_interactions) <- c(  "comparison", names(log_odds_ratios_superclusters))
colnames(superclusters_expected_counts) <- c(  "comparison", names(log_odds_ratios_superclusters))

length(unique(superclusters_listed$comparison)) #105
head(superclusters_listed)
tail(superclusters_listed)


write.table.tab(superclusters_listed, file = "../tables/superclusters_listed.txt")




#### #### #### #### #### #### #### #### #### 

## Modelling and Testing for superclusters

hist(log10(rowMeans(superclusters_interactions[,2:ncol(superclusters_interactions)])))

#filtering of the data
counts_cutoff <- 3 #at last 3 counts
number_of_samples_cutoff <- 3 #in at least 3 samples
means_cutoff <- 3
keep <- (rowSums(superclusters_interactions[,2:ncol(superclusters_interactions)] >= counts_cutoff) >= number_of_samples_cutoff) & ( rowMeans(superclusters_interactions[,2:ncol(superclusters_interactions)]) >= 3)
summary(keep) # False 18 True 151 excluded clusters with no interactions

superclusters_interactions[!keep,] #display which are not kept, none at the moment
superclusters_interactions_sel <- superclusters_interactions[keep,]

superclusters_logodds_sel <- superclusters_logodds[keep,]
superclusters_expected_counts_sel <- superclusters_expected_counts[keep,]
superclusters_listed_sel <- superclusters_listed[superclusters_listed$comparison %in% superclusters_interactions_sel$comparison,]
superclusters_listed_sel$condition <- factor(superclusters_listed_sel$condition, levels = c("UNT", "Ad-IL12", "Ad-IL12+Ad-CCL5")) #reorder the condition names
length(unique(superclusters_listed$comparison))# 105 before selection
length(unique(superclusters_listed_sel$comparison)) #105 afer selection

####  Setting up the generalized linear models
head(superclusters_listed_sel)
ggplot(superclusters_listed_sel[superclusters_listed_sel$interactions_counts < 1000,], aes(x=interactions_counts)) + geom_histogram()

superclusters_listed_sel %>% group_by(condition) %>% summarize(sum_counts = sum(interactions_counts))


#write some files with values groupe by condition and comparison                                                                
                       
write.table.tab(superclusters_listed_sel %>% group_by(comparison, condition) %>% summarize(sum_counts = sum(interactions_counts)) %>% spread(condition, sum_counts),
            file= "../tables/interactions_counts_byConditionComparison.txt")
write.table.tab(superclusters_listed_sel %>% group_by(comparison, condition) %>% summarize(sum_expected = sum(expected_counts)) %>% spread(condition, sum_expected),
                file= "../tables/expected_counts_byConditionComparison.txt")
write.table.tab(superclusters_listed_sel %>% group_by(comparison, condition) %>% summarize(mean_log_odds_ratio = mean(log_odds_ratio)) %>% spread(condition, mean_log_odds_ratio),
                file= "../tables/mean_log_odds_ratio_byConditionComparison.txt")
write.table.tab(superclusters_listed_sel %>% group_by(comparison) %>% summarize(sum_counts = sum(interactions_counts)),
                file= "../tables/interactions_counts_byComparison.txt")
write.table.tab(superclusters_listed_sel %>% group_by(condition) %>% summarize(sum_counts = sum(interactions_counts)),
                file= "../tables/interactions_counts_byCondition.txt")

pdf(file ="../plots/NegBinomial_Model_QualityPlots.pdf", width = 8)

#### Perform a generalized negative binomial model approach on the total data of superclusters. 
listed_interactions <- superclusters_listed_sel
100*sum(listed_interactions$interactions_counts == 0)/nrow(listed_interactions) #12.15% of the data are 0

#compare three different models: negative binomial, poisson and zero inflated. Use interaction terms to estimte effects of treatments compared to untreated sample for each cell-cell comparison
listed_interactions$condition <- factor(listed_interactions$condition, levels = c("UNT", "Ad-IL12", "Ad-IL12+Ad-CCL5"))
name <- "supercluster"
model.nb <- glm.nb(interactions_counts ~ condition  * comparison + offset(log(listed_interactions$expected_counts)),  data = listed_interactions )
model.poisson <- glm(interactions_counts ~ condition  * comparison + offset(log(listed_interactions$expected_counts)), family = "poisson", data = listed_interactions )

#compare models
summary(model.poisson)
summary(model.nb)
#Neg. Binomial has by far the smallest residual deviance and the highest likelyhood

anova(model.nb) #all terms significant inlcuding interaction
pchisq(2 * (logLik(model.nb) - logLik(model.poisson)), df = 1, lower.tail = FALSE) #neg binomial much better than poisson, p value near zero

## Check for over/underdispersion in the model
E2 <- resid(model.nb, type = "pearson")
N  <- nrow(listed_interactions)
p  <- length(coef(model.nb))   
sum(E2^2) / (N - p) #negative binomial is much better than poisson! Don't use poisson because we have overdispersion. zero inflation makes it worse, so don't use

#plot predicted vs real interaction counts
predicted <- predict(model.nb, type  ="response")
plot(predicted, listed_interactions$interactions_counts)
together <- data.frame(predicted, listed_interactions$interactions_counts, listed_interactions$condition)
ggplot(data= together, aes(x = log(predicted), y = log(listed_interactions.interactions_counts), col = listed_interactions.condition)) +
  geom_point() + ggtitle("log transformed predicted vs real interaction counts") +
  geom_abline(intercept = 0 , slope = 1, col = "black" )
ggplot(data= together, aes(x = predicted, y = listed_interactions.interactions_counts, col = listed_interactions.condition)) +
  geom_point() + ggtitle(" predicted vs real interaction counts") +
  geom_abline(intercept = 0 , slope = 1, col = "black" )

#estimate significance of deviation from model. Sill significant, this probably due to zero values and some outliers
with(model.nb, cbind(res.deviance = deviance, df = df.residual,
                     p = pchisq(deviance, df.residual, lower.tail=FALSE)))

mean.var.plot = function(model.poisson,model.nb){
  xb = predict(model.nb)
  g = cut(xb, breaks=unique(quantile(xb,seq(0,1,0.1))))
  m = tapply(model.poisson$y, g, mean)
  v = tapply(model.poisson$y, g, var)
  pr <- residuals(model.poisson,"pearson")
  phi <- sum(pr^2)/df.residual(model.poisson)
  x = seq(min(m),max(m),length.out = 500)
  line.data = data.frame(x=rep(x,2),y=c(x*phi,x*(1+x/model.nb$theta)),
                         model=c(rep("Q. Poisson",length(x)),rep("Neg. Binomial",length(x))))
  library(ggplot2)
  ggplot() + geom_point(aes(x=log(m),y=log(v))) + 
    geom_line(aes(x=log(x),y=log(y),linetype=model),data=line.data) + 
    theme_bw() + theme(panel.background = element_rect(rgb(.95,.95,.95))) +
    ylab("variance") + xlab("mean") +
    scale_linetype_manual(values = c("solid","dashed")) +
    ggtitle("Mean-Variance Relationship") 
}
mean.var.plot(model.poisson,model.nb) #shows that Neg. Binomial Model fits better to the data than Poission (QuasiPoisson)

dev.off()

pdf(file ="../plots/NegBinomial_Model_Supercluster_interactions_plots.pdf", width = 8, height = 6)
#now extract the results from the neg. binomial model
results_readable <-  coef(summary(model.nb))
results_readable <- data.frame(results_readable)
results_readable$term <- rownames(results_readable)
colnames(results_readable) <- c("estimate", "sd", "zvalue", "pvalue", "term")
tail(results_readable)
results_interactors <- results_readable %>% filter(grepl(":", term))
results_interactors$condition <- gsub("condition", "",  gsub(":.*$" , "",  results_interactors$term)) #rename and extract conditions and comparisons
results_interactors$condition <- gsub("-", "", results_interactors$condition)
results_interactors$comparison <- gsub("^.*:comparison", "", results_interactors$term)
results_interactors <- results_interactors %>% arrange(pvalue) #sort by pvalue
results_interactors$OR <- exp(results_interactors$estimate)
head(results_interactors)
head(results_interactors[order(abs(results_interactors$OR), decreasing = TRUE),])
biggest_OR <- head(results_interactors[order(abs(results_interactors$OR), decreasing = TRUE),])$comparison

superclusters_listed_sel %>% filter(comparison %in% biggest_OR ) %>% group_by(condition, comparison) %>% summarise(log_odds_ratio_mean= mean(log_odds_ratio))

if (!dir.exists("../statistics")){dir.create("../statistics/")}
write.table( results_interactors, file = "../statistics/glm_neg_binomial_interactions_supercluster_B16.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#volcano
ggplot(results_interactors, aes(x=log10(abs(estimate)), y=-log10(pvalue))) + geom_point() + xlab("Absolute Log Odds Ratio")
ggplot(results_interactors, aes(x=estimate, y=-log10(pvalue))) + geom_point() + xlab("Odds Ratio")

#compare the effects of AdV5-empty and AdV5-IL12
pvalues_spread <- results_interactors %>% dplyr::select(pvalue, condition, comparison) %>% spread(condition, pvalue)
colnames(pvalues_spread) <- c("comparison", "AdIL12_p", "AdIL12_AdCCL5_p")
head(pvalues_spread)

p <- ggplot(pvalues_spread, aes(x=-log10(AdIL12_AdCCL5_p), y= -log10(AdIL12_p), label = comparison, color =  -log10(AdIL12_AdCCL5_p) )) + geom_point() + geom_abline(intercept = 0, slope = 1)  +
  geom_label_repel(data = subset(pvalues_spread, AdIL12_AdCCL5_p < 0.05 | AdIL12_p < 0.005),  size = 4.5,  force = 3, box.padding = 0.4) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 3) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("-log10 p-value AdV5-IL12+AdV5-CCL5") +
  ylab("-log10 p-value AdV5-IL12")+
  ggtitle("Change in Interactions - AdV5-IL12+Ad-CCL5 or AdV5-IL12 both vs UNT")+
  theme_bw() + 
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 10),legend.title =  element_text(size = 10) )
p


### Odds Ratios

estimate_spread <- results_interactors %>% dplyr::select(OR, condition, comparison) %>% spread(condition, OR)
colnames(estimate_spread) <- c("comparison", "AdIL12_OR", "AdIL12_AdCCL5_OR")
head(estimate_spread)

ggplot(estimate_spread, aes(x=log10(AdIL12_AdCCL5_OR), y= log10(AdIL12_OR), label = comparison)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 3) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  geom_label_repel() + 
  theme_bw() + ggtitle("Change in Odds Ratios of Interactions")

results_interactors_spread <- left_join(pvalues_spread, estimate_spread, by= "comparison")  
results_interactors_spread


subdata <- results_interactors_spread 
p <- ggplot(subdata, aes( x= log10(AdIL12_AdCCL5_OR),y=log10(AdIL12_OR), label = comparison,size = -log10(AdIL12_AdCCL5_p), color =  -log10(AdIL12_AdCCL5_p) ))  +
  ggtitle(paste(name, " Interaction Analysis -- All")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("Log Odds Ratio AdV5-IL12-AdV5-CCL5 vs untreated") +
  ylab("Log Odds Ratio AdV5-IL12 vs untreated") +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_AdCCL5_p < 0.05 | AdIL12_AdCCL5_OR < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs -log10 p-value\nAdV5-IL12-CCL5 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 10),legend.title =  element_text(size = 10) )
p

subdata <- results_interactors_spread %>% filter(grepl("NK", comparison))
p <- ggplot(subdata, aes( x= AdIL12_AdCCL5_OR,y=AdIL12_OR, label = comparison,size = AdIL12_AdCCL5_OR, color =  -log10(AdIL12_AdCCL5_p) ))  +
  ggtitle(paste(name, " Interaction Analysis -- NK")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 1, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 1, size = 0.3, linetype = 3) +
  xlab("Odds Ratio AdV5-IL12-AdV5-CCL5 vs untreated") +
  ylab("Odds Ratio AdV5-IL12 vs untreated") +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_AdCCL5_p < 0.05 | AdIL12_AdCCL5_OR < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="odds ratio\nAdV5-IL12+AdV5-CCL5 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p

subdata <- results_interactors_spread %>% filter(grepl("CD8", comparison))
p <- ggplot(subdata, aes(x= log10(AdIL12_AdCCL5_OR),y=log10(AdIL12_OR), label = comparison,size = abs(log10(AdIL12_AdCCL5_OR)), color =  -log10(AdIL12_AdCCL5_p) ))  +
  ggtitle(paste(name, " Interaction Analysis -- CD8")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("Log Odds Ratio AdV5-IL12-AdV5-CCL5 vs untreated") +
  ylab("Log Odds Ratio AdV5-IL12 vs untreated") +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_AdCCL5_p < 0.05 |  AdIL12_p < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log10 fold change\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p

subdata <- results_interactors_spread %>% filter(grepl("CD8|NK", comparison))
p <- ggplot(subdata, aes(x= log10(AdIL12_AdCCL5_OR),y=log10(AdIL12_OR), label = comparison,size = abs(log10(AdIL12_AdCCL5_OR)), color =  -log10(AdIL12_AdCCL5_p) ))  +
  ggtitle(paste(name, " Interaction Analysis -- CD8 NK")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("Log Odds Ratio AdV5-IL12-AdV5-CCL5 vs untreated") +
  ylab("Log Odds Ratio AdV5-IL12 vs untreated") +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_AdCCL5_p < 0.05 |  AdIL12_p < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log10 fold change\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p

subdata <- results_interactors_spread %>% filter(grepl("DC", comparison))
p <- ggplot(subdata, aes(x= log10(AdIL12_AdCCL5_OR),y=log10(AdIL12_OR), label = comparison,size = abs(log10(AdIL12_AdCCL5_OR)), color =  -log10(AdIL12_AdCCL5_p) ))  +
  ggtitle(paste(name, " Interaction Analysis -- DC")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("Log Odds Ratio AdV5-IL12-AdV5-CCL5 vs untreated") +
  ylab("Log Odds Ratio AdV5-IL12 vs untreated") +
  xlim(-1.5,1.5) +
  ylim(-1.5,1.5) +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_AdCCL5_p < 0.05 |  AdIL12_p < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log10 fold change\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p


dev.off()


###########################################################################################################3
#same now for the contrast IL12-CCL5 vs IL12 only
#

listed_interactions$condition <- factor(listed_interactions$condition, levels = c("Ad-IL12", "UNT", "Ad-IL12+Ad-CCL5"))
name <- "supercluster"
model.nb <- glm.nb(interactions_counts ~ condition  * comparison + offset(log(listed_interactions$expected_counts)),  data = listed_interactions )

summary(model.nb)

pdf(file ="../plots/NegBinomial_Model_Supercluster_interactions_plots_othercontrast.pdf", width = 8, height = 6)
#now extract the results from the neg. binomial model
results_readable <-  coef(summary(model.nb))
results_readable <- data.frame(results_readable)
results_readable$term <- rownames(results_readable)
colnames(results_readable) <- c("estimate", "sd", "zvalue", "pvalue", "term")
tail(results_readable)
results_interactors <- results_readable %>% filter(grepl(":", term))
results_interactors$condition <- gsub("condition", "",  gsub(":.*$" , "",  results_interactors$term)) #rename and extract conditions and comparisons
results_interactors$condition <- gsub("-", "", results_interactors$condition)
results_interactors$comparison <- gsub("^.*:comparison", "", results_interactors$term)
results_interactors <- results_interactors %>% arrange(pvalue) #sort by pvalue
results_interactors$OR <- exp(results_interactors$estimate)
head(results_interactors)
head(results_interactors[order(abs(results_interactors$OR), decreasing = TRUE),])
biggest_OR <- head(results_interactors[order(abs(results_interactors$OR), decreasing = TRUE),])$comparison

superclusters_listed_sel %>% filter(comparison %in% biggest_OR ) %>% group_by(condition, comparison) %>% summarise(log_odds_ratio_mean= mean(log_odds_ratio))

if (!dir.exists("../statistics")){dir.create("../statistics/")}
write.table( results_interactors, file = "../statistics/glm_neg_binomial_interactions_supercluster_othercontrast_B16.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#compare the effects of AdV5-empty and AdV5-IL12
pvalues_spread <- results_interactors %>% dplyr::select(pvalue, condition, comparison) %>% spread(condition, pvalue)
colnames(pvalues_spread) <- c("comparison", "AdIL12_AdCCL5_p", "UNT_p")
head(pvalues_spread)

p <- ggplot(pvalues_spread, aes(x=-log10(AdIL12_AdCCL5_p), y= -log10(UNT_p), label = comparison, color =  -log10(AdIL12_AdCCL5_p) )) + geom_point() + geom_abline(intercept = 0, slope = 1)  +
  geom_label_repel(data = subset(pvalues_spread, AdIL12_AdCCL5_p < 0.05 | UNT_p < 0.005),  size = 4.5,  force = 3, box.padding = 0.4) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 3) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("-log10 p-value AdV5-IL12+AdV5-CCL5 vs AdV5-IL12") +
  ylab("-log10 p-value UNT vs AdV5-IL12")+
  ggtitle("Change in Interactions - AdV5-IL12+Ad-CCL5 or UNT vs AdV5-IL12")+
  theme_bw() + 
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12-AdV5-CCL5 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 10),legend.title =  element_text(size = 10) )
p


### Odds Ratios

estimate_spread <- results_interactors %>% dplyr::select(OR, condition, comparison) %>% spread(condition, OR)
colnames(pvalues_spread) <- c("comparison", "AdIL12_AdCCL5_OR", "UNT_OR")
head(estimate_spread)

dev.off()

