# This code describes the analysis of the EMT6-HER2 CODEX multiple imaging cell interaction counts
# performed by Marcel P. Trefny in 2020/21


# ### Install necessary packages
# #only need to do this once, thus commented out
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# install.packages("circlize")
# install.packages("ggplot2")
# install.packages("devtools")
# install.packages("multcomp")
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

# Load Packages
library(viridis)
library(devtools)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(circlize)
library(ComplexHeatmap)
# library(flowCore)
# library(ggcyto)
library(tidyr)
library(edgeR)
library(multcomp)
library(gridExtra)
library(RColorBrewer)
library(forcats)

set.seed(123) #define randomness generation 


myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))

#make custom functions
#custom function to write tab deliminated text with the correct settings
write.table.tab <- function(x, file = ""){
  write.table(x,file,sep = "\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
}

### Import Data files
path <- "/Volumes/CIMM$/CIMM/Codex/Mouse-CODEX/Nicole/EMT6_Codex_InteractionAnalysis_GithubDeposition/Interaction_counts_renamed/" #change to path of r script


setwd(path)
interaction_count_files <- list.files(pattern = "interactioncount_mtx.csv")

### extract interaction counts in list format
interactions_counts <- sapply(interaction_count_files, function(x) data.matrix(read.csv(x, row.names = 1, check.names= FALSE)), simplify = FALSE, USE.NAMES = TRUE)

# are column names and row names identical? should be all true if correct
sapply(interactions_counts, function(x) {summary(colnames(x) == rownames(x))}, simplify = TRUE)

#now clean up data to have all clusters in all samples
clusters = unique(colnames(interactions_counts[[1]])) #names of clusters

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


#define which clusters to exclude
#beware they need to be excluded again for the the superclusters
exclusion_clusters <- c("L")
# exclusion_clusters <- c("")
#exclude these clusters
interactions_counts <- sapply(interactions_counts, function(x) {a <- colnames(x) == exclusion_clusters; b <- rownames(x) == exclusion_clusters; x[!b,!a] } , simplify = FALSE) #68 true all



#save files now with all columns and rows even if empty
if (!dir.exists("../interaction_counts_filled")){dir.create("../interaction_counts_filled/")}
sapply(names(interactions_counts),function(x) write.table(interactions_counts[[x]], file = paste0("../interaction_counts_filled/", gsub(".csv", "", x) , "_filled.txt"),sep = "\t", quote = FALSE))

########################################################################################################
supercluster_mapping <- read.table("../input/supercluster_annotation.txt", header = TRUE , sep = "\t")
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


#save the interaction_counts_supercluster files
if (!dir.exists("../interaction_counts_superclusters/")){dir.create("../interaction_counts_superclusters/")}
sapply(names(interactions_counts_superclusters),function(x) write.table(interactions_counts_superclusters[[x]], file = paste0("../interaction_counts_superclusters/", gsub(".csv", "", x) , "_by_supercluster.txt"), sep = "\t",quote = FALSE))

#sapply(interactions_counts_superclusters,function(x) Heatmap((log10(x+1))))
 

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

log_odds_ratios <- sapply(interactions_counts, function(x) generate_log_odds_ratios(x), simplify = FALSE)
log_odds_ratios_superclusters <- sapply(interactions_counts_superclusters, function(x) generate_log_odds_ratios(x), simplify = FALSE) #for each interactions_counts_superclusters file generate its odds ratios

#save files to different folders
if (!dir.exists("../log_odds_ratios_clusters/")){dir.create("../log_odds_ratios_clusters/")}
sapply(names(log_odds_ratios),function(x) write.table(log_odds_ratios[[x]], file = paste0("../log_odds_ratios_clusters/", gsub(".csv", "", x) , "log_odds_ratios.txt"),sep = "\t", quote = FALSE))
if (!dir.exists("../log_odds_ratios_superclusters/")){dir.create("../log_odds_ratios_superclusters/")}
sapply(names(log_odds_ratios_superclusters),function(x) write.table(log_odds_ratios_superclusters[[x]], file = paste0("../log_odds_ratios_superclusters/", gsub(".csv", "", x) , "log_odds_ratios.txt"), sep = "\t", quote = FALSE))

#Generate Heatmaps for individual odds_ratios
if (!dir.exists("../plots/")){dir.create("../plots/")}
pdf("../plots/Heatmaps_log_odds_individually.pdf")
sapply(names(log_odds_ratios),function(x) {Heatmap((log_odds_ratios[[x]]), cluster_rows = TRUE, cluster_columns = TRUE , column_title = gsub("_interactioncount_mtx.csv","", x), name = "log10 odds ratio" )})
sapply(names(log_odds_ratios_superclusters),function(x) {Heatmap((log_odds_ratios_superclusters[[x]]), cluster_rows = TRUE, cluster_columns = TRUE , column_title = gsub("_interactioncount_mtx.csv","", x), name = "log10 odds ratio" )})
dev.off()



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

expected_counts <- sapply(interactions_counts, function(x) generate_expected_counts(x), simplify = FALSE)
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

norm_interactions<- sapply(interactions_counts, function(x) generate_normalized_interactions(x), simplify = FALSE)
norm_interactions_superclusters <- sapply(interactions_counts_superclusters, function(x) generate_normalized_interactions(x), simplify = FALSE) #for each interactions_counts_superclusters file generate its odds ratios


#####################################################################################################################


#load the mapping of which condition is which sample
condition_mapping <- read.table("../input/condition_mapping.txt", sep = "\t", header = TRUE)
condition_mapping

#create empty dataframe and then melt all the three comparisons of superclusters together 
superclusters_listed <- data.frame(matrix(0, nrow=0, ncol=3))

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
   condition <- condition_mapping$condition[grep(gsub("_interactioncount_mtx.csv", "", names(log_odds_ratios_superclusters)[i]), condition_mapping$filename) ] 
   #this generates a list of all comparisons for all conditions, this is important for the steps below. 
   superclusters_listed <- rbind(superclusters_listed,  data.frame( a$comparison, a$log_odds_ratio,b$norm_interactions_counts, c$interactions_counts, d$expected_counts, condition ))

   }
colnames(superclusters_listed) <- c(  "comparison", "log_odds_ratio","norm_interactions","interactions_counts", "expected_counts", "condition")
colnames(superclusters_logodds) <- c(  "comparison",names(log_odds_ratios_superclusters) )
colnames(superclusters_interactions) <- c(  "comparison", names(log_odds_ratios_superclusters))
colnames(superclusters_expected_counts) <- c(  "comparison", names(log_odds_ratios_superclusters))

length(unique(superclusters_listed$comparison)) #91


#create empty dataframe and then melt all the three comparisons of large number of clusters together 
clusters_listed <- data.frame(matrix(0, nrow=0, ncol=3))
#same with conditions spread out as columns and not melted into one dataframe with other variables
a <- setNames( melt(log_odds_ratios[[1]]) , c("cell1", "cell2", "log_odds_ratio"))
clusters_logodds <- data.frame("comparison" = paste0(a$cell1,"_vs_", a$cell2))
clusters_interactions <- data.frame("comparison" = paste0(a$cell1,"_vs_", a$cell2))
clusters_expected_counts <- data.frame("comparison" = paste0(a$cell1,"_vs_", a$cell2))

nrow(clusters_interactions)

for (i in 1:length(log_odds_ratios)) { #for each condition
  a <- setNames( melt(log_odds_ratios[[i]]) , c("cell1", "cell2", "log_odds_ratio"))
  b <- setNames( melt(norm_interactions[[i]]) , c("cell1", "cell2", "norm_interactions_counts"))
  c <- setNames(melt(interactions_counts[[i]]), c("cell1", "cell2", "interactions_counts")  )
  d <- setNames(melt(expected_counts[[i]]), c("cell1", "cell2", "expected_counts")  )
  
  clusters_logodds[,i+1] <- a$log_odds_ratio 
  clusters_interactions[,i+1] <- c$interactions_counts
  clusters_expected_counts[,i+1] <- d$expected_count
  
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
  a <- a[keep,] #now the conditions with the different directions e.g. N_vs_CD8 is the same as CD8_vs_N now -> delete those duplicated rows
  b <- b[keep,]
  c <- c[keep,]
  d <- d[keep,] 
  
  condition <- condition_mapping$condition[grep(gsub("_interactioncount_mtx.csv", "", names(log_odds_ratios)[i]), condition_mapping$filename) ] 
  clusters_listed <- rbind(clusters_listed, data.frame(paste0(a$cell1,"_vs_", a$cell2), a$log_odds_ratio,b$norm_interactions_counts, c$interactions_counts,d$expected_counts,condition) )
   
}

colnames(clusters_listed) <- c( "comparison", "log_odds_ratio","norm_interactions","interactions_counts","expected_counts", "condition")
colnames(clusters_logodds) <- c( "comparison", gsub("_interactioncount_mtx.csv", "", names(log_odds_ratios)))
colnames(clusters_interactions) <-c( "comparison", gsub("_interactioncount_mtx.csv", "", names(log_odds_ratios)))
colnames(clusters_expected_counts) <- c( "comparison", gsub("_interactioncount_mtx.csv", "", names(log_odds_ratios)))

length(unique(clusters_listed$comparison)) #2682
head(clusters_listed)

# clusters_listed %>% group_by(comparison) %>% tally() 

if (!dir.exists("../listed_data/")){dir.create("../listed_data/")}
write.table(superclusters_listed, file = paste0("../listed_data/", "superclusters_listed.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(clusters_listed, file = paste0("../listed_data/", "clusters_listed.txt"),  sep = "\t", quote = FALSE, row.names = FALSE)




#### #### #### #### #### #### #### #### #### 

## Modelling and Testing for only for superclusters


#filtering of the data
mean_cutoff <- -1 # we want all clusters for now
keep <- (rowMeans(superclusters_interactions[,2:ncol(superclusters_interactions)]) > mean_cutoff)
summary(keep) # False 0 True 169

superclusters_interactions[!keep,] #display which are not kept, none at the moment
superclusters_interactions_sel <- superclusters_interactions[keep,]

superclusters_logodds_sel <- superclusters_logodds[keep,]
superclusters_expected_counts_sel <- superclusters_expected_counts[keep,]
superclusters_listed_sel <- superclusters_listed[superclusters_listed$comparison %in% superclusters_interactions_sel$comparison,]
superclusters_listed_sel$condition <- factor(superclusters_listed_sel$condition, levels = c("Ad-IL12", "Ad-empty", "untreated")) #reorder the condition names
length(unique(superclusters_listed$comparison))# 91 before selection
length(unique(superclusters_listed_sel$comparison)) #90 afer selection

####  Setting up the generalized linear models
head(superclusters_listed_sel)
ggplot(superclusters_listed_sel[superclusters_listed_sel$interactions_counts < 1000,], aes(x=interactions_counts)) + geom_histogram()




pdf(file ="../plots/NegBinomial_Model_Supercluster_Plots.pdf", width = 8)

#### Perform a generalized negative binomial model approach on the total data of superclusters. 
listed_interactions <- superclusters_listed_sel
100*sum(listed_interactions$interactions_counts == 0)/nrow(listed_interactions) #1.34% of the data are 0

#compare three different models: negative binomial, poisson and zero inflated. Use interaction terms to estimte effects of treatments compared to untreated sample for each cell-cell comparison
listed_interactions$condition <- factor(listed_interactions$condition, levels=c("untreated", "Ad-empty", "Ad-IL12") )
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
sum(E2^2) / (N - p) #negative binomial is much better than poisson! Don't use poisson because we have overdispersion. 

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

write.table( results_interactors, file = "../statistics/glm_neg_binomial_interactions_supercluster.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#volcano
ggplot(results_interactors, aes(x=OR, y=-log10(pvalue))) + geom_point()

#compare the effects of AdV5-empty and AdV5-IL12
pvalues_spread <- results_interactors %>% dplyr::select(pvalue, condition, comparison) %>% spread(condition, pvalue)
colnames(pvalues_spread) <- c("comparison", "Adempty_p", "AdIL12_p")

ggplot(pvalues_spread, aes(x=-log10(Adempty_p), y= -log10(AdIL12_p), label = comparison)) + geom_point() + geom_abline(intercept = 0, slope = 1)  +
 geom_label_repel(data = subset(pvalues_spread, Adempty_p < 0.005 | AdIL12_p < 0.005)) +
  geom_hline(yintercept = -log10(0.005)) + geom_vline(xintercept = -log10(0.005)) +
  ggtitle("Change in Interactions - Significance - AdV5-empty vs AdV5-IL12")

estimate_spread <- results_interactors %>% dplyr::select(OR, condition, comparison) %>% spread(condition, OR)
colnames(estimate_spread) <- c("comparison", "Adempty_OR", "AdIL12_OR")
zvalue_spread <- results_interactors %>% dplyr::select(zvalue, condition, comparison) %>% spread(condition, zvalue)
colnames(zvalue_spread) <- c("comparison", "Adempty_z", "AdIL12_z")
head(estimate_spread)

ggplot(estimate_spread, aes(x=Adempty_OR, y= AdIL12_OR, label = comparison)) + geom_point() + geom_abline(intercept = 0, slope = 1) + 
  geom_hline(yintercept = 1) + geom_vline(xintercept = 1) +
  geom_label_repel(data = subset(estimate_spread, AdIL12_OR - Adempty_OR > 2)) + 
  theme_bw()

results_interactors_spread <- left_join(pvalues_spread, estimate_spread,  by= "comparison")  
results_interactors_spread <- left_join(results_interactors_spread, zvalue_spread,  by= "comparison")  

results_interactors_spread

subdata <- results_interactors_spread %>% filter(grepl("CD8|NK", comparison))
p <- ggplot(subdata, aes(y=Adempty_OR, x= AdIL12_OR, label = comparison,size = AdIL12_OR, color =  -log10(AdIL12_p) ))  +
  ggtitle(paste(name, " Interaction Analysis")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 1, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 1, size = 0.3, linetype = 3) +
  xlab("Odds Ratio AdV5-IL12 vs untreated") +
  ylab("Odds Ratio AdV5-control vs untreated") +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_p < 0.05 |  Adempty_p < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log2 fold change\nAdV5-IL12 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p

subdata <- results_interactors_spread %>% filter(grepl("CD8|NK", comparison))
p <- ggplot(subdata, aes(y=log(Adempty_OR), x= log(AdIL12_OR), label = comparison,size = AdIL12_OR, color =  -log10(AdIL12_p) ))  +
  ggtitle(paste(name, " Interaction Analysis")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("Log Odds Ratio AdV5-IL12 vs untreated") +
  ylab("Log Odds Ratio AdV5-control vs untreated") +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_p < 0.05 |  Adempty_p < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log2 fold change\nAdV5-IL12 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p

subdata <- results_interactors_spread %>% filter(grepl("CD8|NK", comparison))
p <- ggplot(subdata, aes(y=Adempty_z, x= AdIL12_z, label = comparison,size = AdIL12_z, color =  -log10(AdIL12_p) ))  +
  ggtitle(paste(name, " Interaction Analysis")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("z-value AdV5-IL12 vs untreated") +
  ylab("z-value AdV5-control vs untreated") +
  geom_point()+
  geom_label_repel(data = subset(subdata, AdIL12_p < 0.05 |  Adempty_p < 0.05  ), size = 4.5,  force = 3, box.padding = 0.4) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log2 fold change\nAdV5-IL12 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p

subdata <- results_interactors_spread 
p <- ggplot(subdata, aes(y=log(Adempty_OR), x= log(AdIL12_OR), label = comparison,size = AdIL12_OR, color =  -log10(AdIL12_p) ))  +
  ggtitle(paste(name, " Interaction Analysis - All Cell Types")) +
  geom_abline(intercept = 0, slope = 1, size = 0.3, linetype = 2) +
  geom_hline(yintercept = 0, size = 0.3, linetype = 3) +
  geom_vline(xintercept = 0, size = 0.3, linetype = 3) +
  xlab("Log Odds Ratio AdV5-IL12 vs untreated") +
  ylab("Log Odds Ratio AdV5-control vs untreated") +
  geom_point()+
  geom_text_repel(data = subset(subdata, AdIL12_p < 0.01 |  Adempty_p < 0.05  ), size = 2.5,  force = 3) +
  scale_color_gradient( low="blue",high="red", space ="Lab" , name = "-log10 p-value\nAdV5-IL12 vs untreated") +
  #scale_color_manual(values=c("#57BFFA", "#DE2902"))
  scale_size(name="abs log2 fold change\nAdV5-IL12 vs untreated") +
  theme_bw() + 
  theme(axis.text = element_text(size= 12), axis.title = element_text(size = 15), legend.text = element_text(size = 12),legend.title =  element_text(size = 12) )
p


dev.off()




