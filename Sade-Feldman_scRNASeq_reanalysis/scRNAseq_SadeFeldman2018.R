# This is the re-analysis of the scRNAseq data by Sade-Feldman et al focusing on NK cells
# original raw data from GSE120575

# # BiocManager::install("TENxPBMCData")
# BiocManager::install("SingleCellExperiment", ask = FALSE)
# BiocManager::install("scater", ask = FALSE)
# BiocManager::install("uwot", ask = FALSE)
# BiocManager::install("scran", ask = FALSE)
# BiocManager::install("Rtsne", ask = FALSE)
# devtools::install_github("Irrationone/cellassign", ask = FALSE)
# BiocManager::install("batchelor", ask = FALSE)
# BiocManager::install("BiocSingular", ask = FALSE)
# BiocManager::install("BiocNeighbors", ask = FALSE)
# BiocManager::install("rsvd", ask = FALSE)
# BiocManager::install("igraph", ask = FALSE)
# BiocManager::install("cellassign", ask = FALSE)
# BiocManager::install("rowr", ask = FALSE)

library(batchelor)
library(BiocSingular)
library(BiocNeighbors)
library(Rtsne)
library(rsvd)
library(tidyr) 
library(dplyr)
library(igraph)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(ggplot2)
# library(rowr)
########################################################################

#Data Loading

path <- "/Volumes/CIMM$/CIMM/Projects_Human_Team/Marcel_Nicole/AdV5_IL12_Project/SadeFeldmann_Cell2018/input data/" #change to file path


setwd(path)

#this takes very long to run, thus saved data in rds file
# files_row_column <- read.table("GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt" , sep = "\t", skip = 2,  check.names= FALSE, header = FALSE, stringsAsFactors = FALSE )

# counts <- files_row_column[,-c(1)] #first column are gene names
# counts <- files_row_column[,-c(1,16292)] #last column is an artefact from the file loading
# 
# files_column_data_header <- read.table("single_cells_TPM_GEO_header.txt",  sep = "\t",  check.names= FALSE, header = FALSE, stringsAsFactors = FALSE )
# files_column_data_header <- data.frame(t(files_column_data_header))
# colnames(files_column_data_header) <- c("title", "patientID")
# head(files_column_data_header)
# dim(counts)
# dim(files_column_data_header)
# colnames(counts) <- files_column_data_header$title
# 
# counts_dg <- as(as.matrix(counts), "dgCMatrix")
# saveRDS(counts_dg, file = "counts_dg.rds")


counts_dg <- readRDS("counts_dg.rds")
head(counts_dg)


dim(counts_dg) #55737 16291

files_column_data <- read.table("GSE120575_patient_ID_single_cells_formatted.txt",  sep = "\t",  check.names= FALSE, header = TRUE, stringsAsFactors = FALSE )
files_column_data <- files_column_data[,1:7]
colnames(files_column_data) <- c("sample", "title", "source", "organism", "patientID", "response", "therapy")
head(files_column_data)
dim(files_column_data) #16291     7

genes <- read.table("genes.txt",  sep = "\t",  check.names= FALSE, header = FALSE, stringsAsFactors = FALSE )
head(genes)

counts_dg[55730:55737,16285:16291] #the last column is NA, thus remove
counts_dg <- counts_dg[,-c(ncol(counts_dg))] 
files_column_data <- files_column_data[-c(nrow(files_column_data)),] #also remove this cell's information


sce <- SingleCellExperiment(assays=SimpleList(counts=counts_dg), rowData = data.frame("GeneSymbol" = genes$V1),
                            colData=files_column_data)
rownames(sce) <- make.names(rowData(sce)$GeneSymbol) #need to make_names because of genes that start with a number
sce

sce$timepoint <- gsub("_.*$","", colData(sce)$patientID)
sce$patient <- gsub("^.*_P","P", colData(sce)$patientID)
sce$timepoint <- factor(sce$timepoint, levels = c("Pre", "Post"))

saveRDS(sce,file = "SadeFeldman_single_cell_experiment_allCells.rds")


dim(assay(sce, "counts"))
head(rowData(sce))
head(colData(sce))
nrow(counts(sce))
ncol(counts(sce))



mean(counts(sce) == 0)
class(counts(sce))
counts(sce)[1:10, 1:10]


## Total number of detected transcripts
nrow(counts(sce)) #55737
sum(rowSums(counts(sce)) > 0)

########################################################################################################
############# Filtering of genes and cells
########################################################################################################
protein_coding_genes <- read.table("ProteinCodingGenes.txt", header = TRUE, sep = "\t")
protein_coding_genes <- protein_coding_genes$GeneSymbol

#focus on protein coding genes as in publication
summary(rowData(sce)$GeneSymbol %in% protein_coding_genes)
# Mode   FALSE    TRUE 
# logical   30146   20367 
sce <- sce[ rowData(sce)$GeneSymbol %in% protein_coding_genes,] 

#exclude all transcripts which were not detected in any cell
sce <- sce[rowSums(counts(sce)) > 0,] #50513

## Genes detected in single cells
summary(colSums(counts(sce) > 0)) #mean 1997 genes per cell

#exclude non-immune cells as in paper, but only with CD45 as cutoff
non_immune <- counts(sce)[rowData(sce)$GeneSymbol == "PTPRC", ] == 0
sce <- sce[,!non_immune] 


housekeeping_genes <- read.table("HousekeepingGenes.txt", header = TRUE, sep = "\t")
housekeeping_genes <- housekeeping_genes$GeneSymbol

housekeeping_mean_per_cell <- colMeans( counts(sce)[rowData(sce)$GeneSymbol %in% housekeeping_genes,])
summary(housekeeping_mean_per_cell > 2.5) #True for almost all cells
plot(sort(housekeeping_mean_per_cell), breaks = 30)
sce <- sce[,housekeeping_mean_per_cell > 2.5] # take only cells with more than 2.5 log2(TPM+1) as in paper

#exclude genes that were not expressed with log2TPM+1 of >4.5 in at least 10 cells
summary(rowSums(counts(sce) > 4.5)>10) 
# Mode   FALSE    TRUE 
# logical    4521   15846 
sce <- sce[rowSums(counts(sce) > 4.5)>10,]

# #for each patient, see whether each gene has at least 4.5 log2 TPM+1 in at least 10 cells
# genecounts_perPatient <- data.frame(genes = rowData(sce)$GeneSymbol)
# 
# for(i in 1:length(unique(colData(sce)$patientID))){
#   patient <- unique(colData(sce)$patientID)[i]
#   genecounts_perPatient <- cbind(genecounts_perPatient, rowSums(counts(sce)[,colData(sce)$patientID == patient] > 2.5) > 5)
# }
# rownames(genecounts_perPatient) <- genecounts_perPatient$genes
# genecounts_perPatient <- genecounts_perPatient[,-c(1)]
# 
# summary(rowSums(genecounts_perPatient)>= length(unique(colData(sce)$patientID))) 



#Dropout plot
# plot(rowSums(counts(sce)), rowMeans(counts(sce) == 0), log = "x")



####################################
#Extract QC parameters and plot

mt <- rownames(sce)[grep("^MT-", rowData(sce)$GeneSymbol)] #mitochondrial genes
mt

sce <- addPerFeatureQC(sce) #scores on all genes
sce <- addPerCellQC(sce, subsets = list(MT = mt)) #scores on mitochondrial genes

plotColData(sce, x = "sum", y="detected", colour_by="patientID") 
colnames(rowData(sce))
colnames(colData(sce))

# scater::plotHighestExprs(sce, n = 15)
# plot(sort(colData(sce)$detected), breaks = 30) #genes detected per cell. did not need to exclude because this is already done above
plot(sort(sce$percent_top_50), breaks = 30, main = "percent_top_50")
plot(sort(sce$subsets_MT_percent), breaks = 30,  main = "percent mitochondrial")

############################################################
#Filtering of low gene detection and high mitochondrial DNA = death contamination

# low_detected <- isOutlier(sce$detected, type = "lower", 
#                           log = TRUE, nmads = 4)
high_mt <- isOutlier(sce$subsets_MT_percent, type = "higher",
                     log = FALSE, nmads = 3) #cells that are dying -> 3 standard deviations from mean mitochondrial percentage away
duplets <- isOutlier(sce$sum, type = "higher",
                     log = FALSE, nmads = 4) #doublets = cells with lots of genes detected

#cells with low number of genes detected are excluded
# plot(rank(-sce$detected), sce$detected, col = low_detected + 1 )
#duplets = very high number of counts
plot(rank(-sce$sum), sce$sum, col = duplets + 1 )
# dying cells = high content of mtDNA
plot(rank(sce$subsets_MT_percent), sce$subsets_MT_percent,
     col = high_mt + 1)

#based on this cut away the unwanted cells
sce$retain <-  !high_mt & !duplets
table(sce$retain) 
## 
#FALSE  TRUE 
# 519 15390 

sce <- sce[, sce$retain] #delete cells that did not pass QC

dim(sce)
# 15846 15390



########################################################################
#### Normalization
#look at the different library sizes of each cell with some patients having much more transcripts per cell sequenced
plot(colSums(counts(sce)), breaks = 30 )

cpms <- assay(sce, "counts")
libsizes <- colSums(cpms)
size.factors <- libsizes/mean(libsizes)
logcounts(sce) <- t(t(cpms)/size.factors) + 1
assayNames(sce)
plot(colSums(logcounts(sce)), breaks = 30 )


## Fit a trend
dec.trend <- modelGeneVar(sce)
fit.trend <- metadata(dec.trend)
plot(fit.trend$mean, fit.trend$var, xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")
curve(fit.trend$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)
head(dec.trend[order(dec.trend$bio, decreasing = TRUE), ])


sce <- runPCA(sce, exprs_values = "logcounts", ncomponents = 50, 
              ntop = 4000)
sce
reducedDimNames(sce)


pdf(file = "../plots/reducedDimensions_PCA_tSNE_UMAP.pdf", width = 8)

plotReducedDim(sce, "PCA", colour_by = "detected")
plotReducedDim(sce, "PCA", colour_by = "patientID")
plotReducedDim(sce, "PCA", colour_by = "therapy")


set.seed(123)

#perform tSNE
sce <- runTSNE(sce, dimred = "PCA",  perplexity = 30)
reducedDimNames(sce)

sce <- runUMAP(sce, dimred = "PCA")
reducedDimNames(sce)

plotReducedDim(sce, "TSNE", colour_by = "subsets_MT_percent")
plotReducedDim(sce, "TSNE", colour_by = "sum")
plotReducedDim(sce, "TSNE", colour_by = "detected")
plotReducedDim(sce, "TSNE", colour_by = "response")
plotReducedDim(sce, "TSNE", colour_by = "patientID")
plotReducedDim(sce, "TSNE", colour_by = "therapy")


#recreate their cluster ids from paper supplentary file
cluster_ID_publication <- read.table("cluster_ID_Paper.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(cluster_ID_publication)
summary(cluster_ID_publication$title %in% colData(sce)$title)

sce$clusterSadeFeldman <-  left_join(data.frame(colData(sce)), cluster_ID_publication, by = "title")$clusterSadeFeldman
sce$clusterSadeFeldman <- as.factor(sce$clusterSadeFeldman)

plotReducedDim(sce, "TSNE", colour_by = "clusterSadeFeldman", text_by = "clusterSadeFeldman")

#their clusters are not so nice, thus do my own clustering:

# Create a graph of shared nearest neighbors:
g <- scran::buildSNNGraph(sce, k=10, use.dimred = 'PCA') # Build SNN graph
clust <- igraph::cluster_louvain(g)$membership # use the louvain method for community detection
table(clust)
sce$clust <- factor(clust)

#plot newly generated clusters
plotReducedDim(sce, "TSNE", colour_by = "clust", text_by = "clust")


#find markers that distinguish these clusters
markers_all <- scran::findMarkers(
  sce, groups = sce$clust, assay.type="logcounts",
  pval.type = "all"
)
scater::plotExpression(sce, features = rownames(head(markers_all[[1]])), 
                       x = "clust")
head(markers_all[[1]], n= 20)
head(rownames(markers_all[[1]]), n = 50)




plotReducedDim(sce, "TSNE", colour_by = "NCR1")
plotReducedDim(sce, "TSNE", colour_by = "NCAM1")
plotReducedDim(sce, "TSNE", colour_by = "CD3E")
plotReducedDim(sce, "TSNE", colour_by = "CD8A")
plotReducedDim(sce, "TSNE", colour_by = "PDCD1")
plotReducedDim(sce, "TSNE", colour_by = "HAVCR2")
plotReducedDim(sce, "TSNE", colour_by = "TIGIT")
plotReducedDim(sce, "TSNE", colour_by = "TCF7")
plotReducedDim(sce, "TSNE", colour_by = "LEF1")
plotReducedDim(sce, "TSNE", colour_by = "CD2")
plotReducedDim(sce, "TSNE", colour_by = "CD14")
plotReducedDim(sce, "TSNE", colour_by = "CD33")
plotReducedDim(sce, "TSNE", colour_by = "FUT4")
plotReducedDim(sce, "TSNE", colour_by = "CLEC4C")
plotReducedDim(sce, "TSNE", colour_by = "THBD")
plotReducedDim(sce, "TSNE", colour_by = "ITGAX")
plotReducedDim(sce, "TSNE", colour_by = "CD4")
plotReducedDim(sce, "TSNE", colour_by = "XBP1")
plotReducedDim(sce, "TSNE", colour_by = "CD163")
plotReducedDim(sce, "TSNE", colour_by = "CD19")
plotReducedDim(sce, "TSNE", colour_by = "CCL5")
plotReducedDim(sce, "TSNE", colour_by = "CCR5")
plotReducedDim(sce, "TSNE", colour_by = "IL12A")
plotReducedDim(sce, "TSNE", colour_by = "FCGR3A")

plotReducedDim(sce, "UMAP", colour_by = "clust", text_by = "clust")
plotReducedDim(sce, "UMAP", colour_by = "FCGR3A", text_by = "clust")

ggcells(sce, mapping=aes_string(x = "clust", y = "CD3E"), exprs_values = "logcounts") +
  geom_violin() + ylab("log normalized counts CD3E") + theme_bw()  + geom_jitter(shape=16,size = 0.7, height= 0.4 )
ggcells(sce, mapping=aes_string(x = "clust", y = "FCGR3A"), exprs_values = "logcounts") +
  geom_violin() + ylab("log normalized counts CD3E") + theme_bw()  + geom_jitter(shape=16,size = 0.7, height= 0.4 )

plotReducedDim(sce, "UMAP", colour_by = "FCGR3A", text_by = "clust")

dev.off()

#NK cluster number 8 is still contaminated with NKT cells that express CD3E. Thus define new cluster called "isNK" which is cluster 8 minus cells that express CD3E
sce$isNK <- (sce$clust == "1") & (counts(sce)[rowData(sce)$GeneSymbol == "CD3E",] == 0)
summary(sce$clust == "1")
summary(sce$isNK)
# Mode   FALSE    TRUE 
# logical   15044     346 
#see where these NK cells are

pdf(file = "../plots/selected_cell_expression_plots.pdf", width = 8)

plotReducedDim(sce, "TSNE", colour_by = "isNK", text_by = "clust")
ggcells(sce[,sce$isNK], mapping=aes_string(x= "response", y = "CCL5", col = "timepoint" ) ) + geom_boxplot()  + theme_bw()

ggcells(sce[,sce$isNK], mapping=aes_string(x= "response", y = "CCL5", fill = "timepoint" ) ) + 
  geom_violin(position=position_dodge(0.8))  +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.4) +
  theme_bw() + 
  ylab("CCL5 log2(TPM+1)")+
  ggtitle("CCL5 Expression in NK cells")

ggcells(sce[,sce$clust == 1], mapping=aes_string(x= "response", y = "CCL5", fill = "timepoint" ) ) + 
  geom_violin(position=position_dodge(0.8))  +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.4) +
  theme_bw() + 
  ylab("CCL5 log2(TPM+1)")+
  ggtitle("CCL5 Expression in NK cells")


dev.off()

##################################################################################33
signatures_files <- list.files("../signatures/")

signatures <- sapply(signatures_files, function(x) read.table(paste0("../signatures/", x) , sep = "\t", stringsAsFactors = FALSE  )$V1, USE.NAMES = TRUE)
signatures

sapply(signatures, function(x) summary(x %in% rowData(sce)$GeneSymbol))
sapply(signatures, function(x) x[!x %in% rowData(sce)$GeneSymbol]) #which are not in the dataset

signatures <- sapply(signatures, function(x) x[x %in% rowData(sce)$GeneSymbol]) #only take those that are found
names(signatures) <- gsub(".txt","", names(signatures)) #give them good names
signatures

# # only done once because too big to do each time
# # for each cell, calculate signature score for all signatures
signatures_byCell <- data.frame(matrix(ncol=length(signatures), nrow= nrow(colData(sce))))

location <- which(colData(sce)$clust == "1")

for (j in 1:length(signatures)){ #for each signature
  for(i in 1:length(location)) { # for each cell
    df <- data.frame(cpms = logcounts(sce)[,location[i]], genes = rowData(sce)$GeneSymbol)
    df <- df[order(df$cpms, decreasing = TRUE),]
    where <- which(df$genes %in% signatures[[j]])
    signatures_byCell[location[i],j]  <- 1-mean(where)/nrow(df) #the score is defined as the mean rank of the NKsignature in the ranked cpm list. The score is inverted towards 1, because intuitively higher score then reflects more NK cells
  }
  print(paste("done with signature", j))
}
colnames(signatures_byCell) <- names(signatures)
dim(signatures_byCell)

saveRDS(signatures_byCell,file = "signatures_byCell.rds")


################################
# cell numbers and frequency 
###############################

#count cells per sample
patientStats <- data.frame(patientID = unique(sce$patientID))

patientStats <-  colData(sce)[,c("response", "therapy","patient",  "timepoint", "patientID")]
patientStats <- patientStats[!duplicated(patientStats$patientID),]
patientStats <- data.frame(patientStats)

#count NK cells per sample
NKcells <- data.frame(table(colData(sce)[,c( "isNK", "patientID")]))
NKcells <- NKcells %>% filter(isNK == TRUE) %>% dplyr::select(-isNK)
colnames(NKcells) <- c("patientID", "NKcells_n")
NKcells

table_nrcells <- data.frame(table(colData(sce)$patientID))
colnames(table_nrcells) <- c("patientID", "cells_n")
table_nrcells

#count all different clusters per patient
table_clust <- data.frame(table(colData(sce)[,c("patientID", "clust")]))
table_clust <- table_clust %>% spread(clust, Freq)
colnames(table_clust) <- c("patientID", paste0( "clust", colnames(table_clust)[2:ncol(table_clust)], "_n" ))
table_clust

#combine all numbers
patientStats <- left_join(patientStats, table_nrcells, by = "patientID")
patientStats <- left_join(patientStats, NKcells, by = "patientID")
patientStats <- left_join(patientStats, table_clust, by = "patientID")
patientStats

#calculate frequencies
patientStats <- patientStats %>% 
  mutate(NKcells_f = NKcells_n/cells_n*100,
         clust1_f = clust1_n/cells_n*100,
         clust2_f = clust2_n/cells_n*100,
         clust3_f = clust3_n/cells_n*100,
         clust4_f = clust4_n/cells_n*100,
         clust5_f = clust5_n/cells_n*100,
         clust6_f = clust6_n/cells_n*100,
         clust7_f = clust7_n/cells_n*100,
         clust8_f = clust8_n/cells_n*100,
         clust9_f = clust9_n/cells_n*100,
         clust10_f = clust10_n/cells_n*100,
         clust11_f = clust11_n/cells_n*100,
         clust12_f = clust12_n/cells_n*100,
         clust13_f = clust13_n/cells_n*100)
patientStats
write.table(patientStats, file = "../lists/patientStats.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#summarize these values for response timeoint and therapy
patientStats_summarised <- patientStats %>% group_by(response, timepoint, therapy) %>%
  summarise(cells_n = sum(cells_n),
            NKcells_n = sum(NKcells_n),
            clust1_n = sum(clust1_n),
            clust2_n = sum(clust2_n),
            clust3_n = sum(clust3_n),
            clust4_n = sum(clust4_n),
            clust5_n = sum(clust5_n),
            clust6_n = sum(clust6_n),
            clust7_n = sum(clust7_n),
            clust8_n = sum(clust8_n),
            clust9_n = sum(clust9_n),
            clust10_n = sum(clust10_n),
            clust11_n = sum(clust11_n),
            clust12_n = sum(clust12_n),
            clust13_n = sum(clust13_n))

#calculate frequencies
patientStats_summarised <- patientStats_summarised %>% 
  mutate(NKcells_f = NKcells_n/cells_n*100,
         clust1_f = clust1_n/cells_n*100,
         clust2_f = clust2_n/cells_n*100,
         clust3_f = clust3_n/cells_n*100,
         clust4_f = clust4_n/cells_n*100,
         clust5_f = clust5_n/cells_n*100,
         clust6_f = clust6_n/cells_n*100,
         clust7_f = clust7_n/cells_n*100,
         clust8_f = clust8_n/cells_n*100,
         clust9_f = clust9_n/cells_n*100,
         clust10_f = clust10_n/cells_n*100,
         clust11_f = clust11_n/cells_n*100,
         clust12_f = clust12_n/cells_n*100,
         clust13_f = clust13_n/cells_n*100)
#save to file
write.table(patientStats_summarised, file = "../lists/patientStats_summarised.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file = "../plots/cell_frequency_plots.pdf", width = 8)

patientStats %>% 
  ggplot(aes(x = response, y = NKcells_f, fill = timepoint))  + 
  geom_violin(position=position_dodge(0.8)) +
  geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.5) + 
  ggtitle("NK cell frequency") +
  ylab("frequency NK cells of all immune cells") +
  facet_wrap(~therapy) + theme_bw()

for(i in 1:13){
  cluster <- paste0("clust", i, "_f")  
  p <- patientStats %>% 
    ggplot(aes_string(x = "response", y = cluster, fill = "timepoint"))  + 
    geom_violin(position=position_dodge(0.8)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.5) + 
    ggtitle(paste("cluster",i,"frequency")) +
    ylab(paste("cluster",i,"frequency of all immune cells")) +
    facet_wrap(~therapy) + theme_bw()
  print(p)
}

for(i in 1:13){
  cluster <- paste0("clust", i, "_f")  
  p <- patientStats %>% 
    ggplot(aes_string(x = "response", y = cluster))  + 
    geom_violin(position=position_dodge(0.8)) +
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), dotsize = 0.5) + 
    ggtitle(paste("cluster",i,"frequency")) +
    ylab(paste("cluster",i,"frequency of all immune cells"))  + theme_bw()
  print(p)
}



patientStats %>% filter(therapy == "anti-PD1") %>%
  ggplot(aes(x = NKcells_f*100, y = clust1_f, col = response, pch = timepoint) ) + geom_point() + 
  ggtitle("Frequency cluster 1 DCs vs NK cells")

cormat <- round(cor(patientStats[,grepl("_f", colnames(patientStats)) ] ),2)
head(cormat)
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + ggtitle("correlation between cluster frequencies") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))

dev.off()

#############
# Signatures by cells
############




signatures_byCell <- readRDS(file = "signatures_byCell.rds")

#add signature score to each cell colData
colData(sce) <- colData(sce)[,!(colnames(colData(sce))  %in% colnames(signatures_byCell))    ]

colData(sce) <- cbind(colData(sce), signatures_byCell)

pdf(file = "../plots/cells_signatures.pdf", width = 8)

sce$UMAP1 <- reducedDim(sce, "UMAP")[,1]
sce$UMAP2 <- reducedDim(sce, "UMAP")[,2]
sce$PCA1 <- reducedDim(sce, "PCA")[,1]
sce$PCA2 <- reducedDim(sce, "PCA")[,2]
sce$TSNE1 <- reducedDim(sce, "TSNE")[,1]
sce$TSNE2 <- reducedDim(sce, "TSNE")[,2]

dat <- data.frame(as.data.frame(colData(sce)), "CCL5" = matrix(logcounts(sce[rowData(sce)$GeneSymbol == "CCL5",])))

p <- ggplot(dat, aes(x = TSNE1, y = TSNE2, col = NK2_Signature_Zilionis2019)) +
  geom_point(  alpha = 0.4) + 
  theme_bw() + 
  scale_color_gradient( low="blue",high="red", space ="Lab" ) +
  ggtitle("PCA NK2")
print(p)
p <- ggplot(dat, aes(x = TSNE1, y = TSNE2, col = NK1_Signature_Zilionis2019)) +
  geom_point(  alpha = 0.4) + 
  theme_bw() + 
  scale_color_gradient( low="blue",high="red", space ="Lab" ) +
  ggtitle("PCA NK1")
print(p)




ggplot(dat, aes(x = NK1_Signature_Zilionis2019, y= NK2_Signature_Zilionis2019 ,col = TissueResidentNK_Marquart2019)) + 
  geom_point() + scale_color_gradient( low="blue",high="red", space ="Lab" ) +
  theme_bw()
ggplot(dat, aes(x = NK1_Signature_Zilionis2019, y= NK2_Signature_Zilionis2019 ,col = CD56bright_CD16neg_Hanna2004)) + 
  geom_point() + scale_color_gradient( low="blue",high="red", space ="Lab" ) +
  theme_bw()
ggplot(dat, aes(x = NK1_Signature_Zilionis2019, y= NK2_Signature_Zilionis2019 ,col = CCL5)) + 
  geom_point() + scale_color_gradient( low="blue",high="red", space ="Lab" ) +
  theme_bw()

ggplot(dat, aes(x = NK1_Signature_Zilionis2019, y= CD56dim_CD16pos_Hanna2004 ,col =  NK2_Signature_Zilionis2019 )) + 
  geom_point() +
  theme_bw()
ggplot(dat, aes(x = CD56bright_CD16neg_Hanna2004, y= NK2_Signature_Zilionis2019 ,col = isNK)) + 
  geom_point() +
  theme_bw()

ggcells(sce, aes(x = NK2_Signature_Zilionis2019, y= CD3E)) + 
  geom_point() +
  theme_bw()
ggcells(sce, aes(x = NK1_Signature_Zilionis2019, y= CD3E)) + 
  geom_point() +
  theme_bw()

ggplot(dat_isNK, aes(x = NK1_Signature_Zilionis2019, y= NK2_Signature_Zilionis2019, col = CD56bright_CD16neg_Hanna2004)) + 
  geom_point() +
  theme_bw()
ggplot(dat_isNK, aes(x = CD56bright_CD16neg_Hanna2004, y= NK2_Signature_Zilionis2019, col = TissueResidentNK_Marquart2019)) + 
  geom_point() +
  theme_bw()

dev.off()

#mean per patients
clust1_signatures_selectedcells <- data.frame(colData(sce)) %>% dplyr::filter(clust == 1) %>% dplyr::select(patient, response, therapy, timepoint, names(signatures))
 
clust1_signatures_perPatient <- clust1_signatures_selectedcells %>%
  group_by(patient, timepoint) %>% 
  mutate(meanNK2 = mean(NK2_Signature_Zilionis2019),
            meanNK1 = mean(NK1_Signature_Zilionis2019),
            meanCD56bright = mean(CD56bright_CD16neg_Hanna2004),
            meanCD56dim = mean(CD56dim_CD16pos_Hanna2004),
            meanTissueRes = mean(TissueResidentNK_Marquart2019)
                                  ) %>%
  ungroup() %>% 
  dplyr::select(patient, response, therapy, timepoint, meanNK2, meanNK1, meanCD56bright, meanCD56dim) %>%
  distinct(patient, timepoint, .keep_all = TRUE)

clust1_signatures_perPatient <- clust1_signatures_perPatient %>% mutate(NK2_vs_NK1 = meanNK2/meanNK1) 
clust1_signatures_perPatient <- clust1_signatures_perPatient %>% mutate(bright_vs_dim = meanCD56bright/meanCD56dim) 
clust1_signatures_perPatient$patientID <- paste0(clust1_signatures_perPatient$timepoint, "_", clust1_signatures_perPatient$patient)
# clust1_signatures_perPatient <- clust1_signatures_perPatient %>% filter(therapy == "anti-PD1" | therapy== "anti-CTLA4+PD1")


meanSignatures <- c("meanNK2", "meanNK1", "meanCD56bright" ,"meanCD56dim",  "bright_vs_dim", "NK2_vs_NK1")

for(i in 1:length(meanSignatures)){
  signat <- meanSignatures[i]
  dat <- clust1_signatures_perPatient %>% filter(timepoint == "Pre") %>% dplyr::select(response, signat)
  colnames(dat) <- c("response", "value")
  t <- t.test(value ~ response, data = dat)$p.value
  w <- wilcox.test(value ~ response, data = dat)$p.value
  p <- ggplot(dat, aes_string(x= "response", y="value")) +
    geom_boxplot() +
    geom_point() + 
    theme_bw() +
    ggtitle(paste(signat, sprintf("\n t-test p-value= %.2e",t ) )) 
  print(p)
}

summary(patientStats$clust1_n > 1)
enoughNKs <- patientStats[patientStats$clust1_n > 1,]
enoughNKs

clust1_signatures_perPatient_include <- clust1_signatures_perPatient[clust1_signatures_perPatient$patientID %in% enoughNKs$patientID,]
for(i in 1:length(meanSignatures)){
  signat <- meanSignatures[i]
  dat <- clust1_signatures_perPatient_include %>% dplyr::select(response, signat, timepoint) 
  colnames(dat) <- c("response", "value", "timepoint")
  t <- t.test(value ~ response, data = dat)$p.value
  w <- wilcox.test(value ~ response, data = dat)$p.value
  p <- ggplot(dat, aes_string(x= "response", y="value")) +
    geom_boxplot() +
    geom_jitter(width= 0.1) + 
    theme_bw() +
    ggtitle(paste(signat, sprintf("\n t-test p-value= %.2e",t), sprintf("\n wilcoxon-test p-value= %.2e", w )) ) 
  print(p)
}

write.table(clust1_signatures_perPatient_include, file = "../lists/clust1_signatures_perPatient_morethan1cell.txt", sep = "\t", row.names= FALSE, quote = FALSE)
