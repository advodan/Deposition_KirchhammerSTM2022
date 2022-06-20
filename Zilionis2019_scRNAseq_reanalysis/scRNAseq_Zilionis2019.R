# This is the re-analysis of scRNAseq data from Zilionis et al. 2019
# Parts of the manuscript of Kirchhammer et al STM 2022. Performed by Marcel P. Trefny in 2021/22



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
# BiocManager::install("virids", ask = FALSE)
#BiocManager::install("Matrix", ask = FALSE)
# BiocManager::install("Seurat", ask = FALSE)
# BiocManager::install("cli", ask = FALSE)
library(Seurat)
library(Matrix)
library(batchelor)
library(BiocSingular)
library(BiocNeighbors)
library(Rtsne)
library(rsvd)
library(tidyr) 
library(dplyr)
# library(cellassign)
library(igraph)
#library(TENxPBMCData)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Rtsne)
library(viridis)


# library(rowr)
########################################################################

#Data Loading

path <- "/Volumes/CIMM$/CIMM/Projects_Human_Team/Marcel_Nicole/AdV5_IL12_Project/Zilionis2019/scRNAseq/" #iMAC Hebelstrasse 
# path <- "/Users/marceltrefny/Documents/work/GSE121861" #macbook 
# 

setwd(path)

cellinfo <- read.table(file = "input/GSE127465_human_cell_metadata_54773x25.tsv", sep = "\t", header = TRUE)
head(cellinfo)
dim(cellinfo)
#54773    25
#####################################
# This whole part is not necessary to run. this was used to make the sparse matrix only with the NK cells
# file GSE127465_human_counts_normalized_54773x41861 can be downloaded from GEO and was not deposited here due to size. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127465
########################################
# NK_cells_index <- which(cellinfo$Minor.subset %in% c("tNK1", "tNK2"))
# length(NK_cells_index)
# 
# counts_dg <- read.table(file = "input/GSE127465_human_counts_normalized_54773x41861.mtx", skip = 3, sep = " ", header = FALSE)
# counts_dg_NK <- counts_dg %>% filter(V1 %in% NK_cells_index)
# colnames(counts_dg_NK) <- c("cellID", "rowID", "value")
# head(counts_dg_NK)
# 
# #clear up the memory
# rm(counts_dg)
# 
counts_dg_NK <- readRDS(file = "input/GSE127465_human_counts_normalized_NK_only.rds")
NK_reindexing <- data.frame("cellID" = unique(counts_dg_NK$cellID), "cellID_new" = seq(1,length(unique((counts_dg_NK$cellID))), 1))
dim(NK_reindexing)
head(NK_reindexing)
cellinfo$cellID <- seq(1,nrow(cellinfo), 1)
head(cellinfo)
cellinfo_NK <- cellinfo[cellinfo$Minor.subset %in% c("tNK1", "tNK2"),]
cellinfo_NK <- left_join(cellinfo_NK, NK_reindexing, by = "cellID")
cellinfo_NK <- cellinfo_NK[order(cellinfo_NK$cellID_new),]
dim(cellinfo_NK)
head(cellinfo_NK)
 
# counts_dg_NK <- left_join(counts_dg_NK, NK_reindexing, by = "cellID")
# head(counts_dg_NK)
# 
# 
# cellinfo_NK$cellID_new <- NK_reindexing$cellID_new
# 
# counts_NK_sparse <- sparseMatrix(i = as.numeric(0), j = as.numeric(0), dims = c(max(counts_dg_NK$rowID), max(counts_dg_NK$cellID_new)))
# counts_NK_sparse <- Matrix(data = 0 , nrow = max(counts_dg_NK$rowID), ncol=  max(counts_dg_NK$cellID_new), sparse = TRUE )
# for(i in 1:nrow(counts_dg_NK)){
#   values <- counts_dg_NK[i,]
#   counts_NK_sparse[values[,"rowID"], values[,"cellID_new"] ] <- values[,"value"]
#   }
# 
# saveRDS(counts_NK_sparse,file = "input/GSE127465_human_counts_normalized_NK_only_sparseMatrix.rds")

counts_NK_sparse <- readRDS(file = "input/GSE127465_human_counts_normalized_NK_only_sparseMatrix.rds")


genes <- read.table("input/GSE127465_gene_names_human_41861.tsv",  sep = "\t",  check.names= FALSE, header = FALSE, stringsAsFactors = FALSE )
head(genes)

sce <- SingleCellExperiment(assays=SimpleList(counts=counts_NK_sparse), rowData = data.frame("GeneSymbol" = genes$V1),
                            colData=cellinfo_NK)
rownames(sce) <- make.names(rowData(sce)$GeneSymbol) #need to make_names because of genes that start with a number
sce

saveRDS(sce,file = "input/GSE127465_single_cell_experiment_NK_only.rds")
sce <- readRDS(sce,file = "input/GSE127465_single_cell_experiment_NK_only.rds")


colData(sce)

dim(assay(sce, "counts"))

head(rowData(sce))
head(colData(sce))
nrow(counts(sce))
ncol(counts(sce))

mean(counts(sce) == 0)
class(counts(sce))
counts(sce)[1:10, 1:10]

## Total number of detected transcripts
#remove all genes that were not detected at all
sum(rowSums(counts(sce)) > 0)
sce <- sce[rowSums(counts(sce)) > 0,]
nrow(counts(sce)) #21994

#no filtering required as this these are prefiltered cells

########################################################################
#### Normalization
#look at the different library sizes of each cell with some patients having much more transcripts per cell sequenced
sce <- computeSumFactors(sce, min.mean = 0.1) #normalization of library sizes 
sce <- logNormCounts(sce)
assayNames(sce)

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


#pdf(file = "../plots/reducedDimensions_PCA_tSNE_UMAP.pdf", width = 8)

colnames(colData(sce))

plotReducedDim(sce, "PCA", colour_by = "Patient")
plotReducedDim(sce, "PCA", colour_by = "Tissue")
plotReducedDim(sce, "PCA", colour_by = "Major.cell.type")
plotReducedDim(sce, "PCA", colour_by = "Minor.subset")
plotReducedDim(sce, "PCA", colour_by = "CCL5")
plotReducedDim(sce, "PCA", colour_by = "Tissue")
set.seed(123)

#perform tSNE
# sce <- runTSNE(sce, dimred = "PCA",  perplexity = 30)
# reducedDimNames(sce)

sce <- runUMAP(sce, dimred = "PCA", name = "UMAP")
reducedDimNames(sce)

#add the UMAP coordinates to the colData
sce$UMAP1 <- reducedDim(sce, "UMAP")[,1]
sce$UMAP2 <- reducedDim(sce, "UMAP")[,2]
colData(sce)

pdf(file = "./plots/Zilionis_UMAP_NKsubsets.pdf", width = 6, height = 4.5)
genes_of_interest <- c("CCL5", "FGFBP2", "FCGR3A",  "ITGAE", "ITGA1", "ENTPD1", "CD69", "FCGR3A")

for(i in 1:length(genes_of_interest)){
  gene <- genes_of_interest[i]
  if(gene %in% rowData(sce)$GeneSymbol){ #only do if the gene is found

    dat <- data.frame(as.data.frame(colData(sce)), "gene_logcounts" = matrix(logcounts(sce[rowData(sce)$GeneSymbol == gene,])))
    p <- ggplot(dat, aes(x = UMAP1, y = UMAP2, size = exp(gene_logcounts),  col = gene_logcounts)) +
      geom_point(  alpha = 0.9) + 
      labs(size = paste(gene, "counts"), col = paste0("log(", gene," counts)")) +
      theme_bw() + 
      scale_color_gradient( low="blue",high="red", space ="Lab" ) +
      ggtitle(paste("UMAP", gene))
    print(p)
  }
}

dev.off()

#This plot was used in Figure 4I
dat <- data.frame(as.data.frame(colData(sce)), "gene_logcounts" = matrix(logcounts(sce[rowData(sce)$GeneSymbol == "CCL5",])))
p <- ggplot(dat, aes(x = UMAP1, y = UMAP2, size = exp(gene_logcounts),  col = gene_logcounts)) +
  geom_point(  alpha = 0.9) + 
  labs(size = "CCL5 counts", col = "log(CCL5 counts)") +
  theme_bw() + scale_color_viridis() +
  ggtitle(paste("UMAP", "CCL5"))
print(p)

#this plot was used in Figure 4H
pdf(file = "./plots/Zilionis_UMAP_NKsubsets_NK1_NK2.pdf", width = 6, height = 5)
  
  p <- ggplot(dat, aes(x = UMAP1, y = UMAP2, col = Minor.subset)) +
    geom_point(size = 2.3,  alpha = 0.7) +
    theme_bw() + scale_color_manual(values=c("#90C82B", "#435D14")) +
   ggtitle(paste("UMAP of NK cells in Zilionis -- ", "NK1 vs NK2"))
  print(p)

dev.off()


#################
signatures_files <- list.files("./signatures/")

signatures <- sapply(signatures_files, function(x) read.table(paste0("./signatures/", x) , sep = "\t", stringsAsFactors = FALSE  )$V1, USE.NAMES = TRUE)
signatures

sapply(signatures, function(x) summary(x %in% rowData(sce)$GeneSymbol))
sapply(signatures, function(x) x[!x %in% rowData(sce)$GeneSymbol]) #which are not in the dataset

signatures <- sapply(signatures, function(x) x[x %in% rowData(sce)$GeneSymbol]) #only take those that are found
names(signatures) <- gsub(".txt","", names(signatures)) #give them good names
signatures
# #this does not have to be done each time, just for new signatures
# #for each cell, calculate signature score for all signatures
# signatures_byCell <- data.frame(matrix(ncol=length(signatures), nrow= nrow(colData(sce))))
# for (j in 1:length(signatures)){ #for each signature
#   for(i in 1:nrow(colData(sce))) { # for each patient
#     df <- data.frame(cpms = logcounts(sce)[,i], genes = rowData(sce)$GeneSymbol)
#     df <- df[order(df$cpms, decreasing = TRUE),]
#     where <- which(df$genes %in% signatures[[j]])
#     signatures_byCell[i,j]  <- 1-mean(where)/nrow(df) #the score is defined as the mean rank of the NKsignature in the ranked cpm list. The score is inverted towards 1, because intuitively higher score then reflects more NK cells 
#   }
# }
# colnames(signatures_byCell) <- names(signatures)
# 
# saveRDS(signatures_byCell,file = "signatures_byCell.rds")

signatures_byCell <- readRDS(file = "signatures_byCell.rds")

dat <- data.frame(  colData(sce), as.matrix(t(logcounts(sce[rowData(sce)$GeneSymbol == "CCL5",]))))
write.table(dat, file = "./lists/cellData_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)


#add signature score to each cell colData
colData(sce) <- cbind(colData(sce), signatures_byCell)
sce.seurat <- as.Seurat(sce, counts = "counts", data = "logcounts")


pdf(file = "./plots/Zilionis_SeuratDotPlots.pdf", width = 13, height = 4)

DotPlot(sce.seurat, features = signatures[["NK2_Signature_Zilionis2019"]], cols = c("blue", "red"), dot.scale = 8, group.by="Minor.subset") +
  RotatedAxis() + theme(axis.text.x = element_text(size = 6))
DotPlot(sce.seurat, features = signatures[["NK1_Signature_Zilionis2019"]], cols = c("blue", "red"), dot.scale = 8, group.by="Minor.subset") +
  RotatedAxis() + theme(axis.text.x = element_text(size = 6))


DotPlot(sce.seurat, features = signatures[["CD56bright_CD16neg_Hanna2004"]], cols = c("blue", "red"), dot.scale = 8, group.by="Minor.subset") +
  RotatedAxis() + theme(axis.text.x = element_text(size = 6))

DotPlot(sce.seurat, features = signatures[["CD56dim_CD16pos_Hanna2004"]], cols = c("blue", "red"), dot.scale = 8, group.by="Minor.subset") +
  RotatedAxis() + theme(axis.text.x = element_text(size = 6))

DotPlot(sce.seurat, features = signatures[["TissueResidentNK_Marquart2019"]], cols = c("blue", "red"), dot.scale = 8, group.by="Minor.subset") +
  RotatedAxis() + theme(axis.text.x = element_text(size = 6))
dev.off()


pdf(file = "./plots/Zilionis_NK_signatures_plots.pdf")

plotReducedDim(sce, "UMAP", colour_by = "CD56bright_CD16neg_Hanna2004")
plotReducedDim(sce, "UMAP", colour_by = "CD56dim_CD16pos_Hanna2004")
plotReducedDim(sce, "UMAP", colour_by = "NK1_Signature_Zilionis2019")
plotReducedDim(sce, "UMAP", colour_by = "NK2_Signature_Zilionis2019")
plotReducedDim(sce, "UMAP", colour_by = "TissueResidentNK_Marquart2019")
plotReducedDim(sce, "UMAP", colour_by = "clust")
plotReducedDim(sce, "UMAP", colour_by = "CCL5")
plotReducedDim(sce, "UMAP", colour_by = "clust")



for(i in 1:length(signatures)){
  signat <- names(signatures)[i]
  dat <- colData(sce)[, c(signat, "Minor.subset") ]
  colnames(dat) <- c("values", "Minor.subset")
  t <- t.test(values ~ Minor.subset, data = dat)$p.value
  
  p <- ggcells(sce, mapping=aes_string(x= "Minor.subset", y = signat, fill ="Minor.subset"  ) ) + 
    geom_violin(position=position_dodge(0.8))  +
    geom_boxplot(width= 0.1) +
    ggtitle(paste(signat, sprintf("\n t-test p-value= %.2e",t ) )) +
    theme_bw() 
  print(p)
}

genes_of_interest <- c("CCL5", "FGFBP2", "ITGAE", "ITGA1", "ENTPD1", "CD69", "FCGR3A")

for(i in 1:length(genes_of_interest)){
  gene <- genes_of_interest[i]
  if(gene %in% rowData(sce)$GeneSymbol){ #only do if the gene is found
  dat <- data.frame( as.matrix(t(logcounts(sce[rowData(sce)$GeneSymbol == gene,]))), "Minor.subset" = colData(sce)[, "Minor.subset"] )
  colnames(dat) <- c("value", "Minor.subset")
  t <- t.test(value ~ Minor.subset, data = dat)$p.value
  
  p <- ggplot(dat, mapping=aes_string(x= "Minor.subset", y = "value", fill ="Minor.subset"  ) ) + 
    geom_violin(position=position_dodge(0.8))  +
    geom_boxplot(width= 0.1) +
    ylab(paste("logcount", gene)) +
    ggtitle(paste(gene, sprintf("\n t-test p-value= %.2e",t ) )) +
    theme_bw() 
  print(p)
  }
}



genes_of_interest <- c("CCL5", "IL12A", "FGFBP2", "ITGAE", "ITGA1", "FCGR3A")

for(i in 1:length(genes_of_interest)){
  gene <- genes_of_interest[i]
  if(gene %in% rowData(sce)$GeneSymbol){ #only do if the gene is found
    dat <- data.frame( as.matrix(t(logcounts(sce[rowData(sce)$GeneSymbol == gene,]))), "clust" = colData(sce)[, "clust"] )
    colnames(dat) <- c("value", "clust")

    p <- ggplot(dat, mapping=aes_string(x= "clust", y = "value", fill ="clust"  ) ) + 
      geom_violin(position=position_dodge(0.8))  +
      geom_boxplot(width= 0.1) +
      ylab(paste("logcount", gene)) +
      ggtitle(paste(gene)) +
      theme_bw() 
    print(p)
  }
}


dev.off()
