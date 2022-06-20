# Script to analyse Riaz et al RNaseq Dataset to find correlation of CCL5 and IL12 to response to PD-1
# performed by Marcel P. Trefny for manuscript of Kirchhammer et al STM 2020

# BiocManager::install(c("GOstats", "limma","TxDb.Hsapiens.UCSC.hg38.knownGene", "NbClust", "edgeR", "locfit","RColorBrewer", "grDevices", "NMF", "org.Hs.eg.db","ComplexHeatmap", "circlize","statmod") )
# BiocManager::install("dplyr")
library(dplyr)
library(limma)
library(edgeR)
library(locfit)
library(RColorBrewer)
library(grDevices)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(statmod)
library(GOstats)
library(NbClust)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(viridis)
library(ggplot2)
library(stringr)
library(circlize)
library(tidyr)
library(GeneOverlap)
library(clinfun)
library(reshape2)
# BiocManager::install("ggfortify")
# BiocManager::install("clinfun")
#install.packages("corrr")
library(corrr)
myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
set.seed(123)
setwd("/Volumes/CIMM$/CIMM/Projects_Human_Team/Marcel_Nicole/AdV5_IL12_Project/Riaz2017_RNAseq_reanalysis_deposition//") #change to file path

# data was downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061

dgelist <- read.csv(file = "input/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv")
genes <- dgelist$X
genes <- data.frame("ENTREZID" =  as.character(genes))
genes.annot <- unique(AnnotationDbi::select(org.Hs.eg.db, keys=as.character(genes$ENTREZID), columns='SYMBOL', keytype='ENTREZID'))
genes <- left_join(genes, genes.annot, by = "ENTREZID")
head(genes)
tail(genes)

counts <- dgelist %>% dplyr::select(-X)
head(dgelist)

l <- str_split(colnames(counts), "_")
df <- data.frame(matrix(unlist(l), nrow=length(l), byrow=T))
df
colnames(df) <- c("Patient", "timepoint", "ID")

table(df$Patient)

patient_charac <- read.table(file = "input/patients_characteristics.txt", header= TRUE, sep = "\t")
patient_charac

samples <- left_join(df, patient_charac, by = "Patient" )
samples


# exclude all patients that do not have a paired pre and post sample
t <- table(samples$Patient, samples$timepoint)
exclude <- ( samples$Patient %in% c( unique(samples$Patient)[ !(t[,"On"] == 1 & t[,"Pre"] == 1 ) ], samples$Patient[samples$Response == "NE"] ))
counts <- counts[,!exclude]
samples<- samples[!exclude,]
head(counts)

#explore the data
colSums(counts)
colMeans(counts)





#prepare data for initial exploration
#pseudocount #for now just 5 -> determine optimal pseudocount
cpms <- edgeR::cpm(counts, prior.count = 8, log=TRUE)

#Quality control Plots
par(mfrow=c(1,1))
boxplot(log2(counts+1),ylab="log2(counts+1)", main = "Counts Unnormalized" )
boxplot(cpms, ylab="cpms", cex.axis = 0.7, main = "Boxplot Rawdata log cpms") # 

par(mfrow=c(2,2))
for (i in 1:(ncol(cpms))) {
  hist(cpms[,i], breaks = 20, main = paste("Histogram log2(cpm) \n", (colnames(cpms))[i]), cex.main=0.8, xlab="log2(cpm)", ylab="number of genes")
}

#Construct the DGE list
y <- DGEList(counts = counts, genes = genes, samples = samples)

#filter out genes that our not expressed in at least 2 samples at 1CPM  ------------------------
keep <- (rowSums(cpms>0) >= 6 & !is.na(genes$SYMBOL))
summary(keep) 

y <- y[keep,]

d <- duplicated(y$genes$SYMBOL)
y <- y[!d,]
nrow(y) #18557

#normalize for library sizey
y <- calcNormFactors(y)
y$samples$group <- as.factor(y$samples$timepoint)
y$samples$timepoint <- factor(y$samples$timepoint, levels = c("Pre", "On"))
y$samples$Patient <- as.factor(y$samples$Patient)
y$samples$benefit <- factor(ifelse(y$samples$Response %in% c("PR", "CR", "SD"), "CB", "PD")  , levels= c("PD", "CB"))  # CB = clinical benefit, PD = no benefit = progressive disease
y$samples 



par(mfrow=c(1,1))
cpms_filtered <- edgeR::cpm(y, prior.count = 8, log=TRUE, normalized.lib.sizes = T) #log cpms corrected
rownames(cpms_filtered) <- y$genes$SYMBOL
plotDensities(cpms_filtered, group = y$samples$SampleGroup)



limma::plotMDS(y, gene.selection = "common", col = myPalette[as.numeric(y$samples$group)], main = "MDS Analysis RNASeq by Treatment" , xlim = c(-3, 3), ylim = c(-1.5, 3.5))
legend("right", legend = levels(y$samples$group),title ="Treatment" , text.col = myPalette[1:4])

#PCA
summary(apply(cpms_filtered, 1, var) == 0) ## Any gene with no variability?
pca1 <- prcomp(t(cpms_filtered), scale = T)
summary(pca1) ## 59% variance on PCs 1:3
pdf(paste0("./plots/PCA_loading.pdf"), width=7)
{
  barplot(summary(pca1)$importance[2,], ylab="Proportion of Variance", main="PCA Components", cex.axis=0.8, cex.names = 0.5)
  dev.off()
}
# loadings: coordinates of PCs in original coordinate system

loadings <- pca1$rotation
## scores: coordinates of scaled data (observations) projected onto the PCs.
scores <- pca1$x

plotPCA <- function(n=1, m=2) {
  col.v <- viridis(length(y$samples$Patient))[as.integer(y$samples$Patient)]
  pchs <- 21
  pdf(paste0("./plots/PCA_", n,"vs", m ,".pdf"), width=7)
  {
    par(mar=c(5,5, 5,5) )
    plot(scores[,n], scores[,m], xlim=c(min(scores[,n])-30, max(scores[,n])+30), ylim=c(min(scores[,m])-30, max(scores[,m])+30), xlab=paste("PC", n, ": ", round(summary(pca1)$importance[2,n],3)*100, "% variance explained", sep=""), ylab=paste("PC", m, ": ", round(summary(pca1)$importance[2,m],3)*100, "% variance explained", sep=""), col="black", bg=col.v, lwd = 0.5, cex= 1.8, cex.lab=1.3, cex.main = 1.5,  pch=pchs, main = "Principal Component Analysis")
    legend("bottomright", levels(y$samples$group), pch=21, cex=1.2, col="black", pt.bg=col.v)
    #legend("bottomright", levels(y$samples$sex), cex=0.95, col="black")}
    dev.off()
  }
}
plotPCA(1,2)
plotPCA(2,3)
plotPCA(1,3)

#####

#define groups to look at
data.frame(Samplename=colnames(y), y$samples$Patient, y$samples$timepoint) # check if all samples are annotated correctly


#Build Design Matrix, make a batch correction for donor
design <- model.matrix(~ y$samples$Patient) #starting with a paired sample model
rownames(design) <- colnames(y)

PD.changes <- y$samples$timepoint == "On" & y$samples$benefit == "PD"
CB.changes <- y$samples$timepoint == "On" & y$samples$benefit == "CB"

design <- cbind(design, PD.changes, CB.changes)
design_nonpaired <- model.matrix(~ y$samples$timepoint + y$samples$benefit)

colnames(design) <- gsub("y\\$samples\\$", "", colnames(design))
colnames(design_nonpaired) <- gsub("y\\$samples\\$", "", colnames(design_nonpaired))
design
design_nonpaired
#Build GLM
y <- estimateDisp(y, design)

plotBCV(y, main="BCV Plot paired") #Plotting Biological Coefficient of Variation



fit <- glmQLFit(y, design, prior.count=8)
fit_glmFit <- glmFit(y, design, prior.count = 8)
pcut=0.3
logfcut = 0


#define contrasts
mycontrasts <- list("PDchanges" =  c(rep(0, ncol(design)-2),1,0),"CBchanges" =  c(rep(0, ncol(design)-2),0,1) , "diffChanges" = c(rep(0, ncol(design)-2),-1,1)  )  #last one means genes that are more or less regulated in clinical benefit than in no benefit

lrt_all <- sapply(names(mycontrasts),function(x) glmQLFTest(fit,contrast=mycontrasts[[x]] ), simplify=FALSE, USE.NAMES = TRUE )


tables_all_unsorted <- sapply(names(lrt_all), function(x)  topTags(lrt_all[[x]],adjust.method="BH", n = Inf, sort.by = "none")$table,simplify=FALSE, USE.NAMES = TRUE)
tables_all <- sapply(names(lrt_all), function(x)  topTags(lrt_all[[x]],adjust.method="BH",  n = Inf, sort.by = "PValue")$table, simplify=FALSE, USE.NAMES = TRUE)
tables_signif <- sapply(names(lrt_all), function(x) tables_all[[x]][abs(tables_all[[x]]$logFC) > logfcut & tables_all[[x]]$FDR < pcut, ],simplify = FALSE, USE.NAMES = TRUE )

sapply(tables_all, function(x) head(x, n=20),simplify=FALSE, USE.NAMES = TRUE)
sapply(names(lrt_all), function(x) table(decideTestsDGE(lrt_all[[x]], p.value = pcut, adjust.method="BH", lfc = logfcut)), simplify = TRUE, USE.NAMES = TRUE)
de_all <- sapply(names(lrt_all), function(x) decideTestsDGE(lrt_all[[x]], p.value = pcut, adjust.method="BH", lfc = logfcut), simplify = FALSE, USE.NAMES = TRUE)

lrt_all_glmFit <- sapply(names(mycontrasts),function(x) glmLRT(fit_glmFit,contrast=mycontrasts[[x]] ), simplify=FALSE, USE.NAMES = TRUE )
tables_all_unsorted_glmFit <- sapply(names(lrt_all_glmFit), function(x)  topTags(lrt_all_glmFit[[x]],adjust.method="BH", n = Inf, sort.by = "none")$table,simplify=FALSE, USE.NAMES = TRUE)
tables_all_glmFit <- sapply(names(lrt_all_glmFit), function(x)  topTags(lrt_all_glmFit[[x]],adjust.method="BH",  n = Inf, sort.by = "PValue")$table, simplify=FALSE, USE.NAMES = TRUE)
sapply(tables_all_glmFit, function(x) head(x, n=20),simplify=FALSE, USE.NAMES = TRUE)


sapply(names(tables_all), function(x) write.table(tables_all[[x]],paste0("./lists/", x,"_dgelist_full.txt"), sep ="\t", quote=FALSE, row.names = FALSE))
sapply(names(tables_all_glmFit), function(x) write.table(tables_all_glmFit[[x]],paste0("./lists/", x,"_dgelist_glm_full.txt"), sep ="\t", quote=FALSE, row.names = FALSE))

pdf(paste0("./plots/MA_Volcanos.pdf"), width=10)
par(mfrow=c(1,2))
par(mar=c(6,4.5, 3,3.5) )
#MA Plots
plots <- lapply(names(lrt_all),function (x) {plotSmear(lrt_all[[x]], de.tags=rownames(y)[as.logical(de_all[[x]])], main=x,font.main=1, sub = paste(sprintf("Down = %i",summary(de_all[[x]])[1]),sprintf("Up = %i",summary(de_all[[x]])[3])))
  abline(h=c(-logfcut, logfcut), col="black")})

par(mfrow=c(1,2))

#Volcano Plots
plots <- lapply(names(tables_all_unsorted),function (x) {plot(tables_all_unsorted[[x]]$logFC,-log10(tables_all_unsorted[[x]]$FDR), xlab="logFC",ylab="-log10(FDR)", pch = 16, cex=0.3, col=ifelse(abs(de_all[[x]])>0, "#FF2600", 1), main=x, font.main=1,  sub = paste(sprintf("Down = %i",summary(de_all[[x]])[1]),sprintf("Up = %i",summary(de_all[[x]])[3])))
  #abline(v=c(-logfcut, logfcut), h=c(-log10(pcut)), col="black")
})

dev.off()


#########
# Boxplots of cpms for selected genes
make_boxplots <- function(name = "empty", selection = c("IL12")) {
pdf(paste0("./plots/Boxplot_", name,".pdf"), width=8)
par(mfrow=c(2,2))
par(mar=c(8,5, 3,3) )
genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% selection,]
for (i in 1:length(selection)) {
  if (length(cpms_filtered[rownames(cpms_filtered) %in% selection[i],]) > 0 ){
    
    genes_of_interest_annot <-  data.frame("cpms" = genes_of_interest[selection[i],], "benefit" = y$samples$benefit, "timepoint" <- y$samples$timepoint, "patient" =  y$samples$Patient, "response" = factor(y$samples$Response, levels= c("PD", "SD", "PR", "CR")))
    
    #Make plots for benefit vs expression pre and post, then switch x and fill
    p <- ggplot(genes_of_interest_annot, aes(x=benefit, y=cpms, fill = timepoint ) ) + geom_boxplot(position=position_dodge(0.8)) + scale_fill_manual(values=c("#FFA100", "#FE4F00")) +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position=position_dodge(0.8)) +
      theme(panel.background = element_rect(fill= "white", linetype = "solid", color = "black", size = 1.5), axis.text.x = element_text(size=18, angle =45, hjust  = 1),  axis.text.y = element_text( size=18), plot.title = element_text(size= 22, hjust= 0.5 ), axis.title.y = element_text(size = 18), legend.text = element_text(size = 18), legend.key.size = unit(1.1, "cm"), legend.title = element_text(size = 18)) +
      labs(title = (paste0("Expression of ", selection[i] )), x ="", y = "log2 cpms", fill = "timepoint")
    
    plot(p)
    
    p <- ggplot(genes_of_interest_annot, aes(x=timepoint, y=cpms, fill = benefit ) ) + geom_boxplot(position=position_dodge(0.8)) + scale_fill_manual(values=c("#03c6fc", "#0a83c9")) +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position=position_dodge(0.8)) +
      theme(panel.background = element_rect(fill= "white", linetype = "solid", color = "black", size = 1.5), axis.text.x = element_text(size=18, angle =45, hjust  = 1),  axis.text.y = element_text( size=18), plot.title = element_text(size= 22, hjust= 0.5 ), axis.title.y = element_text(size = 18), legend.text = element_text(size = 18), legend.key.size = unit(1.1, "cm"), legend.title = element_text(size = 18)) +
      labs(title = (paste0("Expression of ", selection[i] )), x ="", y = "log2 cpms", fill = "timepoint")
    
    plot(p)
    
    
    #Make plots for benefit vs expression pre and post, then switch x and fill
    p <- ggplot(genes_of_interest_annot, aes(x=response, y=cpms, fill = timepoint ) ) + geom_boxplot(position=position_dodge(0.8)) + scale_fill_manual(values=c("#FFA100", "#FE4F00")) +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position=position_dodge(0.8)) +
      theme(panel.background = element_rect(fill= "white", linetype = "solid", color = "black", size = 1.5), axis.text.x = element_text(size=18, angle =45, hjust  = 1),  axis.text.y = element_text( size=18), plot.title = element_text(size= 22, hjust= 0.5 ), axis.title.y = element_text(size = 18), legend.text = element_text(size = 18), legend.key.size = unit(1.1, "cm"), legend.title = element_text(size = 18)) +
      labs(title = (paste0("Expression of ", selection[i] )), x ="", y = "log2 cpms", fill = "timepoint")
    
    plot(p)
    
    p <- ggplot(genes_of_interest_annot, aes(x=timepoint, y=cpms, fill = response ) ) + geom_boxplot(position=position_dodge(0.8)) + scale_color_manual(values=c("#b03a2e", "#dc7633", "#3498db", "#229954"))  +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, position=position_dodge(0.8)) +
      theme(panel.background = element_rect(fill= "white", linetype = "solid", color = "black", size = 1.5), axis.text.x = element_text(size=18, angle =45, hjust  = 1),  axis.text.y = element_text( size=18), plot.title = element_text(size= 22, hjust= 0.5 ), axis.title.y = element_text(size = 18), legend.text = element_text(size = 18), legend.key.size = unit(1.1, "cm"), legend.title = element_text(size = 18)) +
      labs(title = (paste0("Expression of ", selection[i] )), x ="", y = "log2 cpms", fill = "timepoint")
    
    plot(p)
    
    
  } else {print(paste("not found", selection[i]))}
}

dev.off()
}

# selection <- tables_all_glmFit[[1]]$SYMBOL[1:10]
# make_boxplots("glm_PD", selection)
# selection <- tables_all_glmFit[[2]]$SYMBOL[1:10]
# make_boxplots("glm_CB", selection)
# selection <- tables_all_glmFit[[3]]$SYMBOL[1:10]
# make_boxplots("glm_diffChanges", selection)
# 
# selection <- tables_all[[1]]$SYMBOL[1:10]
# make_boxplots("QLF_PD", selection)
# selection <- tables_all[[2]]$SYMBOL[1:10]
# make_boxplots("QLF_CB", selection)
# selection <- tables_all[[3]]$SYMBOL[1:10]
# make_boxplots("QLF_diffChanges", selection)


selection <- c("CCL5", "IL12A", "IL12B", "CXCL9")
make_boxplots("IL12_CCL5", selection) #this data was used in prism to produce Figure 5B for CCL5.

genes_of_interest <- cpms_filtered[rownames(cpms_filtered) %in% c("CCL5", "IL12A", "IL12B", "CXCL9"),]
genes_of_interest_annot <-  data.frame("cpms_CCL5" = genes_of_interest["CCL5",], "benefit" = y$samples$benefit, "timepoint" = y$samples$timepoint, "patient" =  y$samples$Patient, "response" = factor(y$samples$Response, levels= c("PD", "SD", "PR", "CR")))
write.table(genes_of_interest_annot, file = "lists/CCL5_byPatient.txt", sep = "\t", row.names = FALSE, quote = FALSE)

############################################

#Build Design Matrix, make a batch correction for donor
grouped <- factor(paste0(y$samples$timepoint, "_", y$samples$benefit))
design_predictive <- model.matrix(~0+ grouped)#starting with a paired sample model
rownames(design_predictive) <- colnames(y)

colnames(design_predictive) <- gsub("y\\$samples\\$", "", colnames(design_predictive))
design_predictive

#Build GLM
y_pred <- estimateDisp(y, design_predictive)

plotBCV(y_pred, main="BCV Plot paired") #Plotting Biological Coefficient of Variation

fit_pred_glmFit <- glmFit(y_pred, design_predictive, prior.count = 8)
pcut=0.3
logfcut = 0


#define contrasts
mycontrasts2 <- list("CB-Prevs_PD-Pre" = c(0,0,1,-1), "CB-On_vs_PD-On" = c(1,-1,0,0))  #last one means genes that are more or less regulated in clinical benefit than in no benefit

lrt_all_glmFit_pred <- sapply(names(mycontrasts2),function(x) glmLRT(fit_pred_glmFit,contrast=mycontrasts2[[x]] ), simplify=FALSE, USE.NAMES = TRUE )
tables_all_glmFit_pred <- sapply(names(lrt_all_glmFit_pred), function(x)  topTags(lrt_all_glmFit_pred[[x]],adjust.method="BH",  n = Inf, sort.by = "PValue")$table, simplify=FALSE, USE.NAMES = TRUE)
sapply(tables_all_glmFit_pred, function(x) head(x, n=20),simplify=FALSE, USE.NAMES = TRUE)
tables_all_unsorted_glmFit_pred <- sapply(names(lrt_all_glmFit_pred), function(x)  topTags(lrt_all_glmFit_pred[[x]],adjust.method="BH", n = Inf, sort.by = "none")$table,simplify=FALSE, USE.NAMES = TRUE)

sapply(names(tables_all_glmFit_pred), function(x) write.table(tables_all_glmFit_pred[[x]],paste0("./lists/", x,"_dgelist_glm_full_predictive.txt"), sep ="\t", quote=FALSE, row.names = FALSE))




selection <- tables_all_glmFit_pred[[1]]$SYMBOL[1:10]
make_boxplots("QLF_pred_CB-Pre_vs_PD-Pre", selection)
selection <- tables_all_glmFit_pred[[2]]$SYMBOL[1:10]
make_boxplots("QLF_pred_CB-On_vs_On-Pre", selection)


selection <- c(
  "MLC1",
"UBAP1L",
"PCDH11X",
"CLEC17A",
"P2RX5",
"LRMP",
"C4BPA",
"ARSI",
"NELL2",
"CXCL13",
"TNR",
"FAM189A2",
"FAM159A",
"GPR88")
make_boxplots("Riaz_PreDiff", selection)



###################

cpms_tb_full <- tibble(data.frame(cpms_filtered))
cpms_tb_full$SYMBOL <- y$genes$SYMBOL
cpms_tb_full <- gather(cpms_tb_full, "sample", "cpm", -SYMBOL )
samples_annot <- samples %>% mutate("sample" = paste0(Patient,"_", timepoint, "_", ID))
cpms_tb_full <- left_join(cpms_tb_full, samples_annot, by = "sample")
cpms_tb_full$Response <- factor(cpms_tb_full$Response, levels = c("PD", "SD", "PR", "CR"))
cpms_tb_full$benefit <- ifelse(cpms_tb_full$Response == "PD", "PD", "CB")
cpms_tb_full$benefit <- factor(cpms_tb_full$benefit, levels = c("PD", "CB"))
cpms_tb <- cpms_tb_full %>% dplyr::select(SYMBOL, Patient, Response, timepoint, cpm, benefit) %>% spread(timepoint, cpm) %>% mutate(diff = On - Pre)


makeBenefitResponseBoxplots <- function(cpms_tb, gene_of_interest)
{
cpms_tb %>% dplyr::filter(SYMBOL == gene_of_interest) %>%
ggplot( aes(x=benefit, y=diff)) + geom_boxplot() + geom_point() + ylab("CPM On - Pre") + theme_bw() + ggtitle(gene_of_interest)

cpms_tb %>% dplyr::filter(SYMBOL == gene_of_interest) %>%
  ggplot( aes(x=Response, y=diff)) + geom_boxplot() + geom_point() + ylab("CPM On - Pre") + theme_bw() + ggtitle(gene_of_interest)
}

pdf("./plots/Delta_SelectedBoxplots_OnPre.pdf")

makeBenefitResponseBoxplots(cpms_tb, "CCL5")
makeBenefitResponseBoxplots(cpms_tb, "IL12A")
makeBenefitResponseBoxplots(cpms_tb, "CXCL9")
makeBenefitResponseBoxplots(cpms_tb, "CCL4")
makeBenefitResponseBoxplots(cpms_tb, "IL12B")
dev.off()



#name two genes of interest for which the delta will be plotted
genes_of_interest <- c("CCL5", "IL12A")

makeDiffPlot <- function(cpms_tb, genes_of_interest){
cpms_tb %>% dplyr::filter(SYMBOL %in% genes_of_interest ) %>% dplyr::select(Patient, Response, benefit, diff, SYMBOL) %>% spread(SYMBOL, diff) %>%
  ggplot(aes_string(x = genes_of_interest[1], y=genes_of_interest[2], col = "Response"))  + geom_point(aes_string( size = 3))  + 
    theme_bw() + geom_hline(yintercept=0)+ geom_vline(xintercept=0) + ggtitle("Delta On vs Pre Cpms")+
    scale_color_manual(values=c("#b03a2e", "#dc7633", "#3498db", "#229954")) +
    guides(color = guide_legend(override.aes = list(size = 3)))

}

#select one gene and for that gene, show pre levels on x and difference upon treatment on y
makePreDiffPlot <- function(cpms_tb, gene_of_interest)
{
cpms_tb %>% dplyr::filter(SYMBOL == gene_of_interest) %>%
  ggplot( aes(x=Pre, y=diff, col= Response)) + geom_point(size = 3 ) + ylab("CPM On - Pre")  +
  theme_bw() + geom_hline(yintercept=0)+ ggtitle("Delta On vs Pre Cpms")+
  scale_color_manual(values=c("#b03a2e", "#dc7633", "#3498db", "#229954")) +
  guides(color = guide_legend(override.aes = list(size = 3))) + ggtitle(paste(gene_of_interest), "Delta On-Pre VS Pre")
}

pdf("./plots/DeltaPlots_OnPre.pdf")
makeDiffPlot(cpms_tb, c("CCL5", "IL12A"))
makeDiffPlot(cpms_tb, c("CCL5", "CXCL9"))
makeDiffPlot(cpms_tb, c("CCL5", "PTPRC"))
makeDiffPlot(cpms_tb, c("CCL5", "CCL4"))
makeDiffPlot(cpms_tb, c("IL12A", "IL12B"))


makePreDiffPlot(cpms_tb, "CCL5")
makePreDiffPlot(cpms_tb, "IL12A")
makePreDiffPlot(cpms_tb, "CXCL9")
makePreDiffPlot(cpms_tb, "PTPRC")
makePreDiffPlot(cpms_tb, "CCL4")
makePreDiffPlot(cpms_tb, "IL12B")
dev.off()



############################################################################################################
############################################################################################################
#Signatures

## We excluded CCL5 from the T cell signature of Tirosh. Because we want to look at the correlation between CCL5 and the signature genes.
signatures_files <- list.files("./signatures/")

signatures <- sapply(signatures_files, function(x) read.table(paste0("./signatures/", x) , sep = "\t", stringsAsFactors = FALSE  )$V1, USE.NAMES = TRUE)
signatures

sapply(signatures, function(x) summary(x %in% genes$SYMBOL))
sapply(signatures, function(x) x[!x %in% genes$SYMBOL]) #which are not in the dataset

signatures <- sapply(signatures, function(x) x[x %in% genes$SYMBOL]) #only take those that are found
names(signatures) <- gsub(".txt","", names(signatures)) #give them good names
signatures

signatures_byPatient <- data.frame(matrix(ncol=length(signatures), nrow= ncol(cpms_filtered)))
for (j in 1:length(signatures)){ #for each signature
  for(i in 1:(ncol(cpms_filtered))) { # for each patient
    df <- data.frame(cpms = cpms_filtered[,i], genes = y$genes$SYMBOL)
    df <- df[order(df$cpms, decreasing = TRUE),]
    where <- which(df$genes %in% signatures[[j]])
    signatures_byPatient[i,j]  <- 1-mean(where)/nrow(cpms_filtered) #the score is defined as the mean rank of the NKsignature in the ranked cpm list. The score is inverted towards 1, because intuitively higher score then reflects more NK cells 
  }
}
colnames(signatures_byPatient) <- names(signatures)
rownames(signatures_byPatient) <- colnames(cpms_filtered)

signatures_byPatient.TopBottom <- signatures_byPatient
for(i in 1:ncol(signatures_byPatient)){
  signatures_byPatient.TopBottom[,i] <- ifelse(signatures_byPatient[,i] > median(signatures_byPatient[,i] ), "Top", "Bottom"  )
}
signatures_byPatient$sample <- colnames(cpms_filtered )
signatures_byPatient

colnames(signatures_byPatient.TopBottom) <- paste0(colnames(signatures_byPatient.TopBottom), "_TopBottom")
signatures_byPatient.TopBottom$sample <- colnames(cpms_filtered )
signatures_byPatient.TopBottom

cpms_tb_full <- left_join(cpms_tb_full, signatures_byPatient, by = "sample")
cpms_tb_full <- left_join(cpms_tb_full, signatures_byPatient.TopBottom, by = "sample")
cpms_tb_full



#####3

pdf("./plots/genes_of_signatures.pdf")
for(i in 1:length(names(signatures))){
  sign <-  names(signatures)[i]
  listofgenes <- signatures[[i]]
  p <- cpms_tb_full %>% filter(SYMBOL %in% listofgenes & timepoint == "Pre") %>%
    ggplot(aes(y=log(cpm), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("") + facet_wrap(~SYMBOL) + theme_bw() + ggtitle(paste("Repsonse Pre ", sign))
  print(p)
  p <- cpms_tb_full %>% filter(SYMBOL %in% listofgenes & timepoint == "On") %>%
    ggplot(aes(y=log(cpm), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("")  + theme_bw() + facet_wrap(~SYMBOL)+ ggtitle(paste("Repsonse On ", sign))
  print(p)
  }
dev.off()


#Join Signatures to the samples list
samples_tb <- tibble(y$samples)
samples_tb$sample <- rownames(y$samples)

samples_tb <- left_join(samples_tb,signatures_byPatient, by = "sample")
samples_tb <- left_join(samples_tb, signatures_byPatient.TopBottom, by = "sample")
samples_tb

table(samples_tb$Response, samples_tb$timepoint)

#write to files, to be used also in prism
# write.table(cpms_tb_full, file="./lists/expression_signatures_by_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(samples_tb, file="./lists/signatures_by_samples.txt", sep = "\t", quote = FALSE, row.names = FALSE)


pdf(file= "./plots/Signatures_Heatmaps.pdf")
###########################
#Heatmaps
#########################

#first investigate the overlap between the signatures
correlationmap <- cor(signatures_byPatient[grep(".txt",colnames(signatures_byPatient) )])
p <- Heatmap(correlationmap, column_title= "Correlation between signatures within patients")
print(p)


intersect_fraction_a <- function(a,b){
  length(intersect(a,b))/(length(a))
}

intersection_table <- data.frame(matrix(nrow=length(signatures), ncol = length(signatures)))
for(i in 1:length(signatures)){
  for(j in 1:length(signatures)){
    intersection_table[i,j] <- intersect_fraction_a(signatures[[i]], signatures[[j]])
  }
}
colnames(intersection_table) <- names(signatures)
rownames(intersection_table) <- names(signatures)

Heatmap(intersection_table, cluster_rows = FALSE, cluster_columns = FALSE)




#Heatmaps of Signatures by Sample
samples_tb$Patient <- factor(samples_tb$Patient)
samples_tb$Response <- factor(samples_tb$Response, levels = c("PD", "SD", "PR", "CR"))

topAn <- HeatmapAnnotation(
  df = data.frame(Timepoint =samples_tb$timepoint, Response = samples_tb$Response), annotation_name_gp = gpar(fontsize = 10),
  col= list(Response = structure(c("#b03a2e", "#dc7633", "#3498db", "#229954"), names= c("PD","SD","PR", "CR")),
            Timepoint = structure(c("#FFA100", "#FE4F00"), names= c("Pre", "On"))
  )
)
#Heatmap of signatures scaled
a <- samples_tb %>% dplyr::select(c(names(signatures)))  
b <- a[,1:length(signatures)] %>% scale() %>% t() 

ht <- Heatmap(b, col= colorRamp2(c(min(b), 0, max(b)), c("blue", "white", "red")) , column_title = "All Samples  -- Signatures scaled",  top_annotation = topAn, name = "signature\nscaled",  show_column_names = TRUE , row_names_side = "left", cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 9))
draw(ht, heatmap_legend_side = "right",  padding = unit.c(unit(2,"mm"), unit(17,"mm"),   unit(2,"mm"), unit(2,"mm")))

#same without CR
samples_noCR <- samples_tb %>% filter(Response != "CR")
a <- samples_noCR %>% dplyr::select(c(names(signatures)))  
b <- a[,1:length(signatures)] %>% scale() %>% t() 

topAn <- HeatmapAnnotation(
  df = data.frame(Timepoint =samples_noCR$timepoint, Response = samples_noCR$Response), annotation_name_gp = gpar(fontsize = 10),
  col= list(      Response = structure(c("#b03a2e", "#dc7633", "#007700"), names= c("PD","SD","PR")),
                  Timepoint = structure(c("#727EC3", "#011993"), names= c("Pre", "On")))
)

ht <- Heatmap(b, col= colorRamp2(c(min(b), 0, max(b)), c("blue", "white", "red")) , column_title = "All Samples  -- Signatures scaled",  top_annotation = topAn, name = "signature\nscaled",  show_column_names = TRUE , row_names_side = "left", cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 9))
draw(ht, heatmap_legend_side = "right",  padding = unit.c(unit(2,"mm"), unit(17,"mm"),   unit(2,"mm"), unit(2,"mm")))


#Heatmap of signatures not scaled
b <- a[,1:length(signatures)] %>% t() 

ht <- Heatmap(b, col= colorRamp2( seq(min(b), max(b), length.out= 256), viridis(256) ), column_title = "All Samples  -- Signatures",  top_annotation = topAn, name = "signature",  show_column_names = TRUE , row_names_side = "left", cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 9))
draw(ht, heatmap_legend_side = "right",  padding = unit.c(unit(2,"mm"), unit(17,"mm"),   unit(2,"mm"), unit(2,"mm")))

#only Pre Treatment Heatmaps
samples_onlyPre <- samples_tb %>% filter(timepoint == "Pre")
topAn <- HeatmapAnnotation(
  df = data.frame(Timepoint =samples_onlyPre$timepoint, Response = samples_onlyPre$Response), annotation_name_gp = gpar(fontsize = 10),
  col= list(      Response = structure(c("#b03a2e", "#dc7633", "#3498db", "#229954"), names= c("PD","SD","PR", "CR")),
            Timepoint = structure(c("#FFA100", "#FE4F00"), names= c("Pre", "On"))
  )
)

a <- samples_onlyPre %>% dplyr::select(c(names(signatures))) 


#raw
b <- a  %>% t() 
ht <- Heatmap(b, col= colorRamp2( seq(min(b), max(b), length.out= 256), viridis(256) ), column_title = "Pre Treatment  -- Signatures",  top_annotation = topAn, name = "signature",  show_column_names = TRUE , row_names_side = "left", cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 9))
draw(ht, heatmap_legend_side = "right",  padding = unit.c(unit(2,"mm"), unit(17,"mm"),   unit(2,"mm"), unit(2,"mm")))

#scaled
b <- a %>% scale() %>% t() 
ht <- Heatmap(b, col= colorRamp2(c(min(b), 0, max(b)), c("blue", "white", "red")), column_title = "Pre Treatment  -- Signatures scaled",  top_annotation = topAn, name = "signature\nscaled",  show_column_names = TRUE , row_names_side = "left", cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 9))
draw(ht, heatmap_legend_side = "right",  padding = unit.c(unit(2,"mm"), unit(17,"mm"),   unit(2,"mm"), unit(2,"mm")))

dev.off()



#format data into ggplot dormat
samples_melt <- samples_tb %>% dplyr::select(c(names(signatures), sample, timepoint, Response, Patient)) %>%
  gather( "signature", "value", -sample, -timepoint, -Response, -Patient )
samples_melt$value = as.numeric(samples_melt$value)
samples_melt

##################################
#Plots
#################################33
#Response vs Signatures
#Also calculate lm 
perform_signature_testing <- function(samples_melt, tp = "Pre") {
        pdf(paste0("./plots/Signatures_byResponse_", tp, ".pdf"))
      test_signatures <- data.frame(matrix(0, nrow=0, ncol=6)) 
      for(i in 1:length(names(signatures))){
        sign <-  names(signatures)[i]
        
        dat <- samples_melt %>% filter(timepoint == tp & signature == sign) 
        # dat <- dat %>% filter(SampleID != "X100001106.Pre") #Excluded one patient with very high values everywhere
        j <- jonckheere.test(dat$value, as.numeric(dat$Response)) #jonckheere test for overall trend
        
        res <- lm(value ~ Response, data = dat ) #linear model testing
        w <-  wilcox.test(value ~ Response, data = dat %>% filter(Response != "SD" & Response != "CR") )$p.value #wilcoxon only PD to PR
        test_signatures <- rbind(test_signatures, c(coef(summary(res))[3,], w, j$p.value ) ) #save all test results into a new table
        rownames(test_signatures)[i] <- as.character(sign)
        
       
        #plot NK score vs Response vs Timepoints 
        p <- ggplot(samples_tb %>% filter(timepoint == tp), aes_string(x="Response", y=sign)) + 
          geom_boxplot() +
          geom_point() +
          theme_bw()  +
          ggtitle(paste(tp,  "-- Response vs Signature"))  +
          xlab(paste("Response\n", sprintf("p value Wilcoxon PR vs PD %.4f", w), "\n", sprintf("p value Jonckheere Overall Trend %.4f", j$p.value)) )
        print(p)
      
      
      }
      dev.off()
      
      
      colnames(test_signatures) <- c("estimate", "std_error", "zvalue lm", "pvalueLM", "pvalueWilcoxon", "pvalueJonckheere")
      test_signatures[,"adj_pvalue_lm"] <- p.adjust(test_signatures$pvalueLM, method = "BH")
      test_signatures[,"adj_pvalue_wilcoxon"] <- p.adjust(test_signatures$pvalueWilcoxon, method = "BH")
      test_signatures[,"adj_pvalue_Jonckheere"] <- p.adjust(test_signatures$pvalueJonckheere, method = "BH")
      test_signatures$feature <- rownames(test_signatures)
      return(test_signatures)
}

test_signatures <- perform_signature_testing(samples_melt, "Pre")
write.table(test_signatures, file= "./lists/signatures_statistics_Pre.txt", sep = "\t", quote= FALSE, row.names = FALSE)
test_signagtures <- perform_signature_testing(samples_melt, "On")
write.table(test_signatures, file= "./lists/signatures_statistics_On.txt", sep = "\t", quote= FALSE, row.names = FALSE)


pdf("./plots/Signatures_More_Plots.pdf")

#### correlation among all signatures as heatmap
signatures_correlation <- cor(samples_tb[,names(signatures)])
rownames(signatures_correlation) <- gsub(".txt", "", rownames(signatures_correlation) )
colnames(signatures_correlation) <- gsub(".txt", "", colnames(signatures_correlation) )

ht <- Heatmap(signatures_correlation, color = colorRamp2( colorRamp2(c(min(signatures_correlation), 0, max(signatures_correlation)), c("blue", "white", "red"))), column_title = "Correlation among Signatures",  show_column_names = TRUE , row_names_side = "left", cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = 9))
draw(ht, heatmap_legend_side = "right",  padding = unit.c(unit(2,"mm"), unit(17,"mm"),   unit(2,"mm"), unit(2,"mm")))


#For both timepoints
for(i in 1:length(names(signatures))){
  signature <-  names(signatures)[i]
  p <- ggplot(samples_tb, aes_string(x="Response", y=signature)) + geom_boxplot() +geom_point() +theme_bw() + facet_wrap(~timepoint) #plot NK score vs Response vs Timepoints
  print(p)
  
  p <- ggplot(samples_tb, aes_string(x="Response", y=signature, fill = "timepoint")) + 
    geom_boxplot(position = position_dodge(0.8)) +
    geom_point(position = position_dodge(0.8)) +
    theme_bw()  #plot NK score vs Response vs Timepoints
  print(p)
}



# CCL5 vs Signatures or TopBottom
for(i in 1:length(names(signatures))){
  signature <-  names(signatures)[i]
  correlation <- cpms_tb_full %>% filter(SYMBOL == "CCL5") %>% dplyr::select(cpm, signature) %>%
    cor()
  dat <- cpms_tb_full %>% filter(SYMBOL == "CCL5") %>% dplyr::select(cpm, paste0(signature, "_TopBottom")) 
  colnames(dat) <-c("cpm", "TopBottom")
  ttest <- t.test(cpm ~ TopBottom , data = dat)
  p <- ggplot(cpms_tb_full %>% filter(SYMBOL == "CCL5"), aes_string(x="cpm", y=signature)) +geom_point() +theme_bw()  + xlab(paste("Expression CCL5\ncorrelation", sprintf("%.3f", correlation[1,2])) ) + ggtitle("CCl5 Expression vs Signature Score") 
  print(p)
  p <- ggplot(cpms_tb_full %>% filter(SYMBOL == "CCL5"), aes_string(y="cpm", x=paste0(signature, "_TopBottom")))  +geom_boxplot() + geom_point() + ggtitle("Expression CCL5 vs Top Bottom Classification") +
   xlab(paste("Signature Classification ", signature, "\nt-test", sprintf("%.2e", ttest$p.value )) ) 
  print(p)
}

#Save some of this data into separate files
selectedData  <- cpms_tb_full %>% filter(SYMBOL == "CCL5") %>% dplyr::select(cpm, "NK_Cursons", sample, Patient, timepoint, Response, benefit)
write.table(selectedData, file = "./lists/CCL5_vs_NK_Cursons.txt",  sep = "\t", quote = FALSE, row.names = FALSE)

selectedData  <- cpms_tb_full %>% filter(SYMBOL == "CCL5") %>% dplyr::select(cpm, "NK_Barry", sample, Patient, timepoint, Response, benefit)
write.table(selectedData, file = "./lists/CCL5_vs_NK_Barry.txt",  sep = "\t", quote = FALSE, row.names = FALSE)

selectedData  <- cpms_tb_full %>% filter(SYMBOL == "CCL5") %>% dplyr::select(cpm, "cDC1_Barry", sample, Patient, timepoint, Response, benefit)
write.table(selectedData, file = "./lists/CCL5_vs_cDC1_Barry.txt",  sep = "\t", quote = FALSE, row.names = FALSE)



#Show two selected comparisons between signatures
correlate_signature <- function(signature1 = "B_cell_Tirosh", signature2 = "B_cell_Tirosh"){
  correlation <- samples_tb %>% dplyr::select(signature1, signature2) %>% cor()
  p <- ggplot(samples_tb , aes_string(x=signature1, y=signature2)) +geom_point() +theme_bw() + ggtitle("Signature Correlations") + xlab(paste(signature1,"\ncorrelation", sprintf("%.3f", correlation[1,ifelse(signature1 == signature2, 1, 2)])) ) 
  print(p)
}

#plot every possible comparison
for(i in 1:length(signatures)){
  for(j in 1:length(signatures)){
    correlate_signature(names(signatures)[i], names(signatures)[j])
  }
}

dev.off()




































