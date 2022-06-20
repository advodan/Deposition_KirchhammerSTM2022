# this is the analysis of the Nanostring Data from Algazi et al (Oncosec) of patient samples treated with IL12
# parts of Kirchhammer et al STM 2020 performed by Marcel P. Trenfny in 2020

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
library(ggfortify)
library(clinfun)
library(corrr)
myPalette <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"))
set.seed(123) #define randomness generation 


### Import Data files
path <- "/Volumes/CIMM$/CIMM/Projects_Human_Team/Marcel_Nicole/AdV5_IL12_Project/Nanostring_Algazi_IL12_deposition/" #set to file path

setwd(path)
data <- read.table("input/171002_NS_EB_Analysis_txt.txt", sep= "\t", header =  TRUE, stringsAsFactors=F)
head(data)

#not_used <- c( "X102.007.003.C2.NAU", "X102.001.016.C2.7U", "X102.001.020.C2.5U")
#data <- data[,!names(data) %in% not_used]

#extract expression and sample info
expression1 <- data[4:nrow(data),]
genes <- expression1$SampleName
expression <- data.frame(sapply(expression1[,2:ncol(expression1)], as.numeric))



expression$Symbol <- genes
head(expression)
samples <- data[1:3,]
head(samples)
rownames(samples) <- samples$SampleName
samples <- samples %>% 
  dplyr::select( -SampleName)
samples = as_tibble(t(samples), rownames = "SampleID")
samples
samples$Timepoint <- factor(samples$Timepoint, levels= c("Pre", "Post"))
samples$Response <- factor(samples$Response, levels= c("PD", "SD", "PR"))


expression_melt <- as_tibble(melt(expression, id.vars = "Symbol"))
colnames(expression_melt) <- c("Symbol", "SampleID", "Expression")
expression_melt
str(expression_melt)
expression_tb <- left_join(expression_melt, samples, by= c("SampleID"))
head(expression_tb)
expression_tb$Patient <- factor(expression_tb$Patient)
expression_tb$Timepoint <- factor(expression_tb$Timepoint, levels= c("Pre", "Post"))
expression_tb$Response <- factor(expression_tb$Response, levels= c("PD", "SD", "PR"))
expression_tb$Response <- factor(expression_tb$Response)
expression_tb$Symbol <- factor(expression_tb$Symbol)
expression_tb
table(expression_tb %>% dplyr::select(Patient, Timepoint))


#Plotting
#only pretreatemnt = Pre of selected genes
expression_tb %>% 
  ggplot(aes(y=log(Expression), x = Response)) + geom_violin() + geom_point() + ylab("Normalized Expression") + xlab("")  + theme_bw()

expression_tb %>% 
  ggplot(aes(x=log10(Expression+0.00001))) + geom_histogram(bins=100)


listofgenes <- c("CCL5", "CXCL9")
expression_tb %>% filter(Symbol %in% listofgenes & Timepoint == "Pre") %>%
ggplot(aes(y=log(Expression), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("") + facet_wrap(~Symbol) + theme_bw() + ggtitle("Pre")
expression_tb %>% filter(Symbol %in% listofgenes & Timepoint == "Post") %>%
ggplot(aes(y=log(Expression), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("") + facet_wrap(~Symbol) + theme_bw() + ggtitle("Post")
#NCR1
listofgenes <- c("NCR1")
expression_tb %>% filter(Symbol %in% listofgenes) %>%
  ggplot(aes(y=log(Expression), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("") + facet_wrap(~Timepoint) + theme_bw() + ggtitle(paste("Pre"), listofgenes)
#CD45
listofgenes <- c("PTPRC")
expression_tb %>% filter(Symbol %in% listofgenes) %>%
  ggplot(aes(y=log(Expression), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("") + facet_wrap(~Timepoint) + theme_bw() + ggtitle(paste("Pre"), listofgenes)


#
listofgenes <- c("CCL5", "CXCL9")
expression_tb %>% filter(Symbol %in% listofgenes ) %>%
  ggplot(aes(y=log(Expression), x = Timepoint)) + geom_boxplot() + geom_point(size= 2) + geom_line(aes(group = Patient))   + ylab("Normalized Expression") + xlab("") + facet_wrap(~Symbol) + theme_bw() + ggtitle("Pre")

# 

##################

results <- data.frame(matrix(0, nrow=0, ncol=4)) #create new empty dataframe
for (i in 1:length(unique(expression_tb$Symbol))){
  gene <-  unique(expression_tb$Symbol)[i]
  expression_tb_filtered <- expression_tb %>% filter(Timepoint == "Pre", Symbol == gene)
  res <- lm(log(Expression) ~ Response  ,data = expression_tb_filtered  )
  results <- rbind(results, c(coef(summary(res))[3,]))  #thrid row means "PD" vs "PR"
}


results$Symbol <- unique(expression_tb$Symbol)
colnames(results) <- c("estimate", "std_error", "tvalue", "pvalue", "Symbol")


results[,"adj_pvalue"] <- p.adjust(results$pvalue + 1E-300, method = "BH") #add a minimal p-value of 1E-300 because otherwise they are 0 and can't be log transformed
results[,"log2_fc"] <- log2(exp(results$estimate))
results[,"abs_log2_fc"] <- abs(results$log2_fc)
results <- as_tibble(results)
results
summary(results$pvalue)
summary(results$adj_pvalue)


############################################################################################################
############################################################################################################
#Signatures

## We excluded CCL5 from the T cell signature of Tirosh. Because we want to look at the correlation between CCL5 and the signature genes.
signatures_files <- list.files("./signatures_NanostringPanCancerImmunology//")

signatures <- sapply(signatures_files, function(x) read.table(paste0("./signatures_NanostringPanCancerImmunology/", x) , sep = "\t", stringsAsFactors = FALSE, header = FALSE  )$V1, USE.NAMES = TRUE)
signatures

sapply(signatures, function(x) summary(x %in% genes))
#of NK1 signature found 11 of 40 and NK2: 20 of 88
sapply(signatures, function(x) x[!x %in% genes])
sapply(signatures, function(x) x[x %in% genes])

signatures <- sapply(signatures, function(x) x[x %in% genes])
names(signatures) <- gsub(".txt","", names(signatures))
signatures

signatures_byPatient <- data.frame(matrix(ncol=length(signatures), nrow= ncol(expression)-1))
for (j in 1:length(signatures)){ #for each signature
  for(i in 1:(ncol(expression)-1)) { # for each patient
   df <- data.frame(cpms = expression[,i], genes = genes)
   df <- df[order(df$cpms, decreasing = TRUE),]
   where <- which(df$genes %in% signatures[[j]])
   signatures_byPatient[i,j]  <- 1-mean(where)/nrow(expression) #the score is defined as the mean rank of the NKsignature in the ranked cpm list. The score is inverted towards 1, because intuitively higher score then reflects more NK cells 
  }
}
colnames(signatures_byPatient) <- names(signatures)
rownames(signatures_byPatient) <- colnames(expression %>% dplyr::select(-Symbol))

signatures_byPatient.TopBottom <- signatures_byPatient
for(i in 1:ncol(signatures_byPatient)){
  signatures_byPatient.TopBottom[,i] <- ifelse(signatures_byPatient[,i] > median(signatures_byPatient[,i] ), "Top", "Bottom"  )
}

signatures_byPatient$SampleID <- colnames(expression %>% dplyr::select(-Symbol))

colnames(signatures_byPatient.TopBottom) <- paste0(colnames(signatures_byPatient.TopBottom), "_TopBottom")
signatures_byPatient.TopBottom$SampleID <- colnames(expression %>% dplyr::select(-Symbol))
signatures_byPatient.TopBottom

expression_tb <- left_join(expression_tb, signatures_byPatient, by = "SampleID")
expression_tb <- left_join(expression_tb, signatures_byPatient.TopBottom, by = "SampleID")
expression_tb

pdf("./plots/genes_of_signatures.pdf")
for(i in 1:length(names(signatures))){
  sign <-  names(signatures)[i]
listofgenes <- signatures[[i]]
 p <- expression_tb %>% filter(Symbol %in% listofgenes & Timepoint == "Pre") %>%
    ggplot(aes(y=log(Expression), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("") + facet_wrap(~Symbol) + theme_bw() + ggtitle(paste("Pre ", sign))
 print(p)
  p <- expression_tb %>% filter(Symbol %in% listofgenes & Timepoint == "Post") %>%
    ggplot(aes(y=log(Expression), x = Response)) + geom_boxplot() + geom_point() + ylab("log Expression") + xlab("")  + theme_bw() + facet_wrap(~Symbol)+ ggtitle(paste("Post ", sign))
print(p)
}
dev.off()

#Join Signatures to the samples list
samples_original <- samples
samples <- left_join(samples, signatures_byPatient, by = "SampleID")
samples <- left_join(samples, signatures_byPatient.TopBottom, by = "SampleID")
samples

samples_melt <- gather(samples, "signature", "value", -SampleID, -Timepoint, -Response, -Patient )
samples_melt$value = as.numeric(samples_melt$value)
samples_melt


#write to files, to be used in prism to generate plots for figures. one patient 10001106 was excluded later as a strong outlier
write.table(samples, file="./lists/signatures_by_samples.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_tb, file="./lists/expression_signatures_by_gene.txt", sep = "\t", quote = FALSE, row.names = FALSE)



#find which genes correlate the most with NK cells
contrast.df <- data.frame(matrix(1, nrow=0, ncol=4))
colnames(contrast.df) <-  c("estimate", "std_error", "zvalue", "pvalue")
for(i in 1:length(genes)){
  res <- lm(log(Expression) ~ NK_cells_cytotoxic,  data= expression_tb %>% filter(Symbol == genes[i])) 
  contrast.df <- rbind(contrast.df, c(coef(summary(res))[2,])) #save data to dataframe of this contrast, only save column 2 which is defined as the contrast columns above
  rownames(contrast.df)[i] <- as.character(genes[i])
}
colnames(contrast.df) <-  c("estimate", "std_error", "zvalue", "pvalue")
contrast.df$adj_pvalue <- p.adjust(contrast.df$pvalue, method = "BH")

contrast.df <- contrast.df[order(contrast.df$adj_pvalue, decreasing = FALSE),]
contrast.df$symbol = rownames(contrast.df)

write.table(contrast.df, file = "./lists/genes_correlating_NK_cytotoxic_Signature.txt", sep = "\t", quote = FALSE, row.names = FALSE)


