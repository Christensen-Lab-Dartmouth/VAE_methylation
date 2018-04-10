######################
# EWAS for TCGA BRCA
#
# Author: Alexander Titus
# Created: 03/21/2018
# Updated: 03/21/2018
######################

## Setup working environment
library(data.table) # for fread
library(nlme)
library(glmnet)
library(data.table)
library(limma)
library(qvalue)
library(ggpubr)

######################
## Set WD
# Change this to your base WD, and all other code is relative to the folder structure
base.dir = 'C:/Users/atitus/github/VAE_analysis'
setwd(base.dir)


######################
## Read in data
# Methylation data
beta.file = 'TCGA_BRCA_top300kMAD_cpg.tsv'
beta.dir = paste('data', beta.file, sep = '/')
betas = data.frame(fread(beta.dir)) # on my computer takes ~8min
rownames(betas) = betas[,1]
betas = betas[,2:ncol(betas)]

# Covariate data
covs.file = 'BRCAtarget_covariates.csv'
covs.dir = paste('data', covs.file, sep = '/')
covs = data.frame(fread(covs.dir))

# Annotation file with breast specific enhancer information
anno.file = 'Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv'
anno.dir = paste('data', anno.file, sep = '/')
anno = data.frame(fread(anno.dir))
rownames(anno) = anno[, 1]
anno = anno[, 2:ncol(anno)]

######################
## Clean the data
covs.updated = covs
covs.updated = covs.updated[covs.updated$sample.type != 'Metastatic', ]
covs.updated$BasalVother = ifelse(covs.updated$PAM50.RNAseq == "Basal", 1, 0)
covs.updated$NormalVother = ifelse(covs.updated$PAM50.RNAseq == "Normal", 1, 0)
covs.updated$Her2Vother = ifelse(covs.updated$PAM50.RNAseq == "Her2", 1, 0)
covs.updated$LumAVother = ifelse(covs.updated$PAM50.RNAseq == "LumA", 1, 0)
covs.updated$LumBVother = ifelse(covs.updated$PAM50.RNAseq == "LumB", 1, 0)
covs.updated$LumVother = ifelse(covs.updated$PAM50.RNAseq == "LumA" | 
                                     covs.updated$PAM50.RNAseq == "LumB", 1, 0)

covs.updated$sample.typeInt = ifelse(covs.updated$sample.type == 'Solid Tissue Normal', 0, 1)


######################
## EWAS Tumor v. normal
# Create our model matrix
BRCA.covsSub = covs.updated[!is.na(covs.updated$age.Dx), ]
betas.sub = betas[rownames(betas) %in% BRCA.covsSub$Basename,]

XX <- model.matrix(~sample.typeInt + age.Dx, data = BRCA.covsSub)
rownames(XX) <- BRCA.covsSub$Basename
XX = XX[order(rownames(XX), decreasing=T), ]

betas.sub = betas.sub[order(rownames(betas.sub), decreasing=T), ]

#vae = vae[, 1:ncol(vae)-1]
betas.mat <- data.matrix(betas.sub)

# M-values
betas_dmgrM = ifelse(betas.mat >= 1, 1-1E-6, 
                     ifelse(betas.mat <= 0, 1E-6, betas.mat))
betas_dmgrM <- log(betas_dmgrM)-log(1-betas_dmgrM)
betas_dmgrM = t(betas_dmgrM)
betas_dmgrM2 <- betas_dmgrM[, colnames(betas_dmgrM) %in% rownames(XX)]
all(colnames(betas_dmgrM2) == rownames(XX))

# Limma Models - correlation = ar1
subjectList <- BRCA.covsSub$Basename[order(BRCA.covsSub$Basename,decreasing = T)]
corStructure <- corAR1(form = ~ 1 | subjectList)
lf_Null <- eBayes(lmFit(betas_dmgrM2, XX, block = subjectList, correlation = corStructure))

# Analysis of q-values
q.values <- qvalue(lf_Null$p.value[,2], lambda = seq(0.0, 0.55, 0.05))
results <- cbind(q.values$pvalues, q.values$qvalues, lf_Null$coefficients[, 2])

# Name the rows and columns for each model
colnames(results) <- colnames(results) <- c("EWASpvalues", "EWASqvalues", "EWASbeta")
results = data.frame(results)
results = results[order(results$EWASqvalues, decreasing = T), ]
results$CpG = rownames(results)

# Add EWAS results to annotation information with breast specific enhancers
anno.sub = anno[rownames(anno) %in% rownames(results), ]
results.anno = merge(results, anno.sub, by.x = 'CpG', by.y = 'Name')
rownames(results.anno) = results.anno$CpG

# Write file
write.csv(results.anno, file = 'results/TumorvNorm_EWAS.csv')

# print out summary of the adjusted and unadjusted results
summary(q.values)
