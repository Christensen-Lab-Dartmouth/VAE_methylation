######################
# Comparing significant CpGs from EWAS to LaWAS results
#
# Author: Alexander Titus
# Created: 08/20/2018
# Updated: 08/20/2018
######################

#####################
# Set up the environment
#####################
     require(data.table)
     require(limma)
     library(nlme)
     library(glmnet)
     library(qvalue)
     library(ggpubr)

######################
## Set WD
# Change this to your base WD, and all other code is relative to the folder structure
     base.dir = 'C:/Users/atitus/github/VAE_methylation'
     setwd(base.dir)
     
     # Covariate data
     covs.file = 'Full_data_covs.csv'
     covs.dir = paste('data', covs.file, sep = '/')
     covs = data.frame(fread(covs.dir))
     
     covs.updated = covs
     covs.updated = covs.updated[covs.updated$SampleType != 'Metastatic', ]
     covs.updated$BasalVother = ifelse(covs.updated$PAM50 == "Basal", 1, 0)
     covs.updated$NormalVother = ifelse(covs.updated$PAM50 == "Normal", 1, 0)
     covs.updated$Her2Vother = ifelse(covs.updated$PAM50 == "Her2", 1, 0)
     covs.updated$LumAVother = ifelse(covs.updated$PAM50 == "LumA", 1, 0)
     covs.updated$LumBVother = ifelse(covs.updated$PAM50 == "LumB", 1, 0)
     covs.updated$LumVother = ifelse(covs.updated$PAM50 == "LumA" | 
                                          covs.updated$PAM50 == "LumB", 1, 0)
     
     covs.updated$sample.typeInt = ifelse(covs.updated$SampleType == 'Solid Tissue Normal', 0, 1)
     covs.updated$ERpos = ifelse(covs.updated$ER == "Positive", 1, 
                                 ifelse(covs.updated$ER == "Negative", 0, NA))
     
     
     # Methylation data
     beta.file = 'BreastCancerMethylation_top100kMAD_cpg.csv'
     beta.dir = paste('data/raw', beta.file, sep = '/')
     betas = data.frame(fread(beta.dir)) # on my computer takes ~8min
     rownames(betas) = betas[,1]
     betas = betas[,2:ncol(betas)]
     
     betas = betas[rownames(betas) %in% covs.updated$Basename, ]
     betas = betas[order(rownames(betas), decreasing=T), ]
     
     covs.updated = covs.updated[covs.updated$Basename %in% rownames(betas), ]
     covs.updated = covs.updated[order(covs.updated$Basename, decreasing=T), ]
     
     ## Check sample concordance
     all(covs.updated$Basename == rownames(betas))
     
     
     # Annotation file with breast specific enhancer information
     anno.file = 'Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv'
     anno.dir = paste('data', anno.file, sep = '/')
     anno = data.frame(fread(anno.dir))
     rownames(anno) = anno[, 1]
     
     
#####################
# EWAS
#####################   
     ######################
     ## EWAS ERpos v ERneg
     # Create our model matrix
     BRCA.covsSub = covs.updated[!is.na(covs.updated$Age), ]
     BRCA.covsSub = BRCA.covsSub[!is.na(BRCA.covsSub$ERpos), ]
     
     betas.sub = betas[rownames(betas) %in% BRCA.covsSub$Basename,]
     
     XX <- model.matrix(~ERpos + Age, data = BRCA.covsSub)
     rownames(XX) <- BRCA.covsSub$Basename
     XX = XX[order(rownames(XX), decreasing=T), ]
     
     betas.sub = betas.sub[order(rownames(betas.sub), decreasing=T), ]
     
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
     write.csv(results.anno, file = 'results/ERposVERneg_EWAS.csv')
     
     # print out summary of the adjusted and unadjusted results
     summary(q.values)
     