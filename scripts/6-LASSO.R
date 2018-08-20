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
     library(glmnet)


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
# LASSO
#####################   
     temp = cbind(covs.updated$ERpos, betas)
     x=model.matrix(ERpos~.-1,data=temp) 
     y=Hitters$Salary
     