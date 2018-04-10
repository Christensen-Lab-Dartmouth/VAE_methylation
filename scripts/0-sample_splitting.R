#####################################
# Separate Samples by tumor and normal
# Author: Alexander J. Titus
# Date: March 30th, 2018
#####################################

#####################################
# Set up the environment
#####################################
     base.dir = 'C:/Users/atitus/github/VAE_analysis'
     setwd(base.dir)
     library(data.table)
     

#####################################
# Load and manipulate the data
#####################################   
     # Covariates - TCGA
     covs.file = 'data/Full_data_covs.csv'
     covs = data.frame(fread(covs.file))
     covs.tumor = covs[covs$SampleType == 'Primary Tumor', ]
     covs.normal = covs[covs$SampleType == 'Solid Tissue Normal', ]
     
     # Betas - TCGA
     beta.file = 'data/TCGA_BRCA/TCGA_BRCA_Betas.tsv'
     TCGA.betas = data.frame(fread(beta.file))
     rownames(TCGA.betas) = TCGA.betas$sample_id
     TCGA.betas = TCGA.betas[, 2:ncol(TCGA.betas)]
     colnames(TCGA.betas) = substr(colnames(TCGA.betas), 2, 100)
     TCGA.betas = TCGA.betas[order(rownames(TCGA.betas)), ]
     
     # Betas - Ringner
     beta.file = 'data/GSE75067_Ringner/GSE75067_Ringner_Betas.csv'
     Ringner.betas = data.frame(fread(beta.file))
     
     rowsRing = Ringner.betas$Basename
     Ringner.betas = Ringner.betas[, colnames(Ringner.betas) %in% rownames(TCGA.betas)]
     means = rowMeans(Ringner.betas, na.rm = T); means
     col.has.na <- apply(Ringner.betas, 2, function(x){any(is.na(x))}); col.has.na
     Ringner.betas = Ringner.betas[, !col.has.na]
     colsRing = colnames(Ringner.betas)
     
     # Betas - Fleischer
     beta.file = 'data/GSE84207_Fleischer/GSE84207_Fleischer_Betas.csv'
     Fleischer.betas = data.frame(fread(beta.file))
     rows = Fleischer.betas$V1
     cols = colnames(Fleischer.betas)
     cols[1] = 'CpG'
     Fleischer.betas = Fleischer.betas[Fleischer.betas$CpG %in% colnames(Ringner.betas), ]
     means = colMeans(Fleischer.betas[, 2:ncol(Fleischer.betas)], na.rm = T); means
     col.has.na <- apply(Fleischer.betas, 2, function(x){any(is.na(x))}); col.has.na
     Fleischer.betas = Fleischer.betas[, !col.has.na]
     colnames(Fleischer.betas) = cols
     rownames(Fleischer.betas) = Fleischer.betas$CpG
     
     
     # Sort data to make sure it all aligns correctly
     Ringner.betas = Ringner.betas[, colnames(Ringner.betas) %in% Fleischer.betas$CpG]
     TCGA.betas = TCGA.betas[rownames(TCGA.betas) %in% colnames(Ringner.betas), ]
     
     # Transpose and combine
     Ringner.betas = data.frame(t(Ringner.betas))
     colnames(Ringner.betas) = rowsRing
     
     Ringner.betas = Ringner.betas[order(rownames(Ringner.betas)), ]
     Fleischer.betas = Fleischer.betas[order(rownames(Fleischer.betas)), ]
     TCGA.betas = TCGA.betas[order(rownames(TCGA.betas)), ]
     
     all(rownames(TCGA.betas) == rownames(Ringner.betas))
     all(rownames(TCGA.betas) == rownames(Fleischer.betas))
     all(rownames(Fleischer.betas) == rownames(Ringner.betas))
     
     TCGA.betas$CpG = rownames(TCGA.betas)
     Ringner.betas$CpG = rownames(Ringner.betas)
     Fleischer.betas$CpG = rownames(Fleischer.betas)
     
     full.betas = merge(TCGA.betas, Ringner.betas, by = 'CpG')
     full.betas = merge(full.betas, Fleischer.betas, by = 'CpG')
     rownames(full.betas) = full.betas$CpG

     betas.tumor = full.betas[, colnames(full.betas) %in% covs.tumor$Basename]
     betas.normal = full.betas[, colnames(full.betas) %in% covs.normal$Basename]
     
     betas.tumor = cbind('CpG' = rownames(betas.tumor), betas.tumor)
     betas.normal = cbind('CpG' = rownames(betas.normal), betas.normal)

#####################################
# Write out the data
#####################################  
     fwrite(full.betas, 'data/Full_data_betas.csv', sep = ',')
     fwrite(betas.tumor, 'data/Full_data_Tumor_betas.csv', sep = ',')
     fwrite(betas.normal, 'data/Full_data_Normal_betas.csv', sep = ',')
     fwrite(betas.tumor, 'data/Full_data_Tumor_covs.csv', sep = ',')
     fwrite(betas.normal, 'data/Full_data_Normal_covs.csv', sep = ',')
     