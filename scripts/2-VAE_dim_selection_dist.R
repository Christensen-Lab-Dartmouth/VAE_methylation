######################
# Dimension selection for TCGA BRCA
#
# Author: Alexander Titus
# Created: 03/21/2018
# Updated: 03/21/2018
######################

### Setup the environment we need ###
source("https://bioconductor.org/biocLite.R")
# biocLite("qvalue")
# install.packages('glmnet')
# install.packages('data.table')
require(glmnet)
require(data.table)
library(limma)
library(qvalue)
library(robustbase)

######################
## Set WD
# Change this to your base WD, and all other code is relative to the folder structure
base.dir = 'C:/Users/atitus/github/VAE_methylation'
setwd(base.dir)


######################
## Read in data
vae.file = 'encoded_methyl_onehidden_warmup_batchnorm_100k_100.tsv'
vae.dir = paste('results/VAE_training', vae.file, sep = '/')
vae = data.frame(fread(vae.dir))

colnames(vae) = as.character(vae[1, ])
vae = vae[2:nrow(vae), ]
rownames(vae) = vae[, 1]
meltVAE = melt(vae)
colnames(meltVAE) = c('Basename', 'Dimension', 'Activation')

#rownames(vae) = vae[, 1]
#vae = vae[, 2:ncol(vae)]


# Covariate data
covs.file = 'Full_data_covs.csv'
covs.dir = paste('data', covs.file, sep = '/')
covs = data.frame(fread(covs.dir))

######################
## Clean the data 
covs.updated = covs
covs.updated = covs.updated[covs.updated$SampleType != 'Metastatic', ]
covs.updated$TripleNeg = ifelse(covs.updated$ER == 'Negative' &
                                     covs.updated$PR == 'Negative' &
                                     covs.updated$HER2 == 'Negative', 'Triple Negative', 'Not Triple Negative')
merged.data = merge(covs.updated, vae, by.x = 'Basename', by.y = 'CpG')


######################
# Visualize VAE distirbution
#####################
     # ggplot theme
     #####################   
     theme_Publication <- function(base_size=14, base_family="helvetica") {
          library(grid)
          library(ggthemes)
          (theme_foundation(base_size=base_size, base_family=base_family)
               + theme(plot.title = element_text(face = "bold",
                                                 size = rel(1.2), hjust = 0.5),
                       text = element_text(),
                       panel.background = element_rect(colour = NA),
                       plot.background = element_rect(colour = NA),
                       panel.border = element_rect(colour = NA),
                       axis.title = element_text(face = "bold",size = rel(1)),
                       axis.title.y = element_text(angle=90,vjust =2),
                       axis.title.x = element_text(vjust = -0.2),
                       axis.text = element_text(), 
                       axis.line = element_line(colour="black"),
                       axis.ticks = element_line(),
                       panel.grid.major = element_line(colour="#f0f0f0"),
                       panel.grid.minor = element_blank(),
                       plot.margin=unit(c(10,5,5,5),"mm"),
                       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                       strip.text = element_text(face="bold")
               ))
          
     }

     meltVAE = meltVAE[meltVAE$Basename %in% covs.updated$Basename, ]
     
     # ER pos v. neg
     meltVAE$ER = rep('ER', nrow(meltVAE))
     meltVAE$PR = rep('PR', nrow(meltVAE))
     meltVAE$HER2 = rep('HER2', nrow(meltVAE))
     meltVAE$TripleNeg = rep('TripleNeg', nrow(meltVAE))
     meltVAE$PAM50 = rep('PAM50', nrow(meltVAE))
     
     for(i in 1:nrow(meltVAE)){
          temp = covs.updated[covs.updated$Basename == meltVAE[i, ]$Basename, ]
          meltVAE[i, ]$ER = temp$ER
          meltVAE[i, ]$PR = temp$PR
          meltVAE[i, ]$HER2 = temp$HER2
          meltVAE[i, ]$TripleNeg = temp$TripleNeg
          meltVAE[i, ]$PAM50 = temp$PAM50
     }
     
     vaeMedians = colMedians(as.matrix(vae[, 2:ncol(vae)]))
     
     # ER status
     temp = meltVAE[meltVAE$ER != '', ]
     p = ggplot(temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum),  
                             y = Activation, fill = factor(Dimension)))
     
     file.name = 'vaeDim_ERstatus_activation_distribution_100D.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2400, height = 1000, res = 100)
          p + geom_boxplot() + facet_grid(ER ~ .) + # two panel
               labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
               theme_Publication() + theme(legend.position="none")
     dev.off()
     
     # PR status
     temp = meltVAE[meltVAE$PR != '', ]
     temp = temp[!is.na(temp$PR), ]
     p = ggplot(temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum),  
                          y = Activation, fill = factor(Dimension)))
     
     file.name = 'vaeDim_PRstatus_activation_distribution_100D.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2400, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(PR ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme_Publication() + theme(legend.position="none")
     dev.off()
     
     # HER2 status
     temp = meltVAE[meltVAE$HER2 != '', ]
     p = ggplot(temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum),  
                          y = Activation, fill = factor(Dimension)))
     
     file.name = 'vaeDim_HER2status_activation_distribution_100D.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2400, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(HER2 ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme_Publication() + theme(legend.position="none")
     dev.off()
     
     # TripleNeg status
     temp = meltVAE[meltVAE$TripleNeg != '', ]
     p = ggplot(temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum),  
                          y = Activation, fill = factor(Dimension)))
     
     file.name = 'vaeDim_TripleNegstatus_activation_distribution_100D.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2400, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(TripleNeg ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme_Publication() + theme(legend.position="none")
     dev.off()
     
     # PAM50
     temp = meltVAE[meltVAE$PAM50 != '', ]
     temp = temp[!is.na(temp$PAM50), ]
     temp = temp[temp$PAM50 != 'Normal', ]
     p = ggplot(temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum),  
                          y = Activation, fill = factor(Dimension)))
     
     file.name = 'vaeDim_PAM50_activation_distribution_100D.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2000, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(PAM50 ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme(legend.position="none")
     p$labels$fill = 'Sample type'
     dev.off()
     

# Separate ER status
     ERpositive = merged.data[merged.data$ER == 'Positive', ]
     ERnegative = merged.data[merged.data$ER == 'Negative', ]
     
     
# Median dimension activations
     ERPVaeMedians = colMedians(as.matrix(ERpositive[, (ncol(ERpositive)-100+1):ncol(ERpositive)]))
     ERNVaeMedians = colMedians(as.matrix(ERnegative[, (ncol(ERnegative)-100+1):ncol(ERnegative)]))
     medians = data.frame(cbind(vaeMedians, ERPVaeMedians, ERNVaeMedians))
     
     medians = medians[order(medians[,2], medians[,3], medians[,1]), ]

     
# Identify the dimensions reprentative of tumor and normal, or both
     ERPMeds = medians[medians$ERPVaeMedians > 0 & 
                             medians$ERNVaeMedians == 0, ]
     write.csv(normMeds, file = 'results/medianDimActivation_ERpositive.csv')
     
     ERNMeds = medians[medians$ERPVaeMedians == 0 & 
                             medians$ERNVaeMedians > 0, ]
     write.csv(tumorMeds, file = 'results/medianDimActivation_ERnegative.csv')
     
     ERMeds = medians[medians$ERPVaeMedians > 0 & 
                                medians$ERNVaeMedians > 0, ]
     write.csv(normTumMeds, file = 'results/medianDimActivation_ERstatus.csv')
     
     
# Visualize VAE distribution of dimensions reprentative of ER status
     meltVAE.temp = meltVAE[meltVAE$Dimension %in% rownames(ERMeds), ]
     meltVAE.temp = meltVAE.temp[meltVAE.temp$ER != '', ]
     
     # ER status
     p = ggplot(meltVAE.temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum), 
                             y = Activation, fill = factor(Dimension))) # ordered by median

     file.name = 'vaeDim_activation_distribution_ERstatus.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2000, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(ER ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme_Publication() + theme(legend.position="none")
     dev.off()
     
     
# Visualize VAE distribution of dimensions reprentative of ER+ & ER-
     # ER+
     meltVAE.temp = meltVAE[meltVAE$Dimension %in% rownames(ERPMeds), ]
     meltVAE.temp = meltVAE.temp[meltVAE.temp$ER != '', ]
     
     p = ggplot(meltVAE.temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum), 
                                  y = Activation, fill = factor(Dimension))) # ordered by median
     
     file.name = 'vaeDim_activation_distribution_ERpositive.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2000, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(ER ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme_Publication() + theme(legend.position="none")
     dev.off()
     
     # ER-
     meltVAE.temp = meltVAE[meltVAE$Dimension %in% rownames(ERNMeds), ]
     meltVAE.temp = meltVAE.temp[meltVAE.temp$ER != '', ]
     
     p = ggplot(meltVAE.temp, aes(x = reorder(factor(Dimension), Activation, FUN = sum), 
                                  y = Activation, fill = factor(Dimension))) # ordered by median
     
     file.name = 'vaeDim_activation_distribution_ERnegative.png'
     file.dir = paste('results', file.name, sep = '/')
     png(file.dir, width = 2000, height = 1000, res = 100)
     p + geom_boxplot() + facet_grid(ER ~ .) + # two panel
          labs(x = 'VAE latent dimension', y = 'VAE latent dimension activation') +
          theme_Publication() + theme(legend.position="none")
     dev.off()
     
     
     