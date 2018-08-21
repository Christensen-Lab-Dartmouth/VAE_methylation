######################
# Comparing significant CpGs from Random Forest to LaWAS results
#
# Author: Alexander Titus
# Created: 08/20/2018
# Updated: 08/20/2018
######################

#####################
# Set up the environment
#####################
     library(data.table)
     library(glmnet)
     library(ggplot2)
     library(ggpubr)
     library(gridExtra)


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
# Venn Diagram - VAE, EWAS, LASSO
#####################  
     #install.packages('VennDiagram')
     library(VennDiagram)
     
     ewas.file = 'ERposVERneg_EWAS.csv'
     ewas.dir = paste('results', ewas.file, sep = '/')
     ewas = data.frame(fread(ewas.dir)) 
     ewas = ewas[ewas$EWASqvalues < 0.001, ]
     ewas.cpg = ewas$CpG
     
     lasso.file = 'ERposVERneg_LASSO.csv'
     lasso.dir = paste('results', lasso.file, sep = '/')
     lasso = data.frame(fread(lasso.dir)) 
     lasso.cpg = lasso$Name
     
     vae.file = 'ERposVERneg_VAEcorr.csv'
     vae.dir = paste('results', vae.file, sep = '/')
     vae = data.frame(fread(vae.dir)) 
     vae.cpg = vae$CpG
     vae.cpg = unique(vae.cpg)
     
     png(filename = 'results/cpgVennDiagram.png', height = 1600, width = 1600, res = 300) 
     vennPlot <- venn.diagram(list(A = ewas.cpg,
                                   B = lasso.cpg,
                                   C = vae.cpg), 
                              NULL, lwd = 0.1, cat.cex = 2,cex = 3,
                              fill = c('deepskyblue','yellow', 'green'),
                              alpha = c(0.25, 0.25, 0.25),
                              fontface = 4, cat.fontface = 1, cat.dist = 0, margin = 0.05,
                              category.names = c('EWAS', 'LASSO', 'VAE'))
     grid.draw(vennPlot)
     dev.off()    
     
     
     
#####################
# Venn Diagram - VAE, PMD, HMD
#####################  
     vae.file = 'ERposVERneg_VAEcorr.csv'
     vae.dir = paste('results', vae.file, sep = '/')
     vae = data.frame(fread(vae.dir)) 
     vae.cpg = vae$CpG
     vae.cpg = unique(vae.cpg)
     
     pmd.file = 'hm450.comPMD.probes.csv'
     pmd.dir = paste('data', pmd.file, sep = '/')
     pmd = data.frame(fread(pmd.dir)) 
     pmd.cpg = pmd$CpG
     
     hmd.file = 'hm450.comHMD.probes.csv'
     hmd.dir = paste('data', hmd.file, sep = '/')
     hmd = data.frame(fread(hmd.dir)) 
     hmd.cpg = hmd$CpG

     png(filename = 'results/cpgVennDiagramSoloCpG.png', height = 1600, width = 1600, res = 300) 
     vennPlot2 <- venn.diagram(list(A = pmd.cpg,
                                   B = hmd.cpg,
                                   C = vae.cpg), 
                              NULL, lwd = 0.1, cat.cex = 2,cex = 3,
                              fill = c('deepskyblue','yellow', 'green'),
                              alpha = c(0.25, 0.25, 0.25),
                              fontface = 4, cat.fontface = 1, cat.dist = 0, margin = 0.05,
                              category.names = c('PMD', 'HMD', 'VAE'))
     grid.draw(vennPlot2)
     dev.off()    
     
     
#####################
# Volcano plot
##################### 
     ewas$Log10P = -log10(ewas$EWASpvalues)
     qcutHigh = min(ewas$Log10P[ewas$EWASqvalues < 0.001])
     qcutLow = min(ewas$Log10P[ewas$EWASqvalues < 0.01])
     
     cpgHigh = nrow(ewas[ewas$EWASqvalues < 0.001, ])
     cpgLow  = nrow(ewas[ewas$EWASqvalues >= 0.001 & ewas$EWASqvalues < 0.01, ])
     cpgUnder = nrow(ewas[ewas$EWASqvalues >= 0.01,])
     
     cpgHigh + cpgLow + cpgUnder == nrow(ewas)
     ewas$tum_cont_label = factor(ifelse(ewas$CpG %in% vae.cpg, 1, 0))
     ewas <- ewas[order(ewas$tum_cont_label,decreasing = F),]
     
     png(filename = 'Figures/VolcanoPlots/ColoredByTumContComparison/TumVsCont.png',width = 800,height = 800)
     p <- ggscatter(ewas, x = 'EWASbeta', y = 'Log10P', color = 'tum_cont_label',
                    alpha = 0.75,size = 1, title = 'ER-positive vs. ER-negative',
                    xlab = 'Beta Coefficient', ylab = expression('-log'[10]*' p value'))
     
     ggpar(p, legend = 'none',palette = c('black','red3'),
           font.x = 30,font.y = 30,font.main = 34,font.tickslab = 30) +
          # Add lines to indicate Q value cutoffs
          geom_hline(yintercept = qcutHigh, color = "red", linetype = "dashed") + 
          geom_hline(yintercept = qcutLow, color = "black", linetype = "dashed") #+ 
          # Add annotations of number of CpGs in each bin
          annotate(geom = 'text',x = -3.35,y = qcutHigh + 0.5,label = paste(format(cpgHigh,big.mark = ','),'CpGs'),size = 8) +
          annotate(geom = 'text',x = -3.35,y = qcutLow + 0.5,label = paste(format(cpgLow,big.mark = ','),'CpGs'),size = 8) +
          annotate(geom = 'text',x = -3.35,y = qcutLow - 0.5,label = paste(format(cpgUnder,big.mark = ','),'CpGs'),size =8)
     dev.off()
     