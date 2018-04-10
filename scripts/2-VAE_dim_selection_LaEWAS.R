######################
# LaEWAS for TCGA BRCA
#
# Author: Alexander Titus
# Created: 10/15/2017
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

#dir = '/Users/alexandertitus/Documents/github/DNAm_data_generation/code'
dir = 'C:/Users/atitus/Documents/github/VAE_analysis/results'
setwd(dir)


# Load all the data and make variables
### Load our data to work with ###
vae.file = 'encoded_methyl_onehidden_warmup_batchnorm_300K-100.tsv'
BRCA.covFile = '../BRCAtarget_covariates.csv'
tSNE.file = '../results/tSNE/vae_tsne_out_300K-100_3d.tsv'

tSNE_features = data.frame(fread(tSNE.file), row.names=1)
colnames(tSNE_features) = tSNE_features[1, ]
tSNE_features = tSNE_features[2:nrow(tSNE_features), ]
tSNE_features$Basename = rownames(tSNE_features)

vae_features = data.frame(fread(vae.file), row.names=1)
colnames(vae_features) = vae_features[1, ]
vae_features = vae_features[2:nrow(vae_features), ]
vae_features$Basename = rownames(vae_features)

BRCA.covs = data.frame(fread(BRCA.covFile), row.names=1)
BRCA.covs$sample.typeInt = ifelse(BRCA.covs$sample.type == 'Solid Tissue Normal', 0, 1)

full.data = merge(vae_features, BRCA.covs, by='Basename')
full.data = merge(full.data, tSNE_features, by='Basename')
full.data$BasalVother = ifelse(full.data$PAM50.RNAseq == "Basal", 1, 0)
full.data$NormalVother = ifelse(full.data$PAM50.RNAseq == "Normal", 1, 0)
full.data$Her2Vother = ifelse(full.data$PAM50.RNAseq == "Her2", 1, 0)
full.data$LumAVother = ifelse(full.data$PAM50.RNAseq == "LumA", 1, 0)
full.data$LumBVother = ifelse(full.data$PAM50.RNAseq == "LumB", 1, 0)
full.data$LumVother = ifelse(full.data$PAM50.RNAseq == "LumA" | 
                                  full.data$PAM50.RNAseq == "LumB", 1, 0)


## LaWAS
# Create our model matrix
BRCA.covsSub = BRCA.covs[!is.na(BRCA.covs$age.Dx), ]
vae_featuresSub = vae_features[rownames(vae_features) %in% BRCA.covsSub$Basename,]

XX <- model.matrix(~sample.typeInt + age.Dx, data = BRCA.covsSub)
rownames(XX) <- BRCA.covsSub$Basename
XX = XX[order(rownames(XX), decreasing=T), ]

vae = vae_featuresSub[order(rownames(vae_featuresSub), decreasing=T), ]
vae = vae[, 1:ncol(vae)-1]
vae_mat <- data.matrix(vae)

## M-values
betas_dmgrM = vae_mat # ifelse(vae_mat>=1, 1-1E-6, ifelse(vae_mat<=0, 1E-6, vae_mat))
#betas_dmgrM <- log(betas_dmgrM)-log(1-betas_dmgrM)
betas_dmgrM = t(betas_dmgrM)
betas_dmgrM2 <- betas_dmgrM[, colnames(betas_dmgrM) %in% rownames(XX)]
all(colnames(betas_dmgrM2) == rownames(XX))

# Limma models
lf_Null <- eBayes(lmFit(betas_dmgrM2, XX))


# Limma Models - correlation = ar1
library(nlme)
subjectList <- BRCA.covsSub$Basename[order(BRCA.covsSub$Basename,decreasing = T)]
corStructure <- corAR1(form = ~ 1 | subjectList)
lf_Null <- eBayes(lmFit(betas_dmgrM2, XX, block = subjectList, correlation = corStructure))

# Analysis of q-values
q.values <- qvalue(lf_Null$p.value[,2], lambda = seq(0.0, 0.55, 0.05))
results <- cbind(q.values$pvalues, q.values$qvalues, lf_Null$coefficients[, 2])

# Name the rows and columns for each model
colnames(results) <- colnames(results) <- c("pvalues", "qvalues", "beta")
results = data.frame(results)

# print out summary of the adjusted and unadjusted results
summary(q.values)



#### Generate the theme of each ggplot ####
results$negLog10P = -log(results$pvalues)
results$Nodes = rownames(results)
library(ggpubr)

png('Nodes_LaEWAS.png', width = 2000, height = 2000, res = 300)

ggscatter(results, x = "beta", y = "negLog10P",
              color = "black", 
              xlab = 'Beta', 
              ylab = '-log10(P-value)',
              label = 'Nodes',
              repel = T)

dev.off()

library(ggplot2)
overall_theme <-  theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(), 
                        panel.background = element_blank(), 
                        panel.border = element_blank(),
                        axis.line = element_line(size = rel(4), color = "black"), 
                        axis.text = element_text(size = rel(6), color = "black"),
                        axis.ticks = element_line(size = rel(5), color = "black"),
                        axis.ticks.length = unit(0.5, "cm"),
                        axis.title = element_text(size = rel(6.8)),
                        axis.title.y = element_text(vjust = 4.5),
                        axis.title.x = element_text(vjust = -4.5),
                        plot.title = element_text(size = rel(7.5)),
                        legend.key.size = unit(2, "cm"),
                        legend.text = element_text(size = rel(4.5)), 
                        legend.title = element_text(size = rel(4.7)),
                        plot.margin = unit(c(1.2, 2, 2, 1.5), 'cm'))


#### Tumor vs Adj Normal Plots ####
# Append adjusted Q values to the results file
qtmp <- qvalue(results$pvalues, fdr.level = 0.01, lambda = seq(0.0, 0.55, 0.05))

qtmp_high <- qvalue(results$pvalues, fdr.level = 0.05, lambda = seq(0.0, 0.55, 0.05))

# Get Q value cutoffs to draw lines in plots
tmp1 <- results[qtmp$significant == T, ]
tmp2 <- results[qtmp_high$significant == T, ]

# Of these cutoffs, what is the lowest p value?
qcut1 <- min(-log10(as.numeric(paste(tmp1[ ,1])))) 
qcut2 <- min(-log10(as.numeric(paste(tmp2[ ,1]))))

# Only if no sites are significant
#qcut1 = -log10(0.01)
#qcut2 = -log10(0.05)

# Select p Threshold
pThreshold <- 0.001
pComp <- as.data.frame(matrix(nrow = 3,ncol = 3))

#### Plot Results
#install.packages('ggrepel')
library(ggrepel)
BRCA <- ggplot(results, aes(as.numeric(paste(results[ ,3])), 
                                      -log10(as.numeric(paste(results[ ,1]))))) + 
     geom_point(colour = '#989898', fill = '#C7BBC9', pch = 21, size = 25) + 
     scale_color_gradient2(low = "blue", mid="grey", high = "red") +   
     labs(list(x = "Beta Coefficient", 
               y = "-log10 p Value", 
               title = "Tumor vs.  Normal")) + 
     geom_hline(yintercept = qcut1, color = "red", linetype = "dashed", size = 1.8) + 
     geom_hline(yintercept = qcut2, color = "black", linetype = "dashed", size = 1.8) +
     overall_theme +
     geom_text_repel(aes(label=rownames(results)), size=12)



png('VAE_LaWAS.png',width = 1600, height = 2400)
BRCA
dev.off()


#####################
# VAE correlation chord diagram
#####################
     # http://zuguang.de/circlize_book/book/the-chorddiagram-function.html

     #install.packages('dplyr') # Only have to run once.
     library(dplyr)
     library(reshape2)
     library(circlize)
     correlations = cor(vae_mat)
     d_cor_melt <- arrange(melt(correlations), -abs(value))
     d_cor_melt = data.frame(d_cor_melt)
     
     #write.csv(d_cor_melt, '../results/node_correlations.csv')
     
     d_cor_melt$Var1 = as.character(d_cor_melt$Var1)
     d_cor_melt$Var2 = as.character(d_cor_melt$Var2)
     
     cor_thresh = 0.8
     cor.sub = d_cor_melt[(d_cor_melt$value > cor_thresh | d_cor_melt$value < -1*cor_thresh) &
                          (d_cor_melt$value < 1 & d_cor_melt$value > -1), ]
     
     
     png('VAE_ChordPlot.png', width = 2000, height = 2000, res = 300)
     
     col_fun = colorRamp2(range(cor.sub$value), c("yellow", "blue"), 
                          transparency = 0.5)
     
     chordDiagram(cor.sub, 
                  transparency = 0.4,
                  order = as.character(seq(1:100)),
                  col = col_fun,
                  link.sort = TRUE, link.decreasing = TRUE,
                  annotationTrack = c('name', "grid"))
     
     title(paste("Correlations between VAE nodes (corr > abs(", cor_thresh, '))',  sep = ''), cex = 0.8)
     
     circos.clear()
     dev.off()





     