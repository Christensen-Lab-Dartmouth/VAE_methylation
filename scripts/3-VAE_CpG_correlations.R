######################
# Correlation between CpGs and dimensions
#
# Author: Alexander Titus
# Created: 03/21/2018
# Updated: 03/21/2018
######################

#####################
# Set up the environment
#####################
     require(data.table)
     require(limma)

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
     
     ## VAE nodes
     # 2D analyses
     vae.file = 'results/VAE_training/encoded_methyl_onehidden_warmup_batchnorm_100k_2.tsv'
     # 3D analyses
     vae.file = 'results/VAE_training/encoded_methyl_onehidden_warmup_batchnorm_100k_3.tsv'
     # 100D analyses
     vae.file = 'results/VAE_training/encoded_methyl_onehidden_warmup_batchnorm_100k_100.tsv'
     
     vae = data.frame(fread(vae.file))
     colnames(vae) = as.character(vae[1,])
     rownames(vae) = vae[,1]
     vae = vae[2:nrow(vae), 2:ncol(vae)]
     vae = vae[rownames(vae) %in% covs.updated$Basename,]
     
     vae = vae[order(rownames(vae), decreasing=T), ]
     
     ## Check sample concordance
     all(rownames(vae) == rownames(betas))
     all(covs.updated$Basename == rownames(betas))
     
     
#####################
# Correlations
#####################  

     node_corrs = function(node, betas, vae) {
          vaeNode = vae[, as.character(node)]
          
          cor.func = function(x){return(cor(x, vaeNode, method = 'spearman'))}
          
          correlations = apply(betas, 2, cor.func)
          correlations = data.frame(correlations)
          correlations = cbind('CpG' = rownames(correlations), correlations)
          
          nodeLabel = paste('VAE', node, sep = '')
          correlations = cbind(nodeLabel = rep(nodeLabel, nrow(correlations)), correlations)
          results = correlations[order(abs(correlations$correlations), decreasing=T), ]
          
          return(results)
     }

     subset_cors = function(threshold, correlations){
          cor.sub = correlations[(correlations$correlations >= threshold | 
                                       correlations$correlations < -1*threshold), ]
          return(cor.sub)
     }
     
     threshold = 0.5

#####################  
     # 2D analyses
     correlations1 = node_corrs(1, betas, vae)
     correlations1 = subset_cors(threshold, correlations1)
     
     correlations2 = node_corrs(2, betas, vae)
     correlations2 = subset_cors(threshold, correlations2)  
     
     # 3D analyses
     correlations1 = node_corrs(1, betas, vae)
     correlations1 = subset_cors(threshold, correlations1)
     
     correlations2 = node_corrs(2, betas, vae)
     correlations2 = subset_cors(threshold, correlations2)
     
     correlations3 = node_corrs(3, betas, vae)
     correlations3 = subset_cors(threshold, correlations3)
     
     # 100D analyses
     # ER-
     correlations24 = node_corrs(24, betas, vae)
     correlations24 = subset_cors(threshold, correlations24)
     
     correlations35 = node_corrs(35, betas, vae)
     correlations35 = subset_cors(threshold, correlations35)
     
     correlations43 = node_corrs(43, betas, vae)
     correlations43 = subset_cors(threshold, correlations43)
     
     # ER+
     correlations47 = node_corrs(47, betas, vae)
     correlations47 = subset_cors(threshold, correlations47)
     
     correlations91 = node_corrs(91, betas, vae)
     correlations91 = subset_cors(threshold, correlations91)
     
     correlations93 = node_corrs(93, betas, vae)
     correlations93 = subset_cors(threshold, correlations93)
     
     # ER+/-
     correlations63 = node_corrs(63, betas, vae)
     correlations63 = subset_cors(threshold, correlations63)
     
     correlations37 = node_corrs(37, betas, vae)
     correlations37 = subset_cors(threshold, correlations37)
     
     correlations22 = node_corrs(22, betas, vae)
     correlations22 = subset_cors(threshold, correlations22)
     
               
#####################
# Correlation diagram
##################### 
     plot_corr_elbow = function(correlationsN, node, threshold, plot_line = F, line = NA){
          ## Correlation elbow
          cors = correlationsN$correlations
          file.name = paste('results/correlation_elbow_node', node, '.png', sep = '')
          
          png(file.name, width = 1000, height = 1000, res = 100)
          plot(abs(cors), main=paste('Node ', node, sep = ''), 
               ylab = '|Correlation|', xlab='Index', ylim = c(threshold, 0.9))
          if(plot_line == T){
               abline(h = line)   
               abline(v = 1000, col = "black")
               abline(v = 5000, col = "green")
          }
          dev.off()
     }

##################### 
     # 2D analyses
     plot_corr_elbow(correlations1, 1, threshold, plot_line = T, line = 0.6)
     plot_corr_elbow(correlations2, 2, threshold, plot_line = T, line = 0.8)
     
     # 3D analyses
     plot_corr_elbow(correlations1, 1, threshold, plot_line = T, line = 0.65)
     plot_corr_elbow(correlations2, 2, threshold, plot_line = T, line = 0.75)
     plot_corr_elbow(correlations3, 3, threshold, plot_line = T, line = 0.7)
     
     # 100D analyses
     plot_corr_elbow(correlations24, 24, threshold, plot_line = T, line = 0.7)
     plot_corr_elbow(correlations35, 35, threshold, plot_line = T, line = 0.625)
     plot_corr_elbow(correlations43, 43, threshold, plot_line = T, line = 0.7)
     plot_corr_elbow(correlations47, 47, threshold, plot_line = T, line = 0.55)
     plot_corr_elbow(correlations91, 91, threshold, plot_line = T, line = 0.625)
     plot_corr_elbow(correlations93, 93, threshold, plot_line = T, line = 0.6)
     plot_corr_elbow(correlations63, 63, threshold, plot_line = T, line = 0.65)
     plot_corr_elbow(correlations37, 37, threshold, plot_line = T, line = 0.6)
     plot_corr_elbow(correlations22, 22, threshold, plot_line = T, line = 0.65)

     
#####################
# Genomic context
#####################
     ## https://www.bioconductor.org/help/workflows/methylationArrayAnalysis/
     #install.packages('matrixStats')
     #source("https://bioconductor.org/biocLite.R")
     #biocLite('missMethyl')
     library(missMethyl)
     library(Gviz)
     library(minfi)
     library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     
     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     head(ann450k)
     
     go_pathway_analysis = function(correlationsN, node, threshold, annotations){
          correlations = correlationsN
          
          cor.sub = correlations[(correlations$correlations >= threshold | 
                                       correlations$correlations < -1*threshold), ]
          
          cpgs = rownames(cor.sub)
          cpg_set = cor.sub
          cpg_set = rbind(cpg_set, 1, cor.sub)
          
          cpgs = cpg_set$CpG     
          
          ann450k = annotations
          
          anno.sub = ann450k[rownames(ann450k) %in% cpgs, ]
          anno.sub = data.frame(anno.sub)
          anno.sub = cbind('NodeCor' = cor.sub$correlations, anno.sub)
          
          file.name = paste('results/anno450K_node', node, '.csv', sep = '')
          write.csv(anno.sub, file.name)     
          
          
          ## GSA & GO analysis
          # load Broad human curated (C2) gene sets
          # http://bioinf.wehi.edu.au/software/MSigDB/
          load('human_c2_v5p2.rdata')
          
          # analysis
          par(mfrow=c(1,1))
          all = colnames(betas)
          cpgs = as.character(cpgs)
          gst <- gometh(sig.cpg=cpgs, all.cpg=all, plot.bias=TRUE)
          gsa <- gsameth(sig.cpg=cpgs, all.cpg=all, collection=Hs.c2)
          
          # Top 10 GO categories
          go = topGO(gst, number = 50)
          
          # Top 10 gene sets
          gsa = topGSA(gsa, number=50)
          
          gofile.name = paste('results/go_anno_node', node, '.csv', sep = '')
          write.csv(go, gofile.name) 
          
          gsafile.name = paste('results/gsa_anno_node', node, '.csv', sep = '')
          write.csv(gsa, gsafile.name) 
     }
     

#####################
     # 2D analyses
     go_pathway_analysis(correlations1, 1, threshold = 0.6, ann450k)
     go_pathway_analysis(correlations2, 2, threshold = 0.8, ann450k)
     
     # 3D analyses
     go_pathway_analysis(correlations1, 1, threshold = 0.65, ann450k)
     go_pathway_analysis(correlations2, 2, threshold = 0.75, ann450k)
     go_pathway_analysis(correlations3, 3, threshold = 0.7, ann450k)
     
     # 100D analyses
     go_pathway_analysis(correlations24, 24, threshold = 0.7, ann450k)
     go_pathway_analysis(correlations35, 35, threshold = 0.625, ann450k)
     go_pathway_analysis(correlations43, 43, threshold = 0.7, ann450k)
     go_pathway_analysis(correlations47, 47, threshold = 0.55, ann450k)
     go_pathway_analysis(correlations91, 91, threshold = 0.625, ann450k)
     go_pathway_analysis(correlations93, 93, threshold = 0.6, ann450k)
     go_pathway_analysis(correlations63, 63, threshold = 0.65, ann450k)
     go_pathway_analysis(correlations37, 37, threshold = 0.6, ann450k)
     go_pathway_analysis(correlations22, 22, threshold = 0.65, ann450k)
     
     
#####################
# Enhancer calculations
#####################

     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     breast.anno = data.frame(fread('data/Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv'))
     ann450k = merge(ann450k, breast.anno, by = 'Name')
     rownames(ann450k) = ann450k$Name
     ann450k = ann450k[rownames(ann450k) %in% colnames(betas), ]
     
     enhancer_analysis = function(correlations_list, node_list, threshold, annotations, enhancer = 'Enhancer'){
          
          enhancerResults = c('Node', 'Est', 'Conf95low', 'Conf95high', 'Pvalue')
          
          for( i in 1:length(correlations_list) ){
               correlations = correlations_list[[i]]
               
               nodeName = paste('VAE', node_list[i], sep='')

               cor.sub = correlations[(correlations$correlations >= threshold | 
                                            correlations$correlations < -1*threshold), ]
               
               
               cpgs = rownames(cor.sub)
               cpg_set = cor.sub
               cpg_set.unique = cpg_set[rownames(cpg_set) %in% unique(cpg_set$CpG), ]
               
               ann450k = annotations
               
               anno.sub = ann450k[rownames(ann450k) %in% cpg_set.unique[cpg_set.unique$nodeLabel == nodeName, ]$CpG, ]
               anno.sub = data.frame(anno.sub)
               anno.sub = anno.sub[order(rownames(anno.sub)), ]
               
               cpg_set.unique = cpg_set.unique[order(cpg_set.unique$CpG), ]
               
               temp = cpg_set.unique[cpg_set.unique$nodeLabel == nodeName, ]
               all(temp$CpG == rownames(anno.sub))
               
               anno.sub = cbind('NodeCor' = temp$correlations, anno.sub)
               
               if(enhancer == 'Enhancer'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Enhancer == 'TRUE', ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Enhancer == '', ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Enhancer == 'TRUE', ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Enhancer == '', ]); annoNoEnhancer
               } else{
                    nodeEnhancer = nrow(anno.sub[anno.sub[[enhancer]] == 1, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub[[enhancer]] == 0, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k[[enhancer]] == 1, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k[[enhancer]] == 0, ]); annoNoEnhancer
               }
               
               
               
               enhancers <- matrix(c(nodeEnhancer, annoEnhancer, 
                                     nodeNoEnhancer, annoNoEnhancer), nrow = 2,
                                   dimnames =
                                        list(c("NodeRelated", "NotNodeRelated"),
                                             c("Enhancer", "NotEnhancer")))
               
               node.results = fisher.test(enhancers); node.results
               
               
               
               row = c(nodeName, node.results$estimate, 
                       node.results$conf.int[1],
                       node.results$conf.int[2],
                       node.results$p.value)
               
               enhancerResults = rbind(enhancerResults, row)
          }
          
          enhancerResults = data.frame(enhancerResults)
          enhancerResults2 = enhancerResults
          colnames(enhancerResults2) = as.character(unlist(enhancerResults2[1, ]))
          enhancerResults2 = enhancerResults2[2:nrow(enhancerResults2), ]
          enhancerResults2$Est = as.numeric(as.character(enhancerResults2$Est))
          enhancerResults2$Conf95low = as.numeric(as.character(enhancerResults2$Conf95low))
          enhancerResults2$Conf95high = as.numeric(as.character(enhancerResults2$Conf95high))
          enhancerResults2$Pvalue = as.numeric(as.character(enhancerResults2$Pvalue))
          
          return(enhancerResults2)
     }
  
     
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
                       legend.key = element_rect(colour = NA),
                       legend.position = "bottom",
                       legend.direction = "horizontal",
                       legend.key.size= unit(0.2, "cm"),
                       legend.margin = unit(0, "cm"),
                       legend.title = element_text(face="italic"),
                       plot.margin=unit(c(10,5,5,5),"mm"),
                       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                       strip.text = element_text(face="bold")
               ))
          
     }
     
     scale_fill_Publication <- function(...){
          library(scales)
          discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
          
     }
     
     scale_colour_Publication <- function(...){
          library(scales)
          discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
          
     }     
     
#####################
# Enhancer analyses
#####################
     library(ggplot2)
     enhancer.names = colnames(ann450k)[47:63]
     enhancer.names = c('Enhancer', enhancer.names)
     enhancer_results.out = data.frame(Node = 0, Est = 0, Conf95low = 0, 
                                       Conf95high = 0, Pvalue = 0, Mark = 0)
     for(name in enhancer.names){
     # 100D analyses
     correlations_list = list(correlations24, correlations35, correlations43,
                              correlations91, correlations93,
                              correlations63, correlations37, correlations22)
     node_list = c(24, 35, 43, 91, 93, 63, 37, 22)
     
     enhancer = name
     enhancer_results = enhancer_analysis(correlations_list, node_list, threshold, ann450k, 
                                          enhancer = enhancer)
     enhancer_results$Node = factor(enhancer_results$Node, 
                                    levels = enhancer_results$Node[order(as.numeric(substr(enhancer_results$Node, 
                                                                                           4, 100)))])
     enhancer_results$Mark = rep(name, nrow(enhancer_results))
     enhancer_results.out = rbind(enhancer_results.out, enhancer_results)
     
     # 100D analyses
     enhancer_results$Node = factor(enhancer_results$Node, levels = c("VAE24", "VAE35", "VAE43",
                                                                      "VAE91", "VAE93",
                                                                      "VAE22", "VAE37", "VAE63"))

     fp <- ggplot(data=enhancer_results, aes(x=Est, 
                                            y=factor(Node), 
                                            xmin=Conf95low, 
                                            xmax=Conf95high)) +
          geom_point(color = c('black', 'black', 'black',
                               'lightgrey', 'lightgrey',
                               'orange', 'orange', 'orange')) +
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-1) +
          geom_errorbarh(height=.02) +
          ylab('VAE latent dimension') +
          geom_vline(xintercept=1, color='black', linetype='dashed') +
          scale_x_continuous(limits=c(0, max(enhancer_results$Conf95high)), name='Odds ratio w/ 95% CI') +
          ggtitle(paste0(enhancer)) + theme_Publication()
     
     
     # 100D analyses
     file.name = paste0('results/OR_enhancer_by_node_100D_', enhancer, '.png')
     png(file.name, width = 2000, height = 2500, res = 300)
     print(fp)
     dev.off()
     }
     
     enhancer_results.out = enhancer_results.out[2:nrow(enhancer_results.out), ]
     write.csv(enhancer_results.out, file = 'results/enhancer_results_by_mark_100D.csv')

     
#####################
# Visualize enhancer results
#####################
     library(tidyverse)
     library(ComplexHeatmap)
     library(ggdendro)
     enhancer_results.out$Est = as.numeric(enhancer_results.out$Est)
     heat.data = enhancer_results.out[, c('Mark', 'Node', 'Est')]
    
     # Elaboration of heatmap (white - steelblue)
     # Run clustering
     heat.data2 = heat.data %>% spread(Mark, Est) 
     rownames(heat.data2) = rowsnames = heat.data2[,1]
     heat.data2 = heat.data2[, 2:ncol(heat.data2)]
     heat.data2 = as.matrix(heat.data2)
     heat.data2 = apply(heat.data2, 2, as.numeric)
     rownames(heat.data2) = rowsnames
     
     
     file.name = paste0('results/OR_heatmap_100D.png')
     png(file.name, width = 2000, height = 2000, res = 300)
     Heatmap(t(heat.data2), col = colorRamp2(c(0, 1, 6), c("red", "white", "black")),
             heatmap_legend_param = list(color_bar = 'continuous', title = 'OR'), 
             cluster_rows = T)
     dev.off()
     

     
#####################
# Visualize genomic context
#####################
     #biocLite('cirlize')
     library(circlize)
     
     ## Visualization
     #http://zuguang.de/circlize_book/book/high-level-genomic-functions.html#genomic-heatmap

     cpg_set = rbind(correlations63, correlations37, correlations22)
     cpg_set.unique = cpg_set[rownames(cpg_set) %in% unique(cpg_set$CpG), ]
     
     anno.sub = ann450k[which(rownames(ann450k) %in% cpg_set.unique$CpG), ]
     anno.sub = data.frame(anno.sub)
     
     anno.sub = anno.sub[order(rownames(anno.sub)), ]
     cpg_set.unique = cpg_set.unique[order(rownames(cpg_set.unique)), ]
     all(rownames(anno.sub) == rownames(cpg_set.unique))
     
     bed = cbind(anno.sub$chr.x, 
                 anno.sub$pos.x, 
                 anno.sub$pos.x + 10,
                 cpg_set.unique$correlations,
                 rownames(anno.sub),
                 rownames(cpg_set.unique),
                 as.character(cpg_set.unique$nodeLabel))
     bed = data.frame(bed)
     colnames(bed) = c('chr', 'start', 'end', 'value1', 'AnnoCpG', 'CorrCpG', 'Node')
     all(bed$AnnoCpG == bed$CorrCpG)
     
     bed$start = as.numeric(as.character(bed$start))
     bed$end = as.numeric(as.character(bed$end))
     bed$value1 = as.numeric(as.character(bed$value1))
     
     bed = bed[with(bed, order(bed$chr)), ]
     
     bed_list = list('VAE63' = bed[bed$Node == 'VAE63', ],
                     'VAE37' = bed[bed$Node == 'VAE37', ],
                     'VAE22' = bed[bed$Node == 'VAE22', ])
     
     circlize_plot = function() {
          circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22))
          circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("lightsteelblue", 
                                                                        "orangered1",
                                                                        "#0000FF80"))
          
          circos.genomicDensity(bed[bed$Node == 'VAE63', ], col = c("lightsteelblue"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE37', ], col = c("orangered1"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE22', ], col = c("#0000FF80"), track.height = 0.1)
          circos.clear()
     }
     
     
     ## Corr > 0.6
     png('results/genomic_context_corr_ERboth_100D.png', width = 2000, height = 2000, res = 300)
     circlize_plot()
     dev.off()
     
     out.data = merge(anno.sub, bed, by.x = 'Name', by.y = 'AnnoCpG')
     write.csv(out.data, file = 'results/genomic_context_anno_corr_3D.csv')
   
     
     
#####################
# multinomial regression
#####################
     # node_list = c(24, 35, 43, 47, 91, 93, 63, 37, 22)
     library("nnet")
     
     set.seed(100)
     
     vae.temp = vae
     vae.temp$Basename = rownames(vae)
     temp = merge(vae.temp, covs.updated[, c("Basename", "ER")], by = 'Basename')
     temp = temp[temp$ER != '', ]
     temp$ER = ifelse(temp$ER == 'Positive', 1, 0)
     
     temp.pos = temp[temp$ER == 1, ]
     temp.neg = temp[temp$ER == 0, ]
     
     trainingRows <- sample(1:nrow(temp.pos), 0.7*nrow(temp.pos))
     training.pos <- temp.pos[trainingRows, ]
     test.pos <- temp.pos[-trainingRows, ]
     
     trainingRows <- sample(1:nrow(temp.neg), 0.7*nrow(temp.neg))
     training.neg <- temp.neg[trainingRows, ]
     test.neg <- temp.neg[-trainingRows, ]
     
     training = rbind(training.neg, training.pos)
     test = rbind(test.neg, test.pos)
     
     train <- multinom(ER ~ `24` + `35` + `43` + `47` + `91` + `93`, data = training)

     summary(train)
     z <- summary(train)$coefficients/summary(train)$standard.errors
     z
     p <- (1 - pnorm(abs(z), 0, 1))*2
     p
     exp(coef(train))
     
     pred <- predict (train, test, "probs") # predict on new data
     pred
     pred_class <- predict(train, test)
     pred_class
     
     table(pred_class, test$ER)
     mean(as.character(pred_class) == as.character(test$ER))
     
     
#####################
# 
#####################
     # http://r-statistics.co/Logistic-Regression-With-R.html
     fit = glm(ER ~ `24` + `35` + `43` + `47` + `91` + `93`, 
               data = training, family = binomial(link = "logit"))
     summary(fit)
     predicted = predict(fit, test, type="response")
     
     library(InformationValue)
     optCutOff <- optimalCutoff(test$ER, predicted)[1] 
     
     misClassError(test$ER, predicted, threshold = optCutOff)
     roc = plotROC(test$ER, predicted, returnSensitivityMat = T)
     
     write.csv(roc, 'results/ER_status_classification_ROC.csv')
     
     Concordance(test$ER, predicted)
     
     sensitivity(test$ER, predicted, threshold = optCutOff)
     specificity(test$ER, predicted, threshold = optCutOff)
     confusionMatrix(test$ER, predicted, threshold = optCutOff)
     
     