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
     #vae.file = 'results/VAE_training/encoded_methyl_onehidden_warmup_batchnorm_100k_2.tsv'
     # 3D analyses
     #vae.file = 'results/VAE_training/encoded_methyl_onehidden_warmup_batchnorm_100k_3.tsv'
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
     #correlations1 = node_corrs(1, betas, vae)
     #correlations1 = subset_cors(threshold, correlations1)
     
     #correlations2 = node_corrs(2, betas, vae)
     #correlations2 = subset_cors(threshold, correlations2)  
     
     # 3D analyses
     #correlations1 = node_corrs(1, betas, vae)
     #correlations1 = subset_cors(threshold, correlations1)
     
     #correlations2 = node_corrs(2, betas, vae)
     #correlations2 = subset_cors(threshold, correlations2)
     
     #correlations3 = node_corrs(3, betas, vae)
     #correlations3 = subset_cors(threshold, correlations3)
     
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
          plot(abs(cors), main=paste('Dimension ', node, sep = ''), 
               ylab = '|Correlation|', xlab='Index', ylim = c(threshold, 0.9))
          if(plot_line == T){
               abline(h = line)   
               abline(v = 1000, col = "green")
          }
          dev.off()
     }

##################### 
     dimensions = c(24, 35, 43, 91, 93, 22, 63)
     thresholds = c(0.7, 0.6, 0.7, 0.625, 0.6, 0.625, 0.675)
     #thresholds = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6)
     
     # 2D analyses
     #plot_corr_elbow(correlations1, 1, threshold, plot_line = T, line = 0.6)
     #plot_corr_elbow(correlations2, 2, threshold, plot_line = T, line = 0.8)
     
     # 3D analyses
     #plot_corr_elbow(correlations1, 1, threshold, plot_line = T, line = 0.65)
     #plot_corr_elbow(correlations2, 2, threshold, plot_line = T, line = 0.75)
     #plot_corr_elbow(correlations3, 3, threshold, plot_line = T, line = 0.7)
     
     # 100D analyses
     plot_corr_elbow(correlations24, 24, threshold, plot_line = T, line = thresholds[1])
     plot_corr_elbow(correlations35, 35, threshold, plot_line = T, line = thresholds[2])
     plot_corr_elbow(correlations43, 43, threshold, plot_line = T, line = thresholds[3])
     plot_corr_elbow(correlations91, 91, threshold, plot_line = T, line = thresholds[4])
     plot_corr_elbow(correlations93, 93, threshold, plot_line = T, line = thresholds[5])
     plot_corr_elbow(correlations22, 22, threshold, plot_line = T, line = thresholds[6])
     plot_corr_elbow(correlations63, 63, threshold, plot_line = T, line = thresholds[7])
     
     
     
     
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
     #go_pathway_analysis(correlations1, 1, threshold = 0.6, ann450k)
     #go_pathway_analysis(correlations2, 2, threshold = 0.8, ann450k)
     
     # 3D analyses
     #go_pathway_analysis(correlations1, 1, threshold = 0.65, ann450k)
     #go_pathway_analysis(correlations2, 2, threshold = 0.75, ann450k)
     #go_pathway_analysis(correlations3, 3, threshold = 0.7, ann450k)
     
     # 100D analyses
     go_pathway_analysis(correlations24, dimensions[1], threshold = thresholds[1], ann450k)
     go_pathway_analysis(correlations35, dimensions[2], threshold = thresholds[2], ann450k)
     go_pathway_analysis(correlations43, dimensions[3], threshold = thresholds[3], ann450k)
     go_pathway_analysis(correlations91, dimensions[4], threshold = thresholds[4], ann450k)
     go_pathway_analysis(correlations93, dimensions[5], threshold = thresholds[5], ann450k)
     go_pathway_analysis(correlations22, dimensions[6], threshold = thresholds[6], ann450k)
     go_pathway_analysis(correlations63, dimensions[7], threshold = thresholds[7], ann450k)
     
     
#####################
# Enhancer calculations
#####################

     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     breast.anno = data.frame(fread('data/Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv'))
     ann450k = merge(ann450k, breast.anno, by = 'Name')
     rownames(ann450k) = ann450k$Name
     ann450k = ann450k[rownames(ann450k) %in% colnames(betas), ]
     
     enhancer_analysis = function(correlations_list, node_list, thresholds, annotations, enhancer = 'Enhancer'){
          
          enhancerResults = c('Node', 'Est', 'Conf95low', 'Conf95high', 'Pvalue')
          
          for( i in 1:length(correlations_list) ){
               correlations = correlations_list[[i]]
               
               nodeName = paste('VAE', node_list[i], sep='')
               threshold = thresholds[i]
               
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
               } else if (enhancer == 'OpenSea'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island == enhancer, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island != enhancer, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Relation_to_Island == enhancer, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Relation_to_Island != enhancer, ]); annoNoEnhancer
               } else if (enhancer == 'Island'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island == enhancer, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island != enhancer, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Relation_to_Island == enhancer, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Relation_to_Island != enhancer, ]); annoNoEnhancer
               } else if (enhancer == 'S_Shore'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island == enhancer, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island != enhancer, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Relation_to_Island == enhancer, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Relation_to_Island != enhancer, ]); annoNoEnhancer
               } else if (enhancer == 'N_Shore'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island == enhancer, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island != enhancer, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Relation_to_Island == enhancer, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Relation_to_Island != enhancer, ]); annoNoEnhancer
               } else if (enhancer == 'N_Shelf'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island == enhancer, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island != enhancer, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Relation_to_Island == enhancer, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Relation_to_Island != enhancer, ]); annoNoEnhancer
               } else if (enhancer == 'S_Shelf'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island == enhancer, ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$Relation_to_Island != enhancer, ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$Relation_to_Island == enhancer, ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$Relation_to_Island != enhancer, ]); annoNoEnhancer
               } else if (enhancer == 'DHS'){
                    nodeEnhancer = nrow(anno.sub[anno.sub$DHS == 'TRUE', ]); nodeEnhancer
                    nodeNoEnhancer = nrow(anno.sub[anno.sub$DHS == '', ]); nodeNoEnhancer
                    annoEnhancer = nrow(ann450k[ann450k$DHS == 'TRUE', ]); annoEnhancer
                    annoNoEnhancer = nrow(ann450k[ann450k$DHS == '', ]); annoNoEnhancer
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
                       #legend.margin = unit(0, "cm"),
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
     enhancer.names = colnames(ann450k)[46:64]
     enhancer.names = c('Enhancer', 'OpenSea', 'Island', 'N_Shore', 'S_Shore', 'N_Shelf', 'S_Shelf', 'DHS',
                        enhancer.names)
     enhancer_results.out = data.frame(Node = 0, Est = 0, Conf95low = 0, 
                                       Conf95high = 0, Pvalue = 0, Mark = 0)
     
     
#####################
# Plots by mark
#####################
     for(name in enhancer.names){
          # 100D analyses
          correlations_list = list(correlations24, correlations35, correlations43,
                                   correlations91, correlations93, 
                                   correlations22, correlations63)
          dimensions = dimensions
          thresholds = thresholds
          
          enhancer = name
          enhancer_results = enhancer_analysis(correlations_list, dimensions, thresholds, ann450k, 
                                               enhancer = enhancer)
          
          enhancer_results$Node = factor(enhancer_results$Node, 
                                         levels = enhancer_results$Node[order(as.numeric(substr(enhancer_results$Node, 
                                                                                                4, 100)))])
          enhancer_results$Mark = rep(name, nrow(enhancer_results))
          enhancer_results.out = rbind(enhancer_results.out, enhancer_results)
          
          # 100D analyses
          enhancer_results$Node = factor(enhancer_results$Node, levels = c("VAE24", "VAE35", "VAE43",
                                                                           "VAE91", "VAE93",
                                                                           "VAE22", "VAE63"))
          # # Log-odds
          # fp <- ggplot(data=enhancer_results, aes(x=log(Est), 
          #                                        y=factor(Node), 
          #                                        xmin=log(Conf95low), 
          #                                        xmax=log(Conf95high))) +
          #      geom_point() +
          #      geom_text(aes(label=format(round(log(Est), 2), nsmall = 2)), hjust=0, vjust=-0.5) +
          #      geom_errorbarh(height=.02) +
          #      ylab('VAE latent dimension') +
          #      geom_vline(xintercept=log(1), color='black', linetype='dashed') +
          #      scale_x_continuous(limits=c(min(-log(1.5), max(-log(70), min(log(enhancer_results$Conf95high)))), 
          #                                  max(log(1.5), min(log(70), max(log(enhancer_results$Conf95high))))),
          #                         name='Log-odds ratio w/ 95% CI') +
          #      ggtitle(paste0(enhancer)) + theme_Publication()
          fp <- ggplot(data=enhancer_results, aes(x = factor(Node), 
                                                  y = Est, 
                                                  ymin = Conf95low, 
                                                  ymax = Conf95high)) +
               geom_pointrange() +
               geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-0.5) +
               geom_hline(aes(fill=factor(Node)), yintercept =1, linetype=2)+
               xlab('VAE latent dimension')+ ylab("Log-Odds Ratio (95% Confidence Interval)")+
               geom_errorbar(aes(ymin=Conf95low, ymax=Conf95high), width=0.1, cex=1) + 
               coord_flip() + theme_Publication() + ggtitle(paste0(enhancer)) +
               scale_y_log10(breaks=c(0.5,1,2),position="top") 
          
          # 100D analyses
          file.name = paste0('results/OR_enhancer_by_mark_100D_', enhancer, '.png')
          png(file.name, width = 2500, height = 3000, res = 300)
          print(fp)
          dev.off()
     }
     
     enhancer_results.out = enhancer_results.out[2:nrow(enhancer_results.out), ]
     write.csv(enhancer_results.out, file = 'results/enhancer_results_by_mark_100D.csv')

     
#####################
# Plots by latent dimension
#####################
     
     for(latdim in unique(enhancer_results.out$Node)){
          
          temp = enhancer_results.out[enhancer_results.out$Node == latdim, ]
          temp$Mark = factor(temp$Mark, levels = rev(c('N_Shelf', 'N_Shore', 'Island', 'S_Shore', 'S_Shelf', 
                                                   'OpenSea', 'Enhancer', 'DHS', 'Br_myo_H3K9me3',
                                                   'Br_myo_H3K4me3', 'Br_myo_H3K36me3', 'Br_myo_H3K9ac',
                                                   
                                                   'Br_myo_H3K27me3', 'Br_myo_H3K4me1',
                                                   'HMEC_H3K4me1', 'HMEC_H3K4me2', 'HMEC_H3K4me3',
                                                   'HMEC_H3K9ac', 'HMEC_H3K27ac', 'HMEC_H3K36me3', 
                                                   'HMEC_H3K79me2', 
                                                   
                                                   'HMEC_H3K20me1', 'HMEC_H3K9me3', 'HMEC_H3K27me3', 
                                                     
                                                   'HMEC_H2A.Z',
                                                   'HMEC_SuperEnhancer', 'MCF7_SuperEnhancer')))
          
          temp2 = temp[temp$Mark == 'HMEC_SuperEnhancer' |
                       temp$Mark == 'MCF7_SuperEnhancer', ]
          
          temp = temp[temp$Mark != 'HMEC_SuperEnhancer' &
                      temp$Mark != 'MCF7_SuperEnhancer', ]
          
          
          fp <- ggplot(data=temp, aes(x = factor(Mark), 
                                                  y = Est, 
                                                  ymin = Conf95low, 
                                                  ymax = Conf95high)) +
               geom_pointrange() +
               geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-0.5) +
               geom_hline(aes(fill=factor(Mark)), yintercept =1, linetype=2)+
               xlab('Genomic Context')+ ylab("Log-Odds Ratio (95% Confidence Interval)")+
               geom_errorbar(aes(ymin=Conf95low, ymax=Conf95high), width=0.1, cex=1) + 
               coord_flip() + theme_Publication() + ggtitle(latdim) +
               scale_y_log10(breaks=c(0, 0.25, 0.5, 1, 2, 4),position="top") 
          
          # By dimension
          file.name = paste0('results/OR_enhancer_by_dimension_100D_', latdim, '.png')
          png(file.name, width = 2500, height = 3000, res = 300)
          print(fp)
          dev.off()
          
          
          ############################
          # by dimension more granular
          ############################
          # Island Context
          gran = 'IslandContext'
          islandContexts = c('N_Shelf', 'N_Shore', 'Island', 'S_Shore', 'S_Shelf', 'OpenSea')
          temp = temp[temp$Mark %in% islandContexts, ]
          
          fp <- ggplot(data=temp, aes(x = factor(Mark), 
                                      y = Est, 
                                      ymin = Conf95low, 
                                      ymax = Conf95high)) +
               geom_pointrange() +
               geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-0.5) +
               geom_hline(aes(fill=factor(Mark)), yintercept =1, linetype=2)+
               xlab('CpG Island Context')+ ylab("Log-Odds Ratio (95% Confidence Interval)")+
               geom_errorbar(aes(ymin=Conf95low, ymax=Conf95high), width=0.1, cex=1) + 
               coord_flip() + theme_Publication() + ggtitle(latdim) +
               scale_y_log10(breaks=c(0, 0.25, 0.5, 1, 2, 4),position="top") 
          
          file.name = paste0('results/OR_enhancer_by_dimension_100D_', latdim, '_', gran, '.png')
          png(file.name, width = 2500, height = 3000, res = 300)
          print(fp)
          dev.off()
     }
     
     
     # Together
     gran = 'IslandContext'
     enhancer_results.out$Node = factor(enhancer_results.out$Node, levels = c("VAE24", "VAE35", "VAE43",
                                                                              "VAE91", "VAE93",
                                                                              "VAE22", "VAE63"))
     enhancer_results.out2 = enhancer_results.out
     enhancer_results.out = enhancer_results.out[enhancer_results.out$Est > 0, ]
     
     islandContexts = c('N_Shelf', 'N_Shore', 'Island', 'S_Shore', 'S_Shelf', 'OpenSea')
     temp = enhancer_results.out[enhancer_results.out$Mark %in% islandContexts, ]
     
     fp1 <- ggplot(data=temp, aes(x = Mark, 
                                                  y = Est, 
                                                  ymin = Conf95low, 
                                                  ymax = Conf95high)) +
          geom_pointrange(aes(col=Node)) +
          geom_hline(aes(fill=Mark), yintercept =1, linetype=2)+
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-0.5) +
          xlab('CpG Island Context')+ ylab("Log-Odds Ratio (95% Confidence Interval)")+
          geom_errorbar(aes(ymin=Conf95low, ymax=Conf95high, col=Node), width=0.1, cex=1) + 
          facet_wrap(~Node, strip.position="left", nrow=7, scales = "free_y") +
          coord_flip() + theme_Publication() + theme(legend.position="none", text = element_text(size = 20)) +
          scale_y_log10(breaks=c(0, 0.25, 0.5, 1, 2, 4),position="top")
     
     file.name = paste0('results/OR_enhancer_by_dimension_100D_', gran, '.png')
     png(file.name, width = 3000, height = 3000, res = 300)
     print(fp1)
     dev.off() 
     
     # Transcriptional Repression
     enhancer_results.out$Mark = factor(enhancer_results.out$Mark, levels = rev(c('N_Shelf', 'N_Shore', 'Island', 'S_Shore', 'S_Shelf', 
                                                                                  'OpenSea', 'Enhancer', 'DHS', 'Br_myo_H3K9me3',
                                                                                  'Br_myo_H3K4me3', 'Br_myo_H3K36me3', 'Br_myo_H3K9ac',
                                                                                  
                                                                                  'Br_myo_H3K27me3', 'Br_myo_H3K4me1',
                                                                                  'HMEC_H3K4me1', 'HMEC_H3K4me2', 'HMEC_H3K4me3',
                                                                                  'HMEC_H3K9ac', 'HMEC_H3K27ac', 'HMEC_H3K36me3', 
                                                                                  'HMEC_H3K79me2', 
                                                                                  
                                                                                  'HMEC_H3K20me1', 'HMEC_H3K9me3', 'HMEC_H3K27me3', 
                                                                                  
                                                                                  'HMEC_H2A.Z',
                                                                                  'HMEC_SuperEnhancer', 'MCF7_SuperEnhancer')))
     temp.colInActivity = c('HMEC_H3K20me1', 'HMEC_H3K9me3', 'HMEC_H3K27me3')
     temp.colInActivity2 = c('Br_myo_H3K4me3', 'Br_myo_H3K36me3', 'Br_myo_H3K4me1')
     context = c(temp.colInActivity, temp.colInActivity2)

     gran = 'Repressive'
     enhancer_results.out$Node = factor(enhancer_results.out$Node, levels = c("VAE24", "VAE35", "VAE43",
                                                                              "VAE91", "VAE93",
                                                                              "VAE22", "VAE63"))
     temp = enhancer_results.out[enhancer_results.out$Mark %in% context, ]
     
     fp2 <- ggplot(data=temp, aes(x = Mark, 
                                 y = Est, 
                                 ymin = Conf95low, 
                                 ymax = Conf95high)) +
          geom_pointrange(aes(col=Node)) +
          geom_hline(aes(fill=Mark), yintercept =1, linetype=2)+
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-0.5) +
          xlab('Transcriptional Repressive Context')+ ylab("Log-Odds Ratio (95% Confidence Interval)")+
          geom_errorbar(aes(ymin=Conf95low, ymax=Conf95high, col=Node), width=0.1, cex=1) + 
          facet_wrap(~Node, strip.position="left", nrow=7, scales = "free_y") +
          coord_flip() + theme_Publication() + theme(legend.position="none", text = element_text(size = 20)) +
          scale_y_log10(breaks=c(0.25, 0.5, 1, 2, 4),position="top")
     
     file.name = paste0('results/OR_enhancer_by_dimension_100D_', gran, '.png')
     png(file.name, width = 3000, height = 3000, res = 300)
     print(fp2)
     dev.off()

     
     # Transcriptional Activation
     temp.colActivity = c('HMEC_H3K4me1', 'HMEC_H3K4me2', 'HMEC_H3K4me3',
                          'HMEC_H3K9ac', 'HMEC_H3K27ac', 'HMEC_H3K36me3', 'HMEC_H3K79me2')
     temp.colActivity2 = c('Br_myo_H3K9me3','Br_myo_H3K27me3', 'Br_myo_H3K9ac')
     context = c(temp.colActivity, temp.colActivity2)
     
     gran = 'Activation'
     enhancer_results.out$Node = factor(enhancer_results.out$Node, levels = c("VAE24", "VAE35", "VAE43",
                                                                              "VAE91", "VAE93",
                                                                              "VAE22", "VAE63"))
     temp = enhancer_results.out[enhancer_results.out$Mark %in% context, ]
     
     fp3 <- ggplot(data=temp, aes(x = Mark, 
                                 y = Est, 
                                 ymin = Conf95low, 
                                 ymax = Conf95high)) +
          geom_pointrange(aes(col=Node)) +
          geom_hline(aes(fill=Mark), yintercept =1, linetype=2)+
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-0.5) +
          xlab('Transcriptional Activation Context')+ ylab("Log-Odds Ratio (95% Confidence Interval)")+
          geom_errorbar(aes(ymin=Conf95low, ymax=Conf95high, col=Node), width=0.1, cex=1) + 
          facet_wrap(~Node, strip.position="left", nrow=7, scales = "free_y") +
          coord_flip() + theme_Publication() + theme(legend.position="none", text = element_text(size = 20)) +
          scale_y_log10(breaks=c(0, 0.25, 0.5, 1, 2, 4),  position="top")
     
     file.name = paste0('results/OR_enhancer_by_dimension_100D_', gran, '.png')
     png(file.name, width = 3000, height = 2500, res = 300)
     print(fp3)
     dev.off()
     
     
     ######################
     # Save combined file
     library(gridExtra)
     
     lay <- rbind(c(1,3),
                  c(2,3))
     
     g = grid.arrange(fp1, fp2, fp3, ncol = 2, layout_matrix = lay)
     print(g)
     
     ggsave(g, file="results/OR_enhancer_by_dimension_100D_combined.png", width=18, height=24)
     
     
#####################
# Visualize enhancer results
#####################
     library(tidyverse)
     library(ComplexHeatmap)
     library(ggdendro)
     library(circlize)
     enhancer_results.out$Mark = factor(enhancer_results.out$Mark, levels = rev(c('N_Shelf', 'N_Shore', 'Island', 'S_Shore', 'S_Shelf', 
                                                  'OpenSea', 'Enhancer', 'DHS', 'Br_myo_H3K9me3',
                                                  'Br_myo_H3K4me3', 'Br_myo_H3K36me3', 'Br_myo_H3K9ac',
                                                  
                                                  'Br_myo_H3K27me3', 'Br_myo_H3K4me1',
                                                  'HMEC_H3K4me1', 'HMEC_H3K4me2', 'HMEC_H3K4me3',
                                                  'HMEC_H3K9ac', 'HMEC_H3K27ac', 'HMEC_H3K36me3', 
                                                  'HMEC_H3K79me2', 
                                                  
                                                  'HMEC_H3K20me1', 'HMEC_H3K9me3', 'HMEC_H3K27me3', 
                                                  
                                                  'HMEC_H2A.Z',
                                                  'HMEC_SuperEnhancer', 'MCF7_SuperEnhancer')))
     
     
     enhancer_results.out$Est = as.numeric(enhancer_results.out$Est)
     enhancer_results.out$Significant = factor(ifelse(enhancer_results.out$Pvalue < 0.05, 'Yes', 'No'))
     enhancer_results.out$Direction = ifelse(enhancer_results.out$Significant == 'Yes' & 
                                             enhancer_results.out$Est > 1, 
                                             'Significantly Enriched', 
                                        ifelse(enhancer_results.out$Significant == 'Yes' & 
                                               enhancer_results.out$Est < 1, 
                                               'Significantly Depleted', 'Not Significant'))
     
     heat.data = enhancer_results.out[, c('Mark', 'Node', 'Est')]
    
     temp.ERneg = c('VAE24', 'VAE35', 'VAE43')
     temp.ERpos = c('VAE91', 'VAE93')
     temp.ERboth = c('VAE22', 'VAE63')
     vae.annos = c(temp.ERboth, temp.ERneg, temp.ERpos)
     vae.annos = ifelse(vae.annos %in% temp.ERneg, 'ER-negative', 
                        ifelse(vae.annos %in% temp.ERpos, 'ER-positive', 
                               ifelse(vae.annos %in% temp.ERboth, 'ER-both', 'N/A')))
     
     ta_anno =  HeatmapAnnotation(df = data.frame(VAE = vae.annos),
                                  col = list(VAE = c('ER-negative' = '#fd7b5f', 
                                                     'ER-positive' = '#fddb77', 
                                                     'ER-both' = '#c9c4fa')))
     
     # Elaboration of heatmap (white - steelblue)
     # Run clustering
     heat.data2 = heat.data %>% spread(Mark, Est) 
     rownames(heat.data2) = rowsnames = heat.data2[,1]
     heat.data2 = heat.data2[, 2:ncol(heat.data2)]
     heat.data2 = as.matrix(heat.data2)
     heat.data2 = apply(heat.data2, 2, as.numeric)
     rownames(heat.data2) = rowsnames
     
     # Heatmap with all marks
     file.name = paste0('results/OR_heatmap_100D.png')
     png(file.name, width = 3000, height = 3000, res = 300)
     h1 = Heatmap(t(heat.data2), col = colorRamp2(c(0, 1, 6), c("#E69F00", "white", "black")),
             heatmap_legend_param = list(color_bar = 'continuous', title = 'OR'), 
             cluster_rows = T, cluster_columns = T,
             top_annotation = ta_anno); h1
     dev.off()
     
     
     # Heatmap with CpG island context
     temp.cols = c('N_Shelf', 'N_Shore', 'Island', 'S_Shore', 'S_Shelf', 'OpenSea', 'Enhancer', 'DHS')
     temp = heat.data2[, colnames(heat.data2) %in% temp.cols]
     
     file.name = paste0('results/OR_heatmap_CpGislandContext_100D.png')
     png(file.name, width = 3000, height = 3000, res = 300)
     h2 = Heatmap(t(temp), col = colorRamp2(c(0, 1, 6), c("#E69F00", "white", "black")),
             heatmap_legend_param = list(color_bar = 'continuous', title = 'OR'), 
             cluster_rows = F, cluster_columns = T,
             top_annotation = ta_anno); h2
     dev.off()
     
     
     # Heatmap with enhancer marks
     temp.cols = c('Enhancer', 'DHS', 'HMEC_SuperEnhancer', 'MCF7_SuperEnhancer')
     temp = heat.data2[, colnames(heat.data2) %in% temp.cols]
     
     file.name = paste0('results/OR_heatmap_enhancers_100D.png')
     png(file.name, width = 3000, height = 3000, res = 300)
     h3 = Heatmap(t(temp), col = colorRamp2(c(0, 1, 6), c("#E69F00", "white", "black")),
             heatmap_legend_param = list(color_bar = 'continuous', title = 'OR'), 
             cluster_rows = T, cluster_columns = T,
             top_annotation = ta_anno); h3
     dev.off()
     
     
     # Heatmap with HMEC marks
     
     temp.cols = c('HMEC_H3K4me1', 'HMEC_H3K4me2', 'HMEC_H3K4me3',
                   'HMEC_H3K9ac', 'HMEC_H3K27ac', 'HMEC_H3K36me3', 'HMEC_H3K79me2', 
                   'HMEC_H3K20me1', 'HMEC_H3K9me3', 'HMEC_H3K27me3')
     temp.cols2 = c('Br_myo_H3K9me3', 'Br_myo_H3K4me3', 
                   'Br_myo_H3K36me3', 'Br_myo_H3K9ac',
                   'Br_myo_H3K27me3', 'Br_myo_H3K4me1')
     
     temp.cols = c(temp.cols, temp.cols2)
     
     temp.colActivity = c('HMEC_H3K4me1', 'HMEC_H3K4me2', 'HMEC_H3K4me3',
                          'HMEC_H3K9ac', 'HMEC_H3K27ac', 'HMEC_H3K36me3', 'HMEC_H3K79me2')
     temp.colActivity2 = c('Br_myo_H3K9me3','Br_myo_H3K27me3', 'Br_myo_H3K9ac')
     temp.colActivity = c(temp.colActivity, temp.colActivity2)
     
     temp.colInActivity = c('HMEC_H3K20me1', 'HMEC_H3K9me3', 'HMEC_H3K27me3')
     temp.colInActivity2 = c('Br_myo_H3K4me3', 'Br_myo_H3K36me3', 'Br_myo_H3K4me1')
     temp.colInActivity = c(temp.colInActivity, temp.colInActivity2)
     
     
     temp.out = heat.data[heat.data$Mark %in% temp.cols, ]
     temp.out = data.frame('Mark' = temp.out$Mark)
     temp = heat.data2[, colnames(heat.data2) %in% temp.cols]
     
     mark.annos = colnames(temp)
     mark.annos = ifelse(mark.annos %in% temp.colActivity, 'Activation', 
                         ifelse(mark.annos %in% temp.colInActivity, 'Repression', 'N/A'))
     
 
     file.name = paste0('results/OR_heatmap_HMEC_100D.png')
     png(file.name, width = 3000, height = 3000, res = 300)
          h4 = Heatmap(t(temp), col = colorRamp2(c(0, 1, 6), c("#E69F00", "white", "black")),
                  heatmap_legend_param = list(color_bar = 'continuous', title = 'OR'), 
                  cluster_rows = T, 
                  cluster_columns = T,
                  row_names_side = "left", 
                  show_row_dend = T,
                  top_annotation = ta_anno) + 
                    Heatmap(mark.annos, col = list('Activation' = 'darkgrey', 'Repression' = 'grey'),
                    name = 'Transcription'); h4
     dev.off()
     
     
     # Bubble chart
     enhancer_results.out$VAE = factor(enhancer_results.out$Node, levels = c('VAE24', 'VAE35', 'VAE43',
                                                                              'VAE91', 'VAE93',
                                                                              'VAE22', 'VAE63'))
          
     file.name = paste0('results/OR_bubble_chart_100D.png')
     png(file.name, width = 3000, height = 3000, res = 300)
     
     group.colors = c('Significantly Enriched' = 'yelow', 
                      'Significantly Depleted' = 'blue', 
                      'Not Significant' = 'green')
     
     ggplot(enhancer_results.out, aes(x = VAE, 
                                      y = Mark,
                                      col = Direction)) +
          geom_point(aes(size=Significant)) + 
          theme_Publication() + 
          ylab('Genomic Context') +
          xlab('VAE Latent Dimension') +
          guides(size=F) +
          scale_color_manual(values = c('grey', 
                                        '#E69F00', 
                                        '#56B4E9'))
     
     dev.off()
     
#####################
# Visualize genomic context
#####################
     #biocLite('cirlize')
     library(circlize)
     
     ## Visualization
     #http://zuguang.de/circlize_book/book/high-level-genomic-functions.html#genomic-heatmap

     cpg_set = rbind(correlations22, correlations63)
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
     
     bed_list = list('VAE22' = bed[bed$Node == 'VAE22', ],
                     'VAE63' = bed[bed$Node == 'VAE63', ])#,
                     #'VAE43' = bed[bed$Node == 'VAE43', ])
     
     circlize_plot = function() {
          circos.initializeWithIdeogram(species = "hg19", chromosome.index = paste0("chr", 1:22))
          circos.genomicRainfall(bed_list, pch = 16, cex = 0.4, col = c("lightsteelblue", 
                                                                        "orangered1"))#,
                                                                        #"#0000FF80"))
          
          circos.genomicDensity(bed[bed$Node == 'VAE22', ], col = c("lightsteelblue"), track.height = 0.1)
          circos.genomicDensity(bed[bed$Node == 'VAE63', ], col = c("orangered1"), track.height = 0.1)
          #circos.genomicDensity(bed[bed$Node == 'VAE43', ], col = c("#0000FF80"), track.height = 0.1)
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
#  Logistic regression
#####################
     library(InformationValue)
     library(ROCR)
     
     # load the library
     library(caret)
     
     # define training control
     train_control <- trainControl(method="LOOCV")
     
     # train the model
     model = train(factor(ER) ~ `24` + `35` + `43` + `47` + `91` + `93`, 
                    data=training, method="glm", family="binomial",
                    trControl=train_control, tuneLength = 5)
     
     pred = predict(model, newdata=test)
     confusionMatrix(data=pred, factor(test$ER))
     
     # summarize results
     print(model)
     
     
     
     # http://r-statistics.co/Logistic-Regression-With-R.html
     # Training
     fit = glm(ER ~ `24` + `35` + `43` + `47` + `91` + `93`, 
               data = training, family = binomial(link = "logit"))
     summary(fit)
     
     
     
     
     # Testing
     predicted = predict(fit, test, type="response")
     
     optCutOff <- optimalCutoff(test$ER, predicted)[1] 
     
     1 - misClassError(test$ER, predicted, threshold = optCutOff)
     roc = plotROC(test$ER, predicted, returnSensitivityMat = T)
     
     write.csv(roc, 'results/ER_status_classification_ROC.csv')
     
     colnames(roc) = c('FP', 'TP', 'Threshold')
     
     predicts = prediction(predicted, test$ER)
     perfs = performance(predicts, measure = 'auc')
     Concordance(test$ER, predicted)
     sense = sensitivity(test$ER, predicted, threshold = optCutOff)
     spec = specificity(test$ER, predicted, threshold = optCutOff)
     confusionMatrix(test$ER, predicted, threshold = optCutOff)
     
     ## Corr > 0.6
     png('results/ER_classifier_performance_100D.png', width = 2000, height = 2000, res = 300)
     
     ggplot(roc, aes(FP, TP)) +
          geom_line(size = 2, alpha = 0.7)+
          labs(x = "False Positive Rate (1-Specificity)", 
               y = "True Positive Rate (Sensitivity)") +
          annotate("text", x = 0.5, y = 0.5, label = paste0("AUC = ", round(as.numeric(perfs@y.values), digits = 2), '\n')) +
          theme_Publication()
     
     dev.off()
     


     