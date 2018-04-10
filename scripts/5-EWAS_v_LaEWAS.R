######################
# Comparing significant CpGs from EWAS to LaWAS results
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
     base.dir = 'C:/Users/atitus/github/VAE_analysis'
     setwd(base.dir)
     

######################
## Read in data
     
     # Methylation data
     beta.file = 'TCGA_BRCA_Betas.tsv'
     beta.dir = paste('data', beta.file, sep = '/')
     betas = data.frame(fread(beta.dir)) # on my computer takes ~8min
     rownames(betas) = betas[,1]
     betas = betas[,2:ncol(betas)]
     betas = betas[order(rownames(betas), decreasing=T), ]
     
######################
     # Tumor v. normal
     # EWAS
     data.ewas = data.frame(fread('results/EWAS/TumorvNorm_EWAS.csv'))
     
     # Latent dimensions
     # node_list = c(1, 2, 67, 86, 28, 33, 84, 95)
     data.1 = data.frame(fread('results/anno450K_node1.csv'))
     data.1$Dim = rep(1, nrow(data.1))
     
     data.2 = data.frame(fread('results/anno450K_node2.csv'))
     data.2$Dim = rep(2, nrow(data.2))
     
     data.67 = data.frame(fread('results/anno450K_node67.csv'))
     data.67$Dim = rep(67, nrow(data.67))
     
     data.86 = data.frame(fread('results/anno450K_node86.csv'))
     data.86$Dim = rep(86, nrow(data.86))
     
     data.28 = data.frame(fread('results/anno450K_node28.csv'))
     data.28$Dim = rep(28, nrow(data.28))
     
     data.33 = data.frame(fread('results/anno450K_node33.csv'))
     data.33$Dim = rep(33, nrow(data.33))
     
     data.84 = data.frame(fread('results/anno450K_node84.csv'))
     data.84$Dim = rep(84, nrow(data.84))
     
     data.95 = data.frame(fread('results/anno450K_node95.csv'))
     data.95$Dim = rep(95, nrow(data.95))
     
     data.all = rbind(data.1, data.2, data.86, data.28, data.33, data.84, data.95)
     write.csv(data.all, file = 'results/data_all.csv')

     
######################
     # Tumor v. normal
     data.1 = data.frame(fread('results/anno450K_node1.csv'))
     data.1$Dim = rep(1, nrow(data.1))
     
     data.56 = data.frame(fread('results/anno450K_node56.csv'))
     data.56$Dim = rep(56, nrow(data.56))
     
     data.65 = data.frame(fread('results/anno450K_node65.csv'))
     data.65$Dim = rep(65, nrow(data.65))
     
     data.85 = data.frame(fread('results/anno450K_node85.csv'))
     data.85$Dim = rep(85, nrow(data.85))
     
     data.16 = data.frame(fread('results/anno450K_node16.csv'))
     data.16$Dim = rep(16, nrow(data.16))
     
     data.58 = data.frame(fread('results/anno450K_node58.csv'))
     data.58$Dim = rep(58, nrow(data.58))
     
     data.70 = data.frame(fread('results/anno450K_node70.csv'))
     data.70$Dim = rep(70, nrow(data.70))
     
     data.98 = data.frame(fread('results/anno450K_node98.csv'))
     data.98$Dim = rep(98, nrow(data.98))
     
     data.ERpos = rbind(data.1, data.56, data.65, data.85)
     write.csv(data.ERpos, file = 'results/dataERpos.csv')
     
     data.ERneg = rbind(data.16, data.58, data.70, data.98)
     write.csv(data.ERneg, file = 'results/dataERneg.csv')
     
     
######################
## Subset the data
     data.ewas = data.ewas[order(data.ewas$EWASqvalues, decreasing = F), ]
     data.ewas = cbind('ID' = seq.int(nrow(data.ewas)), data.ewas)
     data.esub = data.ewas[1:1000, ]
     
     
######################
## Compare CpG sets
     
     # Tumor
     data.1sub = data.1[, 1:2]
     data.1sub = merge(data.1sub, data.ewas, by = 'V1')
     
     data.2sub = data.2[, 1:2]
     data.2sub = merge(data.2sub, data.ewas, by = 'V1')
     
     data.67sub = data.67[, 1:2]
     data.67sub = merge(data.67sub, data.ewas, by = 'V1')
     
     data.86sub = data.86[, 1:2]
     data.86sub = merge(data.86sub, data.ewas, by = 'V1')
     
     # Normal
     data.28sub = data.28[, 1:2]
     data.28sub = merge(data.28sub, data.ewas, by = 'V1')
     
     data.33sub = data.33[, 1:2]
     data.33sub = merge(data.33sub, data.ewas, by = 'V1')
     
     data.84sub = data.84[, 1:2]
     data.84sub = merge(data.84sub, data.ewas, by = 'V1')
     
     data.95sub = data.95[, 1:2]
     data.95sub = merge(data.95sub, data.ewas, by = 'V1')
     
     
#####################
# Enhancer calculations
#####################
     
     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     breast.anno = data.frame(fread('data/Illumina-Human-Methylation-450kilmn12-hg19.annotated.csv'))
     ann450k = merge(ann450k, breast.anno, by = 'Name')
     rownames(ann450k) = ann450k$Name
     ann450k = data.frame(ann450k)
     ann450k = ann450k[rownames(ann450k) %in% rownames(betas), ]
     
     enhancer_analysis = function(correlations_list, node_list, threshold, annotations, mark){
          
          enhancerResults = c('Node', 'Est', 'Conf95low', 'Conf95high', 'Pvalue')
          
          for( i in 1:length(correlations_list) ){
               correlations = correlations_list[[i]]
               
               nodeName = paste('VAE', node_list[i], sep='')
               
               cpgs = correlations$V1
               
               ann450k = annotations
               
               anno.sub = ann450k[rownames(ann450k) %in% cpgs, ]
               anno.sub = data.frame(anno.sub)
               anno.sub = anno.sub[order(rownames(anno.sub)), ]
               
               nodeEnhancer = nrow(anno.sub[anno.sub[[mark]] == 'TRUE', ]); nodeEnhancer
               nodeNoEnhancer = nrow(anno.sub[anno.sub[[mark]] == '', ]); nodeNoEnhancer
               annoEnhancer = nrow(ann450k[ann450k[[mark]] == 'TRUE', ]); annoEnhancer
               annoNoEnhancer = nrow(ann450k[ann450k[[mark]] == '', ]); annoNoEnhancer
               
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
     
     ## lists must be in the same order
     correlations_list = list(data.1, data.2, data.67, data.86,
                              data.28, data.33, data.84, data.95,
                              data.ewas)

     node_list = c(1, 2, 67, 86, 28, 33, 84, 95, 'EWAS')
     
     enhancer_results = enhancer_analysis(correlations_list, node_list, threshold, ann450k, 'Enhancer')
     enhancer_results$Node = factor(enhancer_results$Node, 
                                    levels = enhancer_results$Node[order(as.numeric(substr(enhancer_results$Node, 4, 100)))])
     View(enhancer_results)
     

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
# Enhancer OR
##################### 
     enhancer_results$Node = factor(enhancer_results$Node, levels = c("VAE1","VAE2","VAE67", "VAE86",
                                                                      "VAE28", "VAE33", "VAE84", "VAE95"))
     library(ggplot2)
     fp <- ggplot(data=enhancer_results, aes(x=Est, 
                                             y=factor(Node), 
                                             xmin=Conf95low, 
                                             xmax=Conf95high)) +
          geom_point(color = 'black') +
          geom_text(aes(label=format(round(Est, 2), nsmall = 2)), hjust=0, vjust=-1) +
          geom_errorbarh(height=.02) +
          ylab('Node') +
          geom_vline(xintercept=1, color='black', linetype='dashed') +
          scale_x_continuous(limits=c(0, max(enhancer_results$Conf95high)), name='Odds ratio w/ 95% CI') +
          ggtitle('OR for enhancer in each node set of CpGs') + theme_Publication()
     
     png('results/OR_enhancer_by_node.png', width = 2000, height = 2500, res = 300)
     fp
     dev.off()