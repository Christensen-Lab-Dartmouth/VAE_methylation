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
library(ggpubr)
library(scatterplot3d)


#####################################
# Load and manipulate the data
#####################################   
# Covariates
covs.file = 'data/Full_data_covs.csv'
covs = data.frame(fread(covs.file))
covs.tumor = covs[covs$SampleType == 'Primary Tumor', ]
covs.normal = covs[covs$SampleType == 'Solid Tissue Normal', ]

data = data.frame(fread('../VAE_methylation/data/encoded_methyl_onehidden_warmup_batchnorm_100k_3.tsv'))
colnames(data) = data[1,]
data = data[2:nrow(data), ]
colnames(data) = c('Basename', 'V1', 'V2', 'V3')
data = merge(data, covs.tumor, by.x = 'Basename', by.y = 'Basename')
data$TripNeg = ifelse(data$ER == 'Negative' & 
                           data$PR == 'Negative' &
                           data$HER2 == 'Negative', 1, 0)

data$HER2int = ifelse(data$HER2 == 'Positive', 1, 
                      ifelse(data$HER2 == 'Negative', 0, NA))
data$ERint = ifelse(data$ER == 'Positive', 1, 
                      ifelse(data$ER == 'Negative', 0, NA))
data$PRint = ifelse(data$PR == 'Positive', 1, 
                      ifelse(data$PR == 'Negative', 0, NA))
data$StudyInt = ifelse(data$Study == 'TCGA_BRCA', 0, 
                       ifelse(data$Study == 'GSE84207_Fleischer', 1, 
                              ifelse(data$Study == 'GSE75067_Ringner', 2, NA)))

#####################################
# Visualizing
##################################### 
temp = data[data$Study != '', ]
temp = temp[!is.na(temp$Study), ]
temp = temp[temp$Study != 'Normal', ]
ggboxplot(temp,'Study',  'V1', color = 'Study')
ggboxplot(temp,'Study',  'V2', color = 'Study')
ggboxplot(temp,'Study',  'V3', color = 'Study')
model = lm(V1 ~ Study, data=temp); summary(model)

temp = temp[temp$PAM50 != '', ]
temp = temp[!is.na(temp$PAM50), ]
temp = temp[temp$PAM50 != 'Normal', ]
ggboxplot(temp,'PAM50',  'V1', color = 'PAM50')
ggboxplot(temp,'PAM50',  'V2', color = 'PAM50')
ggboxplot(temp,'PAM50',  'V3', color = 'PAM50')

temp = temp[!is.na(temp$ERint), ]
ggboxplot(temp, 'ER',  'V1', color = 'ER')
ggboxplot(temp, 'ER',  'V2', color = 'ER')
ggboxplot(temp, 'ER',  'V3', color = 'ER')
model = glm(ERint ~ V1 + V2 + V3 + Study, family=binomial(link='logit'), data=data); summary(model)
exp(cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95)))

temp = temp[!is.na(temp$PRint), ]
ggboxplot(temp, 'PR',  'V1', color = 'PR')
ggboxplot(temp, 'PR',  'V2', color = 'PR')
ggboxplot(temp, 'PR',  'V3', color = 'PR')
model = glm(PRint ~ V1 + V2 + V3 + Study, family=binomial(link='logit'), data=data); summary(model)
exp(cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95)))

temp = temp[!is.na(temp$HER2int), ]
ggboxplot(temp, 'HER2',  'V1', color = 'HER2')
ggboxplot(temp, 'HER2',  'V2', color = 'HER2')
ggboxplot(temp, 'HER2',  'V3', color = 'HER2')
model = glm(PRint ~ V1 + V2 + V3 + Study, family=binomial(link='logit'), data=data); summary(model)
exp(cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95)))

temp = temp[temp$TripNeg != '', ]
ggboxplot(temp, 'TripNeg',  'V1', color = 'TripNeg')
ggboxplot(temp, 'TripNeg',  'V2', color = 'TripNeg')
ggboxplot(temp, 'TripNeg',  'V3', color = 'TripNeg')
model = glm(TripNeg ~ V1 + V2 + V3 + Study, family=binomial(link='logit'), data=data); summary(model)
exp(cbind("Odds ratio" = coef(model), confint.default(model, level = 0.95)))

temp = temp[temp$PathStage != '', ]
ggscatter(temp, 'V1', 'V2', color = 'TripNeg', size = 4)
ggscatter(temp, 'V1', 'V3', color = 'TripNeg', size = 4)
ggscatter(temp, 'V2', 'V3', color = 'TripNeg', size = 4)

scatterplot3d(data$V3, data$V1, data$V2)
