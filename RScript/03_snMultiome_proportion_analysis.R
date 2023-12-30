library(Seurat)
library(Signac)
library(tidyverse)
library(scCustomize)
library(speckle)
library(limma)


setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")

dat <- readRDS("multiome_20230730_filtered_annotated.rds")
dat$log2_age <- log2(dat$Estimated_postconceptional_age_in_days)

#################
#we will first focus on finding cell type proportions that change with age
#transform prop using logit
props <- getTransformedProps(dat$type, dat$dataset, transform="logit")
transformed <- props$TransformedProps
transformed_sample <- data.frame(colnames(transformed))
colnames(transformed_sample) <- "dataset"
#add block information to account for samples come from the same individuals
transformed_sample$block <- c(21,22,23,24,25,26,27,1,1,2,2,3,4,5,6,7,8,9,10,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,20,20)
#add metadata in the same sequence of transformed data
metadata <- data.frame(dat$dataset, dat$region_summary, dat$log2_age)
metadata <- distinct(metadata)
rownames(metadata) <- NULL
colnames(metadata) <- c("dataset", "region", "age")
metadata2 <- transformed_sample %>% left_join(metadata)
#limma
design <- model.matrix(~metadata2$region + metadata2$age)
colnames(design) <- c("Intercept", "PFC", "V1", "age")
#estimate correlation
dupcor <- duplicateCorrelation(transformed, design=design,
                               block=metadata2$block)
#fit model
fit1 <- lmFit(transformed, design=design, block=metadata2$block, 
              correlation=dupcor$consensus.correlation)
fit1 <- eBayes(fit1, robust=TRUE)
#age effect
age_results <- topTable(fit1, number=100, coef = 4, adjust="BH")
#save results
write.csv(age_results, file = "results/proportion/age_results.csv", row.names = T)


##################
#region effect, because this is dependent on age groups and is not a simple interaction,
#we will assess region effect at or after the third trimester when regional differences become more prominent
#subset samples to focus on third trimester or later
dat_sub <- subset(dat, subset = Group %in% c("Third_trimester", "Infancy", "Adolescence"))
dat_sub$type <- droplevels(dat_sub$type)
dat_sub$dataset <- droplevels(dat_sub$dataset)
dat_sub$region_summary <- droplevels(dat_sub$region_summary)

props <- getTransformedProps(dat_sub$type, dat_sub$dataset, transform="logit")
transformed <- props$TransformedProps
transformed_sample <- data.frame(colnames(transformed))
colnames(transformed_sample) <- "dataset" 
#add block information to account for some samples come from the same individuals
transformed_sample$block <- c(1,2,3,4,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14)
#add metadata in the same sequence of transformed data
metadata <- data.frame(dat_sub$dataset, dat_sub$region_summary, dat_sub$log2_age)
metadata <- distinct(metadata)
rownames(metadata) <- NULL
colnames(metadata) <- c("dataset", "region", "age")
metadata2 <- transformed_sample %>% left_join(metadata)
#limma
design <- model.matrix(~metadata2$region + metadata2$age)
colnames(design) <- c("Intercept", "V1", "age")

#estimate correlation
dupcor <- duplicateCorrelation(transformed, design=design,
                               block=metadata2$block)
#fit model
fit2 <- lmFit(transformed, design=design, block=metadata2$block, 
              correlation=dupcor$consensus.correlation)
fit2 <- eBayes(fit2, robust=TRUE)

#region effect
region_results <- topTable(fit2, number=100, coef = 2, adjust="BH")
#save results
write.csv(region_results, file = "results/proportion/region_results.csv", row.names = T)



