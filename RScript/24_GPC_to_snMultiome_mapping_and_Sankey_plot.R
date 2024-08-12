library(Seurat)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(singleCellNet)
library(networkD3)
library(viridis)
library(hrbrthemes)
library(circlize)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/")

#import primary data
dat <- readRDS("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/multiome_20230730_filtered_annotated.rds")
GPC <- readRDS("20221216_GPC_differentiation/results/final_GPC.rds")

#######################
#SingleCellNet classification
#######################
#extract information from primary data
DefaultAssay(dat) <- "RNA"
dat_for_scn <- extractSeurat(dat, exp_slot_name = "counts")
dat_tab <-  dat_for_scn$sampTab
dat_tab <- droplevels(dat_tab)
dat_tab$cell <- rownames(dat_tab)
dat_exp <-  dat_for_scn$expDat
rm(dat) #to save memory

#extract data from GBM Seurat object
GPC_for_scn <- extractSeurat(GPC, exp_slot_name = "counts")
GPC_tab <-  GPC_for_scn$sampTab
GPC_tab <- droplevels(GPC_tab)
GPC_tab$cell <- rownames(GPC_tab)
GPC_exp <-  GPC_for_scn$expDat


#Find genes in common to the data sets and limit analysis to these genes, remove gene with variance = 0
RowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
goodGenes <- rownames(GPC_exp)[RowVars(dat_exp) != 0]
commonGenes <-  intersect(goodGenes, rownames(dat_exp))

dat_exp <-  dat_exp[commonGenes,]
GPC_exp <- GPC_exp[commonGenes,]


##########################
#train a classifier for types from primary data
#Split for training and assessment, and transform training data
set.seed(1) #can be any random seed number
stList = splitCommon(sampTab=dat_tab, ncells=700, dLevel="type")
stTrain = stList[[1]]
saveRDS(stTrain, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/stTrain.rds")
expTrain = dat_exp[,rownames(stTrain)]
saveRDS(expTrain, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/expTrain.rds")

#Train the classifier
#If you increase nTopGenes and nTopGenePairs, you may get a even better classifier performance on query data!
#if nTopGenes = 75 and nTopGenePairs = 125, vector too long for GPC prediction, leading to "Error in predict.randomForest(rfObj, t(expQuery[preds, ]), type = "prob") : long vectors (argument 7) are not supported in .C"
#The performance of the model is only slightly worse with fewer gene pairs.
system.time(class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 60, nRand = 400, nTrees = 1000, nTopGenePairs = 150, dLevel = "type", colName_samp = "cell"))
saveRDS(class_info, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/class_info.rds")

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=165, dLevel="type") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
saveRDS(stTest, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/stTest.rds")
expTest = dat_exp[commonGenes,rownames(stTest)]

#predict testing data
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 165)
saveRDS(classRes_val_all, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/classRes_val_all.rds")

tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "type", classQuery = "type", nRand = 165)
#average AUPRC = 0.827
pdf("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/tm_heldoutassessment_PR_type.pdf")
plot_PRs(tm_heldoutassessment)
dev.off()
pdf("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/tm_heldoutassessment_metrics_type.pdf")
plot_metrics(tm_heldoutassessment)
dev.off()

tm_assess_matrix = tm_heldoutassessment$nonNA_PR
score = 0.2
celltype = "IN-MGE-Immature"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)


#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 165
sla = as.vector(stTest$type)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created
#heatmap
pdf("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/evaluate_hmClass_type.pdf")
sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
dev.off()
#bar graph
pdf("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/evaluate_plot_attr_type.pdf")
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="type", sid="cell")
dev.off()

#apply to query data
nqRand = 400
system.time(classRes_GPC <- scn_predict(class_info[['cnProc']], GPC_exp, nrand=nqRand))
saveRDS(classRes_GPC, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/classRes_GPC.rds")

# heatmap classification result
sgrp = as.vector(GPC_tab$type)
names(sgrp) = as.vector(GPC_tab$cell)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)
pdf("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/query_hmClass_type.pdf")
sc_hmClass(classRes_GPC, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
dev.off()

pdf("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/query_plot_attr_type.pdf")
plot_attr(classRes_GPC, GPC_tab, nrand=nqRand, sid="cell", dLevel="type")
dev.off()

# This classifies a cell with  the category with the highest classification score or higher than a classification score threshold of your choosing (we choose 0.2 here).
GPC_tab2 <- get_cate(classRes = classRes_GPC, sampTab = GPC_tab, dLevel = "type", sid = "cell", nrand = nqRand, cThresh = 0)
GPC_tab3 <- get_cate(classRes = classRes_GPC, sampTab = GPC_tab, dLevel = "type", sid = "cell", nrand = nqRand, cThresh = 0.1)
GPC_tab4 <- get_cate(classRes = classRes_GPC, sampTab = GPC_tab, dLevel = "type", sid = "cell", nrand = nqRand, cThresh = 0.15)
GPC_tab5 <- get_cate(classRes = classRes_GPC, sampTab = GPC_tab, dLevel = "type", sid = "cell", nrand = nqRand, cThresh = 0.2)
sc_violinClass(sampTab = GPC_tab5, classRes = classRes_GPC, sid = "cell", dLevel = "type", addRand = nqRand)
ggsave("20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/sc_violinClass.png", width = 6, height = 30, bg = "white")
GPC$scn_type_category_thresh0 <- GPC_tab2$category
GPC$scn_type_category_thresh0 <- ifelse(GPC$scn_type_category_thresh0 == "rand", "Unmapped", GPC$scn_type_category_thresh0)
GPC$scn_type_category_thresh0.1 <- GPC_tab3$category
GPC$scn_type_category_thresh0.1 <- ifelse(GPC$scn_type_category_thresh0.1 == "rand", "Unmapped", GPC$scn_type_category_thresh0.1)
GPC$scn_type_category_thresh0.15 <- GPC_tab4$category
GPC$scn_type_category_thresh0.15 <- ifelse(GPC$scn_type_category_thresh0.15 == "rand", "Unmapped", GPC$scn_type_category_thresh0.15)
GPC$scn_type_category_thresh0.2 <- GPC_tab5$category
GPC$scn_type_category_thresh0.2 <- ifelse(GPC$scn_type_category_thresh0.2 == "rand", "Unmapped", GPC$scn_type_category_thresh0.2)
GPC$scn_type_score <- GPC_tab2$scn_score

mapping_result <- as.data.frame.matrix(table(GPC$scn_type_category_thresh0.2, GPC$type))
write.csv(mapping_result, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/mapping_result.csv")
#save the seurat object with prediction results
saveRDS(GPC, "20221216_GPC_differentiation/results/All_GPC_to_snMultiome_SingleCellNet/GPC_Seurat_with_scn_prediction.rds")


#save metadata
metadata <- GPC[[]]
metadata$Cell_ID <- rownames(metadata)
metadata_filt <- metadata[,c("Cell_ID", "Sample_ID", "Seeding_cell_type", "Stage", "nCount_RNA", "nFeature_RNA", "percent.mito", "class", "subclass", "type", "RNA_snn_res.6", "scn_type_category_thresh0.2")]
colnames(metadata_filt) <- c("Cell_ID", "Sample_ID", "Seeding_cell_type", "Stage", "nCount_RNA", "nFeature_RNA", "Percentage_in_mitochondrial_genes", "Class", "Subclass", "Type", "Cluster", "SingleCellNet_predictions")
write.csv(metadata_filt, "results/cell_level_metadata_v2.csv", row.names = F)



##################################
#draw Sankey plot
##################################
# Load dataset
data <- read.csv("20221216_GPC_differentiation_scRNA-seq/results/All_GPC_to_snMultiome_SingleCellNet/mapping_result.csv", row.names = 1)

# convert long format
data_long <- data %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)
colnames(data_long) <- c("target", "source", "value")
data_long$source <- ifelse(data_long$source == "Ependymal.cell", "Ependymal cell", data_long$source)
data_long$source <- ifelse(data_long$source == "IPC.EN", "IPC-EN", data_long$source)
data_long$source <- ifelse(data_long$source == "IPC.Glia", "IPC-Glia(Tri-IPC)", data_long$source)
data_long$source <- ifelse(data_long$source == "IPC.IN", "IPC-IN", data_long$source)
#distinguish target and source with the same name
data_long$target <- paste(data_long$target, " ", sep="")

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(as.character(data_long$source), as.character(data_long$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1

# prepare colour scale
ColourScal <- 'd3.scaleOrdinal().domain(["Dividing", "RG", "Ependymal cell", "IPC-EN", "EN", "IPC-Glia(Tri-IPC)", "Astrocyte", "OPC", "IPC-IN", "IN", "Astrocyte-Fibrous", "Astrocyte-Immature", "IN-dLGE-Immature", "IN-MGE-Immature", "IPC-Glia", "RG-oRG", "RG-tRG", "RG-vRG", "Unmapped", "Vascular", "EN-Newborn", "Unknown", "Oligodendrocyte-Immature", "IN-CGE-Immature", "IN-CGE-VIP"]).range(["#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A", "#BDCDFF", "#782AB6", "#1C8356", "#FA0087", "#C075A6", "#FE00FA", "#F6222E", "#5A5156", "#d3d3d3", "#3B00FB", "#3283FE", "#E4E1E3", "#7ED7D1", "#85660D", "#B10DA1"])'

# Make the Network
sankey <- sankeyNetwork(Links = data_long, Nodes = nodes,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight=FALSE,
                        colourScale=ColourScal,
                        nodeWidth=40, fontSize=20, nodePadding=20)

saveNetwork(sankey, file = "20221216_GPC_differentiation_scRNA-seq/results/All_GPC_to_snMultiome_SingleCellNet/sankey_plot.html")