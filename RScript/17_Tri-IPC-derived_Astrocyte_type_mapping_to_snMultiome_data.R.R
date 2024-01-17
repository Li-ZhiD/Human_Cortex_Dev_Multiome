library(Seurat)
library(Signac)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(singleCellNet)
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 300000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/20221216_GPC_differentiation")

#import GPC data and primary data
GPC <- readRDS("results/final_GPC.rds")
GPC_AST <- subset(GPC, subset = Sample_ID %in% c("GW22_IPC-Glia_DIV7", "GW22_IPC-Glia_DIV14") & type == "Astrocyte")

primary_data <-  readRDS("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/multiome_20230730_filtered_annotated.rds")

common_features <- intersect(rownames(GPC_AST), rownames(primary_data))
GPC_AST <- GPC_AST[common_features,]

#########################
#primary data at Infancy stage pre-processing
#########################
primary_AST_Infancy <- subset(primary_data, subset = (Group == "Infancy") & (subclass == "Astrocyte") )
DefaultAssay(primary_AST_Infancy) <- "integrated"
primary_AST_Infancy[["ATAC"]] <- NULL
primary_AST_Infancy <- primary_AST_Infancy[common_features,]

###########################################
#re-annotate primary data
primary_AST_Infancy <- RunUMAP(primary_AST_Infancy, reduction = "pca", dims = 1:50, return.model = TRUE)

#check markers for astrocyte subtypes
#key markers: GFAP, IGFBP5, SPARC for S100a11 lineage and OLIG1, CHRDL1, BTBD17
primary_AST_Infancy.markers.plot.list <- lapply(X = c("OLIG1", "CHRDL1", "BTBD17", "GFAP", "IGFBP5", "SPARC"), FUN = function(x) {
  FeaturePlot_scCustom(primary_AST_Infancy, reduction = "umap", features = x, raster = F)
})
patchwork::wrap_plots(primary_AST_Infancy.markers.plot.list, nrow = 2)
ggsave(filename = "results/integration_with_primary_data/Infancy/primary_AST_Infancy.markers.plot.list.png", width = 15, height = 9)

#re-cluster and annotate the primary dataset
DefaultAssay(primary_AST_Infancy) <-  "integrated"
primary_AST_Infancy <- FindNeighbors(primary_AST_Infancy, dims = 1:50)
primary_AST_Infancy <- FindClusters(primary_AST_Infancy, resolution = 0.8)

#annotate the lineages of individual clusters based on Zhou et al., 2023
annotation <- read.csv("results/integration_with_primary_data/Infancy/annotation.csv")
annotation$integrated_snn_res.0.8 <-  as.factor(annotation$integrated_snn_res.0.8)
metadata <- primary_AST_Infancy[[]]
metadata_annotated <- left_join(metadata, annotation, by = "integrated_snn_res.0.8")
rownames(metadata_annotated) <- rownames(metadata)
primary_AST_Infancy <- AddMetaData(primary_AST_Infancy, metadata_annotated)

#save data
saveRDS(primary_AST_Infancy, "results/integration_with_primary_data/Infancy/primary_AST_Infancy.rds")


########################
#Seurat label transfer
########################
#query dataset pre-processing
#SCTransform individual GPC_IPC_IAO datasets to be compatible with the label transfer, remove cell cycle effect in the same way as the reference dataset
GPC_AST.list <- SplitObject(GPC_AST, split.by = "Sample_ID")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
GPC_AST.list <- lapply(X = GPC_AST.list, FUN = CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes)
GPC_AST.list <- lapply(X = GPC_AST.list, FUN = SCTransform, vst.flavor = "v2", vars.to.regress = c("S.Score", "G2M.Score"))

anchors <- list()
for (i in 1:length(GPC_AST.list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = primary_AST_Infancy,
    query = GPC_AST.list[[i]],
    reference.reduction = "pca",
    dims = 1:50,
  )
}

for (i in 1:length(GPC_AST.list)) {
  GPC_AST.list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = GPC_AST.list[[i]],
    reference = primary_AST_Infancy, 
    refdata = list(lineage = "lineage"),
    reference.reduction = "pca",
    reduction.model = "umap"
  )
}

# Merge the two GPC dataset 
GPC_AST_Infancy <- merge(GPC_AST.list[[1]], GPC_AST.list[2:length(GPC_AST.list)], merge.dr = "ref.umap")
GPC_AST_Infancy$predicted.lineage.score.filtered <- ifelse(GPC_AST_Infancy$predicted.lineage.score > 0.5, GPC_AST_Infancy$predicted.lineage, "Unmapped")
GPC_AST_Infancy_lineage_plot <- DimPlot_scCustom(GPC_AST_Infancy, group.by = "predicted.lineage.score.filtered", reduction = "ref.umap", raster = F, colors_use = c("#e63946", "#1d3557"))
Primary_lineage_plot <- DimPlot_scCustom(primary_AST_Infancy, reduction = "umap", group.by = "lineage", raster = F, colors_use = c("#e63946", "#1d3557"))
Primary_lineage_plot + GPC_AST_Infancy_lineage_plot
ggsave("results/integration_with_primary_data/Infancy/GPC_mapped.png", width = 14, height = 6, units = "in")

#save data
saveRDS(GPC_AST_Infancy, "results/integration_with_primary_data/Infancy/GPC_mapped.rds")

########################
#SingleCellNet classification
########################
#import query data
GPC_AST_Infancy <- readRDS("results/integration_with_primary_data/Infancy/GPC_mapped.rds")
DefaultAssay(GPC_AST_Infancy) <- "RNA"
GPC_AST_Infancy_SCN <- extractSeurat(GPC_AST_Infancy, exp_slot_name = "counts")
GPC_AST_Infancy_tab <-  GPC_AST_Infancy_SCN$sampTab
GPC_AST_Infancy_tab <- droplevels(GPC_AST_Infancy_tab)
GPC_AST_Infancy_tab$cell <- rownames(GPC_AST_Infancy_tab)
GPC_AST_Infancy_exp <-  GPC_AST_Infancy_SCN$expDat

#import training data
primary_AST_Infancy <- readRDS("results/integration_with_primary_data/Infancy/primary_AST_Infancy.rds")
DefaultAssay(primary_AST_Infancy) <- "RNA"
primary_AST_Infancy_SCN <- extractSeurat(primary_AST_Infancy, exp_slot_name = "counts")
primary_AST_Infancy_tab <-  primary_AST_Infancy_SCN$sampTab
primary_AST_Infancy_tab <- droplevels(primary_AST_Infancy_filt_tab)
primary_AST_Infancy_tab$cell <- rownames(primary_AST_Infancy_filt_tab)
primary_AST_Infancy_exp <-  primary_AST_Infancy_SCN$expDat

#Find genes in common to the data sets and remove genes with variance = 1
RowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
goodGenes_GPC_AST_Infancy <- rownames(GPC_AST_Infancy_exp)[RowVars(GPC_AST_Infancy_exp) != 0]
goodGenes_primary_AST_Infancy <- rownames(primary_AST_Infancy_exp)[RowVars(primary_AST_Infancy_exp) != 0]
commonGenes <-  intersect(goodGenes_GPC_AST_Infancy, goodGenes_primary_AST_Infancy)

GPC_AST_Infancy_exp <-  GPC_AST_Infancy_exp[commonGenes,]
primary_AST_Infancy_exp <- primary_AST_Infancy_exp[commonGenes,]

##########################
#train a classifier for lineage
#Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=primary_AST_Infancy_tab, ncells = 698, dLevel="lineage")
stTrain = stList[[1]]
expTrain = primary_AST_Infancy_exp[,rownames(stTrain)]


#Train the classifier
#If you increase nTopGenes and nTopGenePairs, you may get a even better classifier performance on query data!
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 200, nRand = 400, nTrees = 1000, nTopGenePairs = 200, dLevel = "lineage", colName_samp = "cell"))

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=698, dLevel="lineage") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = primary_AST_Infancy_exp[commonGenes, rownames(stTest)]

#predict testing data
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 100)
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "lineage", classQuery = "lineage", nRand = 100)
pdf("results/integration_with_primary_data/singleCellNet/tm_heldoutassessment_PR_type.pdf")
plot_PRs(tm_heldoutassessment)
dev.off()
pdf("results/integration_with_primary_data/singleCellNet/tm_heldoutassessment_metrics_type.pdf")
plot_metrics(tm_heldoutassessment)
dev.off()

tm_assess_matrix = tm_heldoutassessment$nonNA_PR
score = 0.35
celltype = "Olig2"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#"SCN score of 0.35 for cell type Olig2 has precision of 0.929 ~ 0.931 and sensitivity of 0.96 ~ 0.96"
score = 0.35
celltype = "S100a11"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#"SCN score of 0.35 for cell type S100a11 has precision of 0.939 ~ 0.942 and sensitivity of 0.955 ~ 0.954"


#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 100
sla = as.vector(stTest$lineage)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created
#heatmap
pdf("results/integration_with_primary_data/singleCellNet/evaluate_hmClass_type.pdf")
sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
dev.off()
#bar graph
pdf("results/integration_with_primary_data/singleCellNet/evaluate_plot_attr_type.pdf")
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="lineage", sid="cell")
dev.off()

#apply to query data
nqRand = 100
system.time(classRes_GPC_AST_Infancy <- scn_predict(class_info[['cnProc']], GPC_AST_Infancy_exp, nrand=nqRand))

# heatmap classification result
sgrp = as.vector(GPC_AST_Infancy_tab$type)
names(sgrp) = as.vector(GPC_AST_Infancy_tab$type)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

pdf("results/integration_with_primary_data/singleCellNet/query_plot_attr_type.pdf")
plot_attr(classRes_GPC_AST_Infancy, GPC_AST_Infancy_tab, nrand=nqRand, sid="cell", dLevel="type")
dev.off()

# This classifies a cell with  the category with the highest classification score or higher than a classification score threshold of your choosing (we choose 0.35 here).
GPC_AST_Infancy_tab2 <- get_cate(classRes = classRes_GPC_AST_Infancy, sampTab = GPC_AST_Infancy_tab, dLevel = "type", sid = "cell", nrand = nqRand, cThresh = 0.35)
# Olig2    rand S100a11 
# 971     328    1107

#save results
save(class_info, classRes_val_all, classRes_GPC_AST_Infancy, file = "results/integration_with_primary_data/singleCellNet/scn_models_type.RData")

