library(Seurat)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(singleCellNet)
library(future)


options(warn = 1)

plan("multicore", workers = 12)
options(future.globals.maxSize = 256000 * 1024^2)


setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/20221216_GPC_differentiation")
#import GPC data and Shi2021_GE data
experiment.aggregate.final <- readRDS("results/final_GPC.rds")
Shi2021_GE <- readRDS("Shi2021/Shi2021_annotated.rds")
DefaultAssay(Shi2021_GE) <- "integrated"
Shi2021_GE <- RunUMAP(Shi2021_GE, dims = 1:30, reduction = "pca", return.model = TRUE)

#save Shi2021 with UMAP model
saveRDS(Shi2021_GE, "Shi2021/Shi2021_annotated.rds")

#######################
#Seurat label transfer
#######################
#label transfer for inhibitory lineage
GPC_IN <- subset(experiment.aggregate.final, subset = type %in% c("IN"))

#label transfer
anchors <- FindTransferAnchors(reference = Shi2021_GE, query = GPC_IN,
                               dims = 1:30, reference.reduction = "pca")

#UMAP Projection
type_mapping <- MapQuery(anchorset = anchors, reference = Shi2021_GE, query = GPC_IN,
                         refdata = list(type = "type"), reference.reduction = "pca", reduction.model = "umap")
type_mapping$predicted.type.filtered <- ifelse(type_mapping$predicted.type.score > 0.5, type_mapping$predicted.type, "Unknown")

predicted.terminal_type <- TransferData(anchorset = anchors, refdata = Shi2021_GE$inferred_terminal_type, dims = 1:30)
type_mapping$predicted.terminal_type <- predicted.terminal_type$predicted.id
type_mapping$predicted.terminal_type.score <- predicted.terminal_type$prediction.score.max
type_mapping$predicted.terminal_type.filtered <- ifelse(type_mapping$predicted.terminal_type.score > 0.5, type_mapping$predicted.terminal_type, "Unknown")

predicted.origin <- TransferData(anchorset = anchors, refdata = Shi2021_GE$origin, dims = 1:30)
type_mapping$predicted.origin <- predicted.origin$predicted.id
type_mapping$predicted.origin.score <- predicted.origin$prediction.score.max
type_mapping$predicted.origin.filtered <- ifelse(type_mapping$predicted.origin.score > 0.5, type_mapping$predicted.origin, "Unknown")

#plot results
type <- DimPlot_scCustom(Shi2021_GE, group.by = "type", reduction = "umap", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", "#3283FEFF", 
                                                                                           "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", "#2ED9FFFF", 
                                                                                           "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", "#C4451CFF", 
                                                                                           "#1C8356FF", "#85660DFF", "#B10DA1FF"), pt.size = 0.1) + xlim (c(-15, 10)) + ylim (c(-15, 10))
p_type <- DimPlot_scCustom(type_mapping, reduction = "ref.umap", group.by = "predicted.type.filtered", colors_use = c("#325A9BFF", "#85660DFF", "#B10DA1FF", "#AA0DFEFF", "#F8A19FFF" , "#FEAF16FF", "#B00068FF", "#1C8356FF", "lightgrey"), pt.size = 0.1) + xlim (c(-15, 10)) + ylim (c(-15, 10))
type + p_type
ggsave(file="results/integration_with_Shi2021/type_comparison_umap.png", width = 20, height = 7, units = "in")

terminal_type <- DimPlot_scCustom(Shi2021_GE, reduction = "umap", group.by = "inferred_terminal_type", colors_use = c("#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#2ED9FFFF", 
                                                                                                                      "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#C4451CFF"), na.value="lightgrey", pt.size = 0.1) + xlim (c(-15, 10)) + ylim (c(-15, 10))
p_terminal_type <- DimPlot_scCustom(type_mapping, reduction = "ref.umap", group.by = "predicted.terminal_type.filtered", colors_use = c("#B00068FF", "#F8A19FFF", "#AA0DFEFF", "#FEAF16FF", "#3283FEFF", "lightgrey"), pt.size = 0.1) + xlim (c(-15, 10)) + ylim (c(-15, 10))
terminal_type + p_terminal_type
ggsave(file="results/integration_with_Shi2021/terminal_type_comparison_umap.png",width = 25, height = 8, units = "in")

#combine all four plots
type2 <- type + NoLegend() + labs(title = NULL)
p_type2 <- p_type + NoLegend() + labs(title = NULL)
terminal_type2 <- terminal_type + NoLegend() + labs(title = NULL)
p_terminal_type2 <- p_terminal_type + NoLegend() + labs(title = NULL)
type2 + terminal_type2 + p_type2 + p_terminal_type2 + plot_layout(ncol = 2)
ggsave(file="results/integration_with_Shi2021/comparison_umap.png",width = 8, height = 8, units = "in")

#save results
saveRDS(type_mapping, "results/integration_with_Shi2021/GPC_IN_type_mapping.rds")

#######################
#SingleCellNet classification
#######################
#extract data from Seurat object
GPC_IN <- extractSeurat(GPC_IN, exp_slot_name = "counts")
GPC_tab <-  GPC_IN$sampTab
GPC_tab <- droplevels(GPC_tab)
GPC_tab$cell <- rownames(GPC_tab)
GPC_exp <-  GPC_IN$expDat

#subset the training data to focus on IN only
DefaultAssay(Shi2021_GE) <- "RNA"
Shi2021_IN <- subset(Shi2021_GE, subset = (subclass == "IN" & origin != "GE"))
Shi2021_IN <- extractSeurat(Shi2021_IN, exp_slot_name = "counts")
Shi2021_tab <-  Shi2021_IN$sampTab
Shi2021_tab <- droplevels(Shi2021_tab)
Shi2021_tab$cell <- rownames(Shi2021_tab)
Shi2021_exp <-  Shi2021_IN$expDat

#Find genes in common to the data sets and limit analysis to these genes, remove gene with variance = 0
RowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
goodGenes <- rownames(GPC_exp)[RowVars(Shi2021_exp) != 0]
commonGenes <-  intersect(goodGenes, rownames(Shi2021_exp))

Shi2021_exp <-  Shi2021_exp[commonGenes,]
GPC_exp <- GPC_exp[commonGenes,]

##########################
#train a classifier for types
#Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=Shi2021_tab, ncells=400, dLevel="type")
stTrain = stList[[1]]
expTrain = Shi2021_exp[,rownames(stTrain)]

#Train the classifier
#If you increase nTopGenes and nTopGenePairs, you may get a even better classifier performance on query data!
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 200, nRand = 400, nTrees = 1000, nTopGenePairs = 200, dLevel = "type", colName_samp = "cell"))

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="type") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = Shi2021_exp[commonGenes,rownames(stTest)]

#predict testing data
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 100)
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "type", classQuery = "type", nRand = 100)
pdf("results/singleCellNet/tm_heldoutassessment_PR_type.pdf")
plot_PRs(tm_heldoutassessment)
dev.off()
pdf("results/singleCellNet/tm_heldoutassessment_metrics_type.pdf")
plot_metrics(tm_heldoutassessment)
dev.off()

tm_assess_matrix = tm_heldoutassessment$nonNA_PR
score = 0.35
celltype = "LGE/CGE_MEIS2/PAX6"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#"SCN score of 0.35 for cell type CGE/LGE_MEIS2/PAX6 has precision of 0.759 ~ 0.756 and sensitivity of 0.63 ~ 0.62"

#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 100
sla = as.vector(stTest$type)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created
#heatmap
pdf("results/singleCellNet/evaluate_hmClass_type.pdf")
sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
dev.off()
#bar graph
pdf("results/singleCellNet/evaluate_plot_attr_type.pdf")
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="type", sid="cell")
dev.off()

#apply to query data
nqRand = 100
system.time(classRes_GPC_all<-scn_predict(class_info[['cnProc']], GPC_exp, nrand=nqRand))

# heatmap classification result
sgrp = as.vector(GPC_tab$subclass)
names(sgrp) = as.vector(GPC_tab$cell)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)
pdf("results/singleCellNet/query_hmClass_type.pdf")
sc_hmClass(classRes_GPC_all, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
dev.off()

pdf("results/singleCellNet/query_plot_attr_type.pdf")
plot_attr(classRes_GPC_all, GPC_tab, nrand=nqRand, sid="cell", dLevel="subclass")
dev.off()

# This classifies a cell with  the category with the highest classification score or higher than a classification score threshold of your choosing (we choose 0.35 here).
GPC_tab3 <- get_cate(classRes = classRes_GPC_all, sampTab = GPC_tab, dLevel = "type", sid = "cell", nrand = nqRand, cThresh = 0.35)
GPC_IN <- subset(experiment.aggregate.final, subset = subclass == "IN")
GPC_IN$scn_type_category <- GPC_tab3$category
GPC_IN$scn_type_category <- ifelse(GPC_IN$scn_type_category == "rand", "Unknown", GPC_IN$scn_type_category)
GPC_IN$scn_type_score <- GPC_tab3$scn_score

#save results
save(class_info, classRes_val_all, classRes_GPC_all, file = "results/singleCellNet/scn_models_type.RData")

#save the seurat object with prediction results
saveRDS(GPC_IN, "results/singleCellNet/GPC_IN_Seurat_with_scn_prediction.rds")

########################################
#plot proportions
metadata <- GPC_IN[[]]
pt_scn_type <- table(metadata$subclass, metadata$scn_type_category)
pt_scn_type <- as.data.frame(pt_scn_type)
colnames(pt_scn_type) <- c("subclass", "scn_type", "Freq")
pt_scn_type_plot <- ggplot(pt_scn_type, aes(x = subclass, y = Freq, fill = scn_type)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Proportion") +
  xlab("subclass") +
  scale_fill_manual(values = DiscretePalette_scCustomize(num_colors = 36, palette = "polychrome")[c(1,2,12,13,14,20)]) +
  theme(legend.title = element_blank(), legend.position = "none") +
  coord_flip()
ggsave("results/singleCellNet/type_pt.png", width = 5, height = 1.4, units = "in")
