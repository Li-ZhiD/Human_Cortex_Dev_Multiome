library(Seurat)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(singleCellNet)
library(future)

setwd("/kriegsteinlab/data3/LiWang/analysis/GBM/GBmap")
#import GBM data
GBmap <- readRDS("GBmap_extended_malignant_cell.rds")

#change reduction object name
new_reduction <- CreateDimReducObject(embeddings = Embeddings(GBmap, reduction = "Umap"), key = "umap")
GBmap[["Umap"]] <- NULL
GBmap[["umap"]] <- new_reduction

#change metadata column name
metadata <- GBmap[[]]
meta_colnames <- dput(colnames(GBmap[[]]))
colnames(metadata) <- c("author", "donor_id", "assay_ontology_term_id", "cell_type_ontology_term_id", 
                        "development_stage_ontology_term_id", "disease_ontology_term_id", 
                        "self_reported_ethnicity_ontology_term_id", "is_primary_data", 
                        "organism_ontology_term_id", "sex_ontology_term_id", "annotation_level_1", 
                        "annotation_level_2", "annotation_level_3", "gbmap", "method", 
                        "stage", "location", "sector", "celltype_original", "EGFR_status", "MET_status", 
                        "p53_status", "TERT_status", "ATRX_status", "PTEN_status", "MGMT_status", "chr1p19q_status", "PDGFR_status", "suspension_type", 
                        "tissue_ontology_term_id", "tissue_type", "cell_type", "assay", 
                        "disease", "organism", "sex", "tissue", "self_reported_ethnicity", 
                        "development_stage")
GBmap <- AddMetaData(GBmap, metadata)
GBmap$EGFR <- NULL
GBmap$MET <- NULL
GBmap$p53 <- NULL
GBmap$TERT <- NULL
GBmap$ATRX <- NULL
GBmap$PTEN <- NULL
GBmap$MGMT <- NULL
GBmap$chr1p19q <- NULL
GBmap$PDGFR <- NULL

saveRDS(GBmap, "GBmap_extended_malignant_cell_clean.rds")

#import primary data
dat <- readRDS("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/multiome_20230730_filtered_annotated.rds")

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
GBmap_for_scn <- extractSeurat(GBmap, exp_slot_name = "counts")
GBmap_tab <-  GBmap_for_scn$sampTab
GBmap_tab <- droplevels(GBmap_tab)
GBmap_tab$cell <- rownames(GBmap_tab)
GBmap_exp <-  GBmap_for_scn$expDat


#Find genes in common to the data sets and limit analysis to these genes, remove gene with variance = 0
RowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
goodGenes <- rownames(GBmap_exp)[RowVars(GBmap_exp) != 0]
commonGenes <-  intersect(goodGenes, rownames(dat_exp))

dat_exp <-  dat_exp[commonGenes,]
GBmap_exp <- GBmap_exp[commonGenes,]

##########################
#train a classifier for types from primary data
#Split for training and assessment, and transform training data
set.seed(1) #can be any random seed number
stList = splitCommon(sampTab=dat_tab, ncells=700, dLevel="type")
stTrain = stList[[1]]
saveRDS(stTrain, "results/SingleCellNet/stTrain.rds")
expTrain = dat_exp[,rownames(stTrain)]
saveRDS(expTrain, "results/SingleCellNet/expTrain.rds")

#Train the classifier
#If you increase nTopGenes and nTopGenePairs, you may get a even better classifier performance on query data!
#if nTopGenes = 75 and nTopGenePairs = 125, vector too long for GBmap prediction, leading to "Error in predict.randomForest(rfObj, t(expQuery[preds, ]), type = "prob") : long vectors (argument 7) are not supported in .C"
#The performance of the model is only slightly worse with fewer gene pairs.
system.time(class_info <- scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 60, nRand = 400, nTrees = 1000, nTopGenePairs = 150, dLevel = "type", colName_samp = "cell"))
saveRDS(class_info, "results/SingleCellNet/class_info.rds")

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=165, dLevel="type") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
saveRDS(stTest, "results/SingleCellNet/stTest.rds")
expTest = dat_exp[commonGenes,rownames(stTest)]

#predict testing data
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 165)
saveRDS(classRes_val_all, "results/SingleCellNet/classRes_val_all.rds")

tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "type", classQuery = "type", nRand = 165)
#average AUPRC = 0.832
pdf("results/SingleCellNet/tm_heldoutassessment_PR_type.pdf")
plot_PRs(tm_heldoutassessment)
dev.off()
pdf("results/SingleCellNet/tm_heldoutassessment_metrics_type.pdf")
plot_metrics(tm_heldoutassessment)
dev.off()

tm_assess_matrix = tm_heldoutassessment$nonNA_PR
score = 0.15
celltype = "IPC-Glia"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#"SCN score of 0.15 for cell type IPC-Glia has precision of 0.639 ~ 0.658 and sensitivity of 0.933 ~ 0.933"

#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 165
sla = as.vector(stTest$type)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created
#heatmap
pdf("results/SingleCellNet/evaluate_hmClass_type.pdf")
sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
dev.off()
#bar graph
pdf("results/SingleCellNet/evaluate_plot_attr_type.pdf")
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="type", sid="cell")
dev.off()

#apply to query data
nqRand = 400
system.time(classRes_GBmap <- scn_predict(class_info[['cnProc']], GBmap_exp, nrand=nqRand))
saveRDS(classRes_GBmap, "results/SingleCellNet/classRes_GBmap.rds")

# heatmap classification result
sgrp = as.vector(GBmap_tab$annotation_level_3)
names(sgrp) = as.vector(GBmap_tab$cell)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)
pdf("results/SingleCellNet/query_hmClass_type.pdf")
sc_hmClass(classRes_GBmap, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
dev.off()

pdf("results/SingleCellNet/query_plot_attr_type.pdf")
plot_attr(classRes_GBmap, GBmap_tab, nrand=nqRand, sid="cell", dLevel="annotation_level_3")
dev.off()

# This classifies a cell with  the category with the highest classification score or higher than a classification score threshold of your choosing (we choose 0.15 here).
GBmap_tab2 <- get_cate(classRes = classRes_GBmap, sampTab = GBmap_tab, dLevel = "annotation_level_3", sid = "cell", nrand = nqRand, cThresh = 0)
sc_violinClass(sampTab = GBmap_tab2, classRes = classRes_GBmap, sid = "cell", dLevel = "annotation_level_3", addRand = nqRand)
ggsave("results/SingleCellNet/sc_violinClass.png", width = 6, height = 30, bg = "white")
GBmap_tab3 <- get_cate(classRes = classRes_GBmap, sampTab = GBmap_tab, dLevel = "annotation_level_3", sid = "cell", nrand = nqRand, cThresh = 0.1)
GBmap_tab4 <- get_cate(classRes = classRes_GBmap, sampTab = GBmap_tab, dLevel = "annotation_level_3", sid = "cell", nrand = nqRand, cThresh = 0.15)
GBmap_tab5 <- get_cate(classRes = classRes_GBmap, sampTab = GBmap_tab, dLevel = "annotation_level_3", sid = "cell", nrand = nqRand, cThresh = 0.2)
GBmap_tab6 <- get_cate(classRes = classRes_GBmap, sampTab = GBmap_tab, dLevel = "annotation_level_3", sid = "cell", nrand = nqRand, cThresh = 0.3)
GBmap$scn_type_category_thresh0 <- GBmap_tab2$category
GBmap$scn_type_category_thresh0 <- ifelse(GBmap$scn_type_category_thresh0 == "rand", "Unmapped", GBmap$scn_type_category_thresh0)
GBmap$scn_type_category_thresh0.1 <- GBmap_tab3$category
GBmap$scn_type_category_thresh0.1 <- ifelse(GBmap$scn_type_category_thresh0.1 == "rand", "Unmapped", GBmap$scn_type_category_thresh0.1)
GBmap$scn_type_category_thresh0.15 <- GBmap_tab4$category
GBmap$scn_type_category_thresh0.15 <- ifelse(GBmap$scn_type_category_thresh0.15 == "rand", "Unmapped", GBmap$scn_type_category_thresh0.15)
GBmap$scn_type_category_thresh0.2 <- GBmap_tab5$category
GBmap$scn_type_category_thresh0.2 <- ifelse(GBmap$scn_type_category_thresh0.2 == "rand", "Unmapped", GBmap$scn_type_category_thresh0.2)
GBmap$scn_type_category_thresh0.3 <- GBmap_tab6$category
GBmap$scn_type_category_thresh0.3 <- ifelse(GBmap$scn_type_category_thresh0.3 == "rand", "Unmapped", GBmap$scn_type_category_thresh0.3)
GBmap$scn_type_score <- GBmap_tab2$scn_score

#save the seurat object with prediction results
saveRDS(GBmap, "results/SingleCellNet/GBmap_Seurat_with_scn_prediction.rds")
#export cell ID of those Tri-IPC-like cells
Tri_IPC_like <- GBmap[, GBmap$scn_type_category_thresh0.15 == "IPC-Glia"]
Tri_IPC_like_ID <- Cells(Tri_IPC_like)
write.csv(Tri_IPC_like_ID, "Tri_IPC_like_ID.csv", row.names = FALSE)
mean(as.data.frame.matrix(table(GBmap$donor_id, GBmap$scn_type_category_thresh0.15))$`IPC-Glia` != 0)
#0.8654709 samples contain Tri-IPCs
########################################
#plot results
#choose 0.15 as threshold
#set consistent colors
type_colors <-  c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", "#3283FEFF",
                  "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", "#2ED9FFFF",
                  "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", "#C4451CFF",
                  "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", "#1CBE4FFF",
                  "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", "#782AB6FF",
                  "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF", "#1C7F93FF",
                  "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#822E1CFF", "#E4E1E3FF")
names(type_colors) <- c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                        "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                        "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                        "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                        "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown", "Unmapped")
GBmap$scn_type_category_thresh0.15 <- factor(GBmap$scn_type_category_thresh0.15, levels = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                                                                            "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                                                                                            "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                                                                            "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                                                                            "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown", "Unmapped"))
GBmap$scn_type_category_thresh0.15 <- droplevels(GBmap$scn_type_category_thresh0.15)
#plot distribution of cell types
state_plot <- DimPlot_scCustom(GBmap, group.by = "annotation_level_3", reduction = "umap", pt.size = 1, raster = T) + labs(title = NULL)
type_plot <- DimPlot_scCustom(GBmap, group.by = "scn_type_category_thresh0.15", reduction = "umap", colors_use = type_colors, pt.size = 1, raster = T) + labs(title = NULL)
state_plot + type_plot
ggsave(file="results/SingleCellNet/type_and_state_umap.png",width = 16.5, height = 5, units = "in", dpi = 600)
ggsave(file="results/SingleCellNet/type_and_state_umap.pdf",width = 16.5, height = 5, units = "in")
#focus on top 7, which all have more than 2% representation in the dataset
GBmap_top7 <- subset(GBmap, subset = scn_type_category_thresh0.15 %in% c("IPC-Glia", "Vascular", "Oligodendrocyte-Immature", "Astrocyte-Fibrous", "OPC", "IN-dLGE-Immature", "Astrocyte-Immature"))
GBmap_top7$scn_type_category_thresh0.15 <- factor(GBmap_top7$scn_type_category_thresh0.15, levels = c("IPC-Glia", "Astrocyte-Immature", "Astrocyte-Fibrous", "Vascular", "OPC", "Oligodendrocyte-Immature", "IN-dLGE-Immature"))
seven_colors <-  c("#E4E1E3FF", "#5A5156FF", "#F6222EFF", "#FEAF16FF", "#FE00FAFF", "#3283FEFF", "#16FF32FF")
names(seven_colors) <- c("IPC-Glia", "Vascular", "Oligodendrocyte-Immature", "Astrocyte-Fibrous", "OPC", "IN-dLGE-Immature", "Astrocyte-Immature")
top7_plot <- DimPlot_scCustom(GBmap_top7, group.by = "scn_type_category_thresh0.15", reduction = "umap", colors_use = seven_colors, pt.size = 1, raster = F) + labs(title = NULL)
state_plot + type_plot + top7_plot
ggsave(file="results/SingleCellNet/combined_umap.png",width = 24, height = 4, units = "in", dpi = 600)


#plot overall cell type proportions
pt_type <- table(GBmap$scn_type_category_thresh0.15)
pt_type <- as.data.frame(pt_type)
pt_type <- pt_type[order(pt_type$Freq),]
colnames(pt_type) <- c("type", "Freq")
pt_type$type <- factor(pt_type$type, levels = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                                "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                                                "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                                "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                                "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown", "Unmapped"))
pt_type$GBM <- "GBM"
ggplot(pt_type, aes(x = GBM, y = Freq, fill = type)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Proportion") +
  xlab("Overall") +
  scale_fill_manual(values = type_colors)
ggsave("results/SingleCellNet/pt_type.pdf", width = 7, height = 5, units = "in")

#plot cell type proportion per state
pt_state_type <- table(GBmap$annotation_level_3, GBmap$scn_type_category_thresh0.15)
pt_state_type <- as.data.frame(pt_state_type)
colnames(pt_state_type) <- c("state", "type", "Freq")
pt_state_type$type <- factor(pt_state_type$type, levels = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                                            "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                                                            "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                                            "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                                            "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown", "Unmapped"))
pt_scn_type_state_plot <- ggplot(pt_state_type, aes(x = state, y = Freq, fill = type)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Proportion") +
  xlab("State") +
  scale_fill_manual(values = type_colors)
ggsave("results/SingleCellNet/pt_state_type.pdf", width = 9, height = 5, units = "in")

#combine overall proportion with state specific proportion for a single plot
pt_type$state <- "All"
colnames(pt_type) <- c("type", "Freq", "state")
pt_state_type2 <- rbind(pt_type, pt_state_type)
pt_state_type2$type <- factor(pt_state_type2$type, levels = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                                              "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                                                              "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                                              "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                                              "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown", "Unmapped"))
pt_state_type2$state <- factor(pt_state_type2$state, levels = c("All", "AC-like", "MES-like", "NPC-like", "OPC-like"))
pt_state_type2_plot <- ggplot(pt_state_type2, aes(x = state, y = Freq, fill = type)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Proportion") +
  xlab("State") +
  scale_fill_manual(values = type_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("results/SingleCellNet/pt_state_type2.pdf", width = 9, height = 6, units = "in")

#plot cell type proportion per sample
pt_sample_type <- table(GBmap$donor_id, GBmap$scn_type_category_thresh0.15)
#only analyze samples with more than 50 malignant cells
pt_sample_type_filt <- pt_sample_type[rowSums(pt_sample_type) >= 50,]
pt_sample_type_filt <- as.data.frame(pt_sample_type_filt)
colnames(pt_sample_type_filt) <- c("sample", "type", "Freq")
pt_sample_type_filt$type <- factor(pt_sample_type_filt$type, levels = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                                                        "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                                                                        "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                                                        "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                                                        "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown", "Unmapped"))
pt_scn_type_sample_plot <- ggplot(pt_sample_type_filt, aes(x = sample, y = Freq, fill = type)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Proportion") +
  xlab("Sample") +
  scale_fill_manual(values = type_colors) +
  theme(axis.text.x = element_blank())
ggsave("results/SingleCellNet/pt_sample_type.pdf", width = 30, height = 5, units = "in")
