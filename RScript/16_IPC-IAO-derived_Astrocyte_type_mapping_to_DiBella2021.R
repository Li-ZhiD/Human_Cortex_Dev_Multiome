library(Seurat)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(orthogene)
library(singleCellNet)
library(future)

plan("multicore", workers = 12)
options(future.globals.maxSize = 200000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/20221216_GPC_differentiation/")

#import DiBella2021 data (ref) and GPC data (query)
GPC <- readRDS("results/final_GPC.rds")
load("DiBella2021/DiBella2021_scRNAseq_mouse_development.RData")

#########################
#DiBella data 2021 pre-processing
#########################
DiBella_AST_lineage <- subset(experiment.aggregate, subset = New_cellType %in% c("Astrocytes", "Cycling glial cells"))
#identify gene list that are one on one homologs to mouse in both GPC and DiBella dataset
#covert mouse gene in DiBella data to human orthologs
counts <- GetAssayData(GPC, assay = "RNA", slot = "counts")
counts_new_feature_name <- orthogene::convert_orthologs(gene_df = counts,
                                                        gene_input = "rownames", 
                                                        gene_output = "rownames", 
                                                        input_species = "human",
                                                        output_species = "mouse",
                                                        non121_strategy = "drop_both_species",
                                                        method = "gprofiler")
common_features <- intersect(rownames(counts_new_feature_name), rownames(DiBella_AST_lineage))

#only retain genes in the DiBella dataset that are common features
DiBella_AST_lineage <- DiBella_AST_lineage[common_features,]

# Run the standard workflow for visualization and clustering
DefaultAssay(DiBella_AST_lineage) <- "RNA"
DiBella_AST_lineage <- NormalizeData(DiBella_AST_lineage, verbose = FALSE)
DiBella_AST_lineage <- FindVariableFeatures(DiBella_AST_lineage)
#regress out cell cycle effect
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes.mouse <- orthogene::convert_orthologs(s.genes,
                                              gene_output = "columns", 
                                              input_species = "human",
                                              output_species = "mouse",
                                              non121_strategy = "keep_both_species",
                                              method = "gprofiler")
g2m.genes.mouse <- orthogene::convert_orthologs(g2m.genes,
                                                gene_output = "columns", 
                                                input_species = "human",
                                                output_species = "mouse",
                                                non121_strategy = "keep_both_species",
                                                method = "gprofiler")
DiBella_AST_lineage <- CellCycleScoring(DiBella_AST_lineage, s.features = s.genes.mouse$ortholog_gene, g2m.features = g2m.genes.mouse$ortholog_gene)
DiBella_AST_lineage <- ScaleData(DiBella_AST_lineage, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(DiBella_AST_lineage))
DiBella_AST_lineage <- RunPCA(DiBella_AST_lineage, verbose = FALSE)
ElbowPlot(DiBella_AST_lineage, ndims = 50)
ggsave(filename = "results/DiBella_integration/DiBella_AST_lineage_elbow.pdf", width = 10, height = 8)
DiBella_AST_lineage <- RunUMAP(DiBella_AST_lineage, reduction = "pca", dims = 1:20, return.model = TRUE)
DiBella_AST_lineage <- FindNeighbors(DiBella_AST_lineage, dims = 1:20)
DiBella_AST_lineage <- FindClusters(DiBella_AST_lineage, resolution = 0.8)


#check markers for astrocyte subtypes
DefaultAssay(DiBella_AST_lineage) <-  "RNA"
AST_markers <- c('Slc7a10', 'Serpinf1', 'H2az1', 'Sfrp1', 'S100a4', 'S100a11', 'Phlda1', 'Btbd17','Grm5', 'Slc6a11', 'Gabbr2', 'Igfbp2', 'Junb', 'Nr4a2', 'Tpt1', 'Rack1', 'Ybx1', 'Egr2', 'Egr3', 'Enkur',
                 'Prex2', 'Riiad1', 'Fbln2', 'Sparc', 'Vim', 'Cebpd', 'Igfbp5', 'S100a6','Eomes', 'Tnc', 'Dcx', 'Gad1', 'Gad2', 'Pdgfra', 'Plp1', 'Olig2', 'Bcas1', 'Mbp',
                 'Egfr', 'C1qa', 'Gfap', 'S100b', 'Aqp4', 'Sparcl1', 'Id1', 'Cldn5', 'Pecam1', 'Pdgfrb', 'Mki67', 'Ascl1', 'Dlx1', 'Dlx2', 'Dlx6', 'Sox2', 'Pax6', 'Hes5',
                 'Vim', "Nes", 'Hopx', 'Cryab', 'Pdgfd', 'Aldh1l1', 'Apoe', 'Cd38', 'Cd74', 'Ttr', 'Foxj1', 'Olig1', 'Gja1', 'Aldoc', 'Slc1a2', 'F3', 'Id2', 'Id3', 'Id4',
                 'Tshr', 'Hgf', 'CD38', 'Adgrv1', 'Rmst', 'Cd38','Tshr', 'Hgf', 'Slc16a12', 'Cep128', 'Lama1', 'Moxd1', 'Cd44', 'Glul', 'Slc1a2', 'Chrdl1', 'Il33', 'Eogt', 'Hes1', 'TFAP2C')
AST_markers_to_plot <- AST_markers[AST_markers %in% rownames(DiBella_AST_lineage)]
for (i in AST_markers_to_plot) {
  FeaturePlot_scCustom(DiBella_AST_lineage, reduction = "umap", features = i)
  ggsave(file=paste0("results/DiBella_integration/markers/", i,".png"),width = 6, height = 5, units = "in")
}

#annotate the lineages of individual clusters based on Zhou et al., 2023
annotation <- read.csv("results/DiBella_integration/annotation.csv")
annotation$RNA_snn_res.0.8 <-  as.factor(annotation$RNA_snn_res.0.8)
metadata <- DiBella_AST_lineage[[]]
metadata_annotated <- left_join(metadata, annotation, by = "RNA_snn_res.0.8")
rownames(metadata_annotated) <- rownames(metadata)
DiBella_AST_lineage <- AddMetaData(DiBella_AST_lineage, metadata_annotated)

#save data
saveRDS(DiBella_AST_lineage, "results/DiBella_integration/DiBella_AST_lineage.rds")

#filter data to focus on the two astrocyte lineages
DiBella_AST_lineage_filt <- subset(DiBella_AST_lineage, subset = lineage %in% c("Olig2", "S100a11"))
DiBella_AST_lineage_filt <- ScaleData(DiBella_AST_lineage_filt, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(DiBella_AST_lineage_filt))
DiBella_AST_lineage_filt <- RunPCA(DiBella_AST_lineage_filt, verbose = FALSE)
DiBella_AST_lineage_filt <- RunUMAP(DiBella_AST_lineage_filt, reduction = "pca", dims = 1:20, return.model = TRUE)

DimPlot_scCustom(DiBella_AST_lineage_filt, reduction = "umap", group.by = "lineage", raster = FALSE)
ggsave(filename = "results/DiBella_integration/DiBella_AST_lineage_filt_lineage.png", width = 6, height = 4)

#plot key markers: GFAP, IGFBP5, SPARC for S100a11 lineage and OLIG1, CHRDL1, BTBD17
DiBella_AST_lineage_filt.markers.plot.list <- lapply(X = c("Olig1", "Chrdl1", "Btbd17", "Gfap", "Igfbp5", "Sparc"), FUN = function(x) {
  FeaturePlot_scCustom(DiBella_AST_lineage_filt, reduction = "umap", features = x, raster = F)
})
patchwork::wrap_plots(DiBella_AST_lineage_filt.markers.plot.list, nrow = 2)
ggsave(filename = "results/DiBella_integration/DiBella_AST_lineage_filt.markers.plot.list.png", width = 15, height = 9)

#save data
saveRDS(DiBella_AST_lineage_filt, "results/DiBella_integration/DiBella_AST_lineage_filt.rds")

########################
#Seurat label transfer
########################
#focus on the whole astrocyte lineage RG-oRG, RG-tRG, IPC-IAO, and astrocytes derived from IPC-IAO
GPC_AST_lineage <- subset(GPC, subset = (Initial_cell_type == "tRG" & type == "RG") | (Initial_cell_type == "oRG" & type == "RG") | type == "IPC-Glia" | (Sample_ID %in% c("GW22_IPC-Glia_DIV7", "GW22_IPC-Glia_DIV14") & type == "Astrocyte"))

metadata <- GPC_AST_lineage[[]]
metadata <- metadata %>%  mutate(
  subtype = case_when(
    Initial_cell_type == "tRG" & type == "RG" ~ "tRG",
    Initial_cell_type == "oRG" & type == "RG" ~ "oRG",
    type == "IPC-Glia" ~ "IPC-IAO",
    TRUE ~ "Astrocyte"
  )
)
GPC_AST_lineage$subtype <- metadata$subtype

#remove excessive levels
refine_metadata_levels <- function(seurat_data){
  for (i in base::colnames(seurat_data@meta.data)){
    if (base::is.factor(seurat_data@meta.data[[i]])){
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
    }
  }
  return (seurat_data)
}
GPC_AST_lineage <- refine_metadata_levels(GPC_AST_lineage)

#covert human gene in IPC_IAO_AST to mouse orthologs
counts <- GetAssayData(GPC_AST_lineage, assay = "RNA", slot = "counts")
counts_new_feature_name <- orthogene::convert_orthologs(gene_df = counts,
                                                        gene_input = "rownames", 
                                                        gene_output = "rownames", 
                                                        input_species = "human",
                                                        output_species = "mouse",
                                                        non121_strategy = "drop_both_species",
                                                        method = "gprofiler")
#only retain genes in DiBella dataset that are common features
GPC_AST_lineage[["mouse"]] <- CreateAssayObject(counts = counts_new_feature_name[rownames(counts_new_feature_name) %in% common_features,])
DefaultAssay(GPC_AST_lineage)<- "mouse"
GPC_AST_lineage <- NormalizeData(GPC_AST_lineage)
GPC_AST_lineage <- FindVariableFeatures(GPC_AST_lineage)

#label transfer
anchors <- FindTransferAnchors(reference = DiBella_AST_lineage_filt, reference.reduction = "pca", query = GPC_AST_lineage, query.assay = "mouse", dims = 1:20)
GPC_AST_lineage <- MapQuery(anchorset = anchors, reference = DiBella_AST_lineage_filt, query = GPC_AST_lineage,
                            refdata = list(lineage = "lineage"), reference.reduction = "pca", reduction.model = "umap")
GPC_AST_lineage$predicted.lineage.score.filtered <- ifelse(GPC_AST_lineage$predicted.lineage.score > 0.5, GPC_AST_lineage$predicted.lineage, "Unmapped")

#save results
saveRDS(GPC_AST_lineage, "results/DiBella_integration/GPC_AST_lineage.rds")

########################################
#focus on astrocytes derived from IPC-IAO
GPC_AST <- subset(GPC, subset = Sample_ID %in% c("GW22_IPC-Glia_DIV7", "GW22_IPC-Glia_DIV14") & type == "Astrocyte")

#remove excessive levels
GPC_AST <- refine_metadata_levels(GPC_AST)

#covert human gene in IPC_IAO_AST to mouse orthologs
counts <- GetAssayData(GPC_AST, assay = "RNA", slot = "counts")
counts_new_feature_name <- orthogene::convert_orthologs(gene_df = counts,
                                                        gene_input = "rownames", 
                                                        gene_output = "rownames", 
                                                        input_species = "human",
                                                        output_species = "mouse",
                                                        non121_strategy = "drop_both_species",
                                                        method = "gprofiler")

#only retain genes in DiBella dataset that are common features
GPC_AST[["mouse"]] <- CreateAssayObject(counts = counts_new_feature_name[rownames(counts_new_feature_name) %in% common_features,])
DefaultAssay(GPC_AST)<- "mouse"
GPC_AST <- NormalizeData(GPC_AST)
GPC_AST <- FindVariableFeatures(GPC_AST)

#label transfer
anchors <- FindTransferAnchors(reference = DiBella_AST_lineage_filt, reference.reduction = "pca", query = GPC_AST, query.assay = "mouse", dims = 1:20)
GPC_AST <- MapQuery(anchorset = anchors, reference = DiBella_AST_lineage_filt, query = GPC_AST,
                            refdata = list(lineage = "lineage"), reference.reduction = "pca", reduction.model = "umap")
GPC_AST$predicted.lineage.score.filtered <- ifelse(GPC_AST$predicted.lineage.score > 0.5, GPC_AST$predicted.lineage, "Unmapped")

GPC_AST_lineage_plot <- DimPlot_scCustom(GPC_AST, group.by = "predicted.lineage.score.filtered", reduction = "ref.umap", raster = F, colors_use = c("#e63946", "#1d3557")) + xlim (c(-8, 8)) + ylim (c(-7, 6)) + ggtitle("")
DiBella_lineage_plot <- DimPlot_scCustom(DiBella_AST_lineage_filt, reduction = "umap", group.by = "lineage", raster = F, colors_use = c("#e63946", "#1d3557")) + xlim (c(-8, 8)) + ylim (c(-7, 6)) + ggtitle("")

DiBella_lineage_plot + GPC_AST_lineage_plot
ggsave("results/DiBella_integration/GPC_AST_mapped_lineage_ref_umap.png", width = 14, height = 6, units = "in")

#save results
saveRDS(GPC_AST, "results/DiBella_integration/GPC_AST.rds")

########################
#SingleCellNet classification
########################
#import query data
GPC_AST_lineage <- readRDS("results/DiBella_integration/GPC_AST_lineage.rds")
GPC_AST_lineage_SCN <- extractSeurat(GPC_AST_lineage, exp_slot_name = "counts")
GPC_AST_lineage_tab <-  GPC_AST_lineage_SCN$sampTab
GPC_AST_lineage_tab <- droplevels(GPC_AST_lineage_tab)
GPC_AST_lineage_tab$cell <- rownames(GPC_AST_lineage_tab)
GPC_AST_lineage_tab_exp <-  GPC_AST_lineage_SCN$expDat

#import training data
DiBella_AST_lineage_filt <- readRDS("results/DiBella_integration/DiBella_AST_lineage_filt.rds")
DefaultAssay(DiBella_AST_lineage_filt) <- "RNA"
DiBella_AST_lineage_filt_SCN <- extractSeurat(DiBella_AST_lineage_filt, exp_slot_name = "counts")
DiBella_AST_lineage_filt_tab <-  DiBella_AST_lineage_filt_SCN$sampTab
DiBella_AST_lineage_filt_tab <- droplevels(DiBella_AST_lineage_filt_tab)
DiBella_AST_lineage_filt_tab$cell <- rownames(DiBella_AST_lineage_filt_tab)
DiBella_AST_lineage_filt_exp <-  DiBella_AST_lineage_filt_SCN$expDat

#Find genes in common to the data sets and remove genes with variance = 1
RowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
goodGenes_GPC_AST_lineage <- rownames(GPC_AST_lineage_tab_exp)[RowVars(GPC_AST_lineage_tab_exp) != 0]
goodGenes_DiBella_AST_lineage_filt <- rownames(DiBella_AST_lineage_filt_exp)[RowVars(DiBella_AST_lineage_filt_exp) != 0]
commonGenes <-  intersect(goodGenes_GPC_AST_lineage, goodGenes_DiBella_AST_lineage_filt)

GPC_AST_lineage_tab_exp <-  GPC_AST_lineage_tab_exp[commonGenes,]
DiBella_AST_lineage_filt_exp <- DiBella_AST_lineage_filt_exp[commonGenes,]

##########################
#train a classifier for lineage
#Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=DiBella_AST_lineage_filt_tab, ncells = 800, dLevel="lineage")
stTrain = stList[[1]]
expTrain = DiBella_AST_lineage_filt_exp[,rownames(stTrain)]

#Train the classifier
#If you increase nTopGenes and nTopGenePairs, you may get a even better classifier performance on query data!
system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 200, nRand = 400, nTrees = 1000, nTopGenePairs = 200, dLevel = "lineage", colName_samp = "cell"))

#validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=800, dLevel="lineage") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = DiBella_AST_lineage_filt_exp[commonGenes,rownames(stTest)]

#predict testing data
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 100)
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "lineage", classQuery = "lineage", nRand = 100)
pdf("results/DiBella_integration/singleCellNet/tm_heldoutassessment_PR_type.pdf")
plot_PRs(tm_heldoutassessment)
dev.off()
pdf("results/DiBella_integration/singleCellNet/tm_heldoutassessment_metrics_type.pdf")
plot_metrics(tm_heldoutassessment)
dev.off()

tm_assess_matrix = tm_heldoutassessment$nonNA_PR
score = 0.35
celltype = "Olig2"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#"SCN score of 0.35 for cell type Olig2 has precision of 0.967 ~ 0.967 and sensitivity of 0.966 ~ 0.966"
score = 0.35
celltype = "S100a11"
calibration = scn_calibration(score = score, celltype = celltype, matrix=tm_assess_matrix)
#"SCN score of 0.35 for cell type S100a11 has precision of 0.961 ~ 0.961 and sensitivity of 0.985 ~ 0.984"


#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names
nrand = 100
sla = as.vector(stTest$lineage)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created
#heatmap
pdf("results/DiBella_integration/singleCellNet/evaluate_hmClass_type.pdf")
sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)
dev.off()
#bar graph
pdf("results/DiBella_integration/singleCellNet/evaluate_plot_attr_type.pdf")
plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="lineage", sid="cell")
dev.off()

#apply to query data
nqRand = 100
system.time(classRes_GPC_all<-scn_predict(class_info[['cnProc']], GPC_AST_lineage_tab_exp, nrand=nqRand))

# heatmap classification result
sgrp = as.vector(GPC_AST_lineage_tab$subtype)
names(sgrp) = as.vector(GPC_AST_lineage_tab$cell)
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)
pdf("results/DiBella_integration/singleCellNet/query_hmClass_type.pdf")
sc_hmClass(classRes_GPC_all, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)
dev.off()

pdf("results/DiBella_integration/singleCellNet/query_plot_attr_type.pdf")
plot_attr(classRes_GPC_all, GPC_AST_lineage_tab, nrand=nqRand, sid="cell", dLevel="subtype")
dev.off()

# This classifies a cell with  the category with the highest classification score or higher than a classification score threshold of your choosing (we choose 0.35 here).
GPC_AST_lineage_tab <- get_cate(classRes = classRes_GPC_all, sampTab = GPC_AST_lineage_tab, dLevel = "subtype", sid = "cell", nrand = nqRand, cThresh = 0.35)

#save results
save(class_info, classRes_val_all, classRes_GPC_all, file = "results/DiBella_integration/singleCellNet/scn_models_type.RData")




