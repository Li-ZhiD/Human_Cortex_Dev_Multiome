#remotes::install_github("satijalab/seurat", "seurat5")
library(Seurat)
library(future)
library(tidyverse)
library(scCustomize)
library(ComplexHeatmap)
library(patchwork)
library(cowplot)
plan("multisession", workers = 4)
options(future.globals.maxSize = 2000000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/MERFISH/analysis")

#load data
ARKFrozen62PFC <- readRDS("ARKFrozen62PFC.rds")
NIH4365BA10 <- readRDS("NIH4365BA10.rds")
UCSF2018003MFG <- readRDS("UCSF2018003MFG.rds")
ARKFrozen65V1 <- readRDS("ARKFrozen65V1.rds")
NIH5900BA17 <- readRDS("NIH5900BA17.rds")
NIH4392BA17 <- readRDS("NIH4392BA17.rds")

#merge data
seurat.obj.list <- c(ARKFrozen62PFC, NIH4365BA10, UCSF2018003MFG, ARKFrozen65V1, NIH5900BA17, NIH4392BA17)
seurat.obj.list.names <- c("ARKFrozen62PFC", "NIH4365BA10", "UCSF2018003MFG", "ARKFrozen65V1", "NIH5900BA17", "NIH4392BA17")
for (i in 1:6) {
  VlnPlot(seurat.obj.list[[i]], features = c("nCount_Vizgen", "nFeature_Vizgen", "volume"), ncol = 3, pt.size = 0, log = T)
  ggsave(file=paste0("results/QC/",seurat.obj.list.names[i],".png"),width = 6, height = 6, units = "in")
}

#filter out cells that are too small (volume < 10 um3), with too few or too many counts (nCount_Vizgen < 25 or > 2000), or with too few genes (nFeature_Vizgen < 10)
seurat.obj.list.filt <- lapply(seurat.obj.list, FUN = subset, subset = nCount_Vizgen >=25 & nCount_Vizgen <= 2000 & nFeature_Vizgen >= 10 & volume >= 10)
names(seurat.obj.list.filt) <- c("ARKFrozen62PFC", "NIH4365BA10", "UCSF2018003MFG", "ARKFrozen65V1", "NIH5900BA17", "NIH4392BA17")
#save data
saveRDS(seurat.obj.list.filt, "filtered_datasets.rds")

#data normalization and integration
seurat.obj.list.filt <- lapply(X = seurat.obj.list.filt, FUN = SCTransform, assay = "Vizgen", clip.range = c(-10, 10))
features <- rownames(seurat.obj.list.filt[[1]])
seurat.obj.list.filt <- PrepSCTIntegration(object.list = seurat.obj.list.filt, anchor.features = features)
seurat.obj.list.filt <- lapply(X = seurat.obj.list.filt, FUN = RunPCA, features = features)
anchors <- FindIntegrationAnchors(object.list = seurat.obj.list.filt, normalization.method = "SCT",
                                  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
saveRDS(anchors, "anchors.rds")
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
saveRDS(combined.sct, "integrated_datasets.rds")

#add metadata
#rename sample ID
combined.sct$Sample_ID <- case_when(
  Cells(combined.sct) %in% Cells(seurat.obj.list.filt[[1]]) ~ "ARKFrozen-62-PFC",
  Cells(combined.sct) %in% Cells(seurat.obj.list.filt[[2]]) ~ "NIH-4365-BA10",
  Cells(combined.sct) %in% Cells(seurat.obj.list.filt[[3]]) ~ "UCSF2018-003-MFG",
  Cells(combined.sct) %in% Cells(seurat.obj.list.filt[[4]]) ~ "ARKFrozen-65-V1",
  Cells(combined.sct) %in% Cells(seurat.obj.list.filt[[5]]) ~ "NIH-5900-BA17",
  Cells(combined.sct) %in% Cells(seurat.obj.list.filt[[6]]) ~ "NIH-4392-BA17",
)
metadata <- read.csv("sample_information.csv")
old_metadata <- combined.sct[[]]
new_metadata <- left_join(old_metadata, metadata)
rownames(new_metadata) <- rownames(old_metadata)
head(new_metadata)
combined.sct <- AddMetaData(combined.sct, new_metadata)
combined.sct$Sample_ID <- factor(combined.sct$Sample_ID, levels = c("ARKFrozen-62-PFC", "NIH-4365-BA10", "UCSF2018-003-MFG", "ARKFrozen-65-V1", "NIH-5900-BA17", "NIH-4392-BA17"))
combined.sct$Group <- factor(combined.sct$Group, levels = c("Second_trimester", "Third_trimester", "Infancy"))
combined.sct$Region <- factor(combined.sct$Region, levels = c("PFC", "V1"))

#dimension reduction and clustering
DefaultAssay(combined.sct) <- "integrated"
combined.sct <- RunPCA(combined.sct)
ElbowPlot(combined.sct, ndims = 50, reduction = "pca")
ggsave("results/pre_filtering/pca_elbpw.png")
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = 2.8)
combined.sct$integrated_snn_res.2.8 <- factor(combined.sct$integrated_snn_res.2.8, levels = 0:97)

DimPlot_scCustom(combined.sct, group.by = "integrated_snn_res.2.8", reduction = "umap", label = TRUE, label.size = 6) + theme(legend.text = element_text(size = 18))
ggsave(file = "results/pre_filtering/umap_cluster_integrated_snn_res.2.8.png", width = 12, height = 8, units = "in")

#correct SCT counts to be comparable between datasets
DefaultAssay(combined.sct) <- "SCT"
combined.sct <- PrepSCTFindMarkers(combined.sct)

saveRDS(combined.sct, "integrated_clustered_sct_corrected.rds")

#add annotation
metadata <- read.csv("results/MERFISH_annotation.csv")
metadata$integrated_snn_res.2.8 <- factor(metadata$integrated_snn_res.2.8, levels = 0:97)
metadata$type <- factor(metadata$type, levels = c("RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                                    "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-L5-ET", "EN-L5_6-NP",
                                                    "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                                    "IN-MGE-PV", "IPC-Glia", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                                    "Oligodendrocyte-Immature", "Oligodendrocyte", "Microglia", "Vascular", "Low quality"))
old_metadata <- combined.sct[[]]
new_metadata <- left_join(old_metadata, metadata)
rownames(new_metadata) <- rownames(old_metadata)
head(new_metadata)
combined.sct <- AddMetaData(combined.sct, new_metadata)

#remove low quality cells
combined.filt <- subset(combined.sct, subset = type != 'Low quality')

#save annotated data
saveRDS(combined.filt, "annotated_filtered.rds")

#spatial plots for individual cell types in each sample
#ARKFrozen62PFC
ARKFrozen62PFC_type <- ImageDimPlot(combined.filt, fov = "ARKFrozen62PFC", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                       "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                       "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                       "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                       "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                       "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                       "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = F, border.size = 0, dark.background = F)

ARKFrozen62PFC_type_split <- ImageDimPlot(combined.filt, fov = "ARKFrozen62PFC", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                             "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                             "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                             "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                             "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                             "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                             "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, split.by = "type", dark.background = F) + NoLegend()
#delete other fovs, because otherwise we cannot subset the object
tmp <- combined.filt
tmp[["NIH4365BA10"]] <-  NULL
tmp[["UCSF2018003MFG"]] <-  NULL
tmp[["ARKFrozen65V1"]] <-  NULL
tmp[["NIH5900BA17"]] <-  NULL
tmp[["NIH4392BA17"]] <-  NULL
ARKFrozen62PFC <- subset(tmp, subset = Sample_ID == "ARKFrozen-62-PFC")
ARKFrozen62PFC.niche <- BuildNicheAssay(object = ARKFrozen62PFC, fov = "ARKFrozen62PFC", group.by = "type",
                                     niches.k = 8, neighbors.k = 50)


#NIH4365BA10
NIH4365BA10_type <- ImageDimPlot(combined.filt, fov = "NIH4365BA10", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                 "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                 "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                 "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                 "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                 "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                 "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = F, border.size = 0, dark.background = F)
NIH4365BA10_type_split <- ImageDimPlot(combined.filt, fov = "NIH4365BA10", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                       "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                       "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                       "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                       "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                       "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                       "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, split.by = "type", dark.background = F) + NoLegend()
#delete other fovs, because otherwise we cannot subset the object
tmp <- combined.filt
tmp[["ARKFrozen62PFC"]] <-  NULL
tmp[["UCSF2018003MFG"]] <-  NULL
tmp[["ARKFrozen65V1"]] <-  NULL
tmp[["NIH5900BA17"]] <-  NULL
tmp[["NIH4392BA17"]] <-  NULL
NIH4365BA10 <- subset(tmp, subset = Sample_ID == "NIH-4365-BA10")
NIH4365BA10.niche <- BuildNicheAssay(object = NIH4365BA10, fov = "NIH4365BA10", group.by = "type",
                                        niches.k = 6, neighbors.k = 50)


#UCSF2018003MFG
UCSF2018003MFG_type <- ImageDimPlot(combined.filt, fov = "UCSF2018003MFG", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                       "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                       "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                       "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                       "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                       "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                       "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = F, border.size = 0, dark.background = F)
UCSF2018003MFG_type_split <- ImageDimPlot(combined.filt, fov = "UCSF2018003MFG", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                             "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                             "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                             "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                             "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                             "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                             "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, split.by = "type", dark.background = F) + NoLegend()
#delete other fovs, because otherwise we cannot subset the object
tmp <- combined.filt
tmp[["ARKFrozen62PFC"]] <-  NULL
tmp[["NIH4365BA10"]] <-  NULL
tmp[["ARKFrozen65V1"]] <-  NULL
tmp[["NIH5900BA17"]] <-  NULL
tmp[["NIH4392BA17"]] <-  NULL
UCSF2018003MFG <- subset(tmp, subset = Sample_ID == "UCSF2018-003-MFG")
UCSF2018003MFG.niche <- BuildNicheAssay(object = UCSF2018003MFG, fov = "UCSF2018003MFG", group.by = "type",
                                     niches.k = 8, neighbors.k = 50)


#ARKFrozen65V1
ARKFrozen65V1_type <- ImageDimPlot(combined.filt, fov = "ARKFrozen65V1", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                     "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                     "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                     "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                     "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                     "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                     "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = F, border.size = 0, dark.background = F)
ARKFrozen65V1_type_split <- ImageDimPlot(combined.filt, fov = "ARKFrozen65V1", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                           "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                           "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                           "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                           "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                           "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                           "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, split.by = "type", dark.background = F) + NoLegend()
#delete other fovs, because otherwise we cannot subset the object
tmp <- combined.filt
tmp[["ARKFrozen62PFC"]] <-  NULL
tmp[["NIH4365BA10"]] <-  NULL
tmp[["UCSF2018003MFG"]] <-  NULL
tmp[["NIH5900BA17"]] <-  NULL
tmp[["NIH4392BA17"]] <-  NULL
ARKFrozen65V1 <- subset(tmp, subset = Sample_ID == "ARKFrozen-65-V1")
ARKFrozen65V1.niche <- BuildNicheAssay(object = ARKFrozen65V1, fov = "ARKFrozen65V1", group.by = "type",
                                        niches.k = 7, neighbors.k = 50)

#NIH5900BA17
NIH5900BA17_type <- ImageDimPlot(combined.filt, fov = "NIH5900BA17", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                 "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                 "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                 "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                 "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                 "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                 "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = F, border.size = 0, dark.background = F)
NIH5900BA17_type_split <- ImageDimPlot(combined.filt, fov = "NIH5900BA17", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                       "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                       "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                       "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                       "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                       "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                       "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, split.by = "type", dark.background = F) + NoLegend()
#delete other fovs, because otherwise we cannot subset the object
tmp <- combined.filt
tmp[["ARKFrozen62PFC"]] <-  NULL
tmp[["NIH4365BA10"]] <-  NULL
tmp[["UCSF2018003MFG"]] <-  NULL
tmp[["ARKFrozen65V1"]] <-  NULL
tmp[["NIH4392BA17"]] <-  NULL
NIH5900BA17 <- subset(tmp, subset = Sample_ID == "NIH-5900-BA17")
NIH5900BA17.niche <- BuildNicheAssay(object = NIH5900BA17, fov = "NIH5900BA17", group.by = "type",
                                       niches.k = 7, neighbors.k = 50)


#NIH4392BA17
NIH4392BA17_type <- ImageDimPlot(combined.filt, fov = "NIH4392BA17", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                 "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                 "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                 "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                 "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                 "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = F, border.size = 0, dark.background = F)
NIH4392BA17_type_split <- ImageDimPlot(combined.filt, fov = "NIH4392BA17", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                       "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                       "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                       "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                       "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                       "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, split.by = "type", dark.background = F) + NoLegend()
#delete other fovs, because otherwise we cannot subset the object
tmp <- combined.filt
tmp[["ARKFrozen62PFC"]] <-  NULL
tmp[["NIH4365BA10"]] <-  NULL
tmp[["UCSF2018003MFG"]] <-  NULL
tmp[["ARKFrozen65V1"]] <-  NULL
tmp[["NIH5900BA17"]] <-  NULL
NIH4392BA17 <- subset(tmp, subset = Sample_ID == "NIH-4392-BA17")
NIH4392BA17.niche <- BuildNicheAssay(object = NIH4392BA17, fov = "NIH4392BA17", group.by = "type",
                                     niches.k = 6, neighbors.k = 50)


#combine all niche assay together for whole dataset kmeans analysis
ARKFrozen62PFC.niche.assay <- GetAssay(ARKFrozen62PFC.niche, assay = "niche")
NIH4365BA10.niche.assay <- GetAssay(NIH4365BA10.niche, assay = "niche")
UCSF2018003MFG.niche.assay <- GetAssay(UCSF2018003MFG.niche, assay = "niche")
ARKFrozen65V1.niche.assay <- GetAssay(ARKFrozen65V1.niche, assay = "niche")
NIH5900BA17.niche.assay <- GetAssay(NIH5900BA17.niche, assay = "niche")
NIH4392BA17.niche.assay <- GetAssay(NIH4392BA17.niche, assay = "niche")

combined.niche.assay <- merge(ARKFrozen62PFC.niche.assay, y = c(NIH4365BA10.niche.assay, UCSF2018003MFG.niche.assay, ARKFrozen65V1.niche.assay, NIH5900BA17.niche.assay, NIH4392BA17.niche.assay))
combined.filt[["niche"]] <- combined.niche.assay
DefaultAssay(combined.filt) <- "niche"
combined.filt <- ScaleData(combined.filt)
results <- kmeans(x = t(combined.filt[["niche"]]@scale.data), centers = 10, 
                  nstart = 50, iter.max = 30)
combined.filt$niches <- results[["cluster"]]
combined.filt$niches <- factor(combined.filt$niches, levels = 1:10)
combined.filt$type <- droplevels(combined.filt$type)


######################
#import manual niche annotation
metadata <- read.csv("results/niche_annotation.csv")
metadata$niches <- factor(metadata$niches, levels = 1:10)
old_metadata <- combined.filt[[]]
new_metadata <- left_join(old_metadata, metadata)
rownames(new_metadata) <- rownames(old_metadata)
head(new_metadata)
combined.filt <- AddMetaData(combined.filt, new_metadata)
combined.filt$niche_name <- factor(combined.filt$niche_name, levels = c("Cortical layer 1", "Cortical upper layer in development", "Cortical layer 2 and 3", "Cortical layer 4", "Cortical layer 5",
                                                                        "Cortical layer 6 and subplate", "Intermediate zone", "White matter and Meninge", "Ventricular zone and subventricular zone", "Dorsal lateral ganglionic eminence"))

#save annotated data
saveRDS(combined.filt, "annotated_filtered_niche_identified.rds")
as.data.frame.matrix(table(combined.filt$type, combined.filt$niche_name))

#save cell-level metadata
meta <- combined.filt[[]]
meta$Cell_ID <- rownames(meta)
meta_filt <- meta[,c("Cell_ID", "Sample_ID", "Estimated_postconceptional_age_in_days", "Group", "Region", "nCount_Vizgen", "nFeature_Vizgen", "volume", "class", "subclass", "type", "seurat_clusters", "niche_name")]
colnames(meta_filt) <- c("Cell_ID", "Sample_ID", "Estimated_postconceptional_age_in_days", "Group", "Region", "nCount_Vizgen", "nFeature_Vizgen", "Volume", "Class", "Subclass", "Type", "Cluster", "Niche")
write.csv(meta_filt, "cell_level_metadata.csv", row.names = F)

##########################
#plots in manuscript
#heatmap of proportion of cell types in each niche
mat <- table(combined.filt$niche_name, combined.filt$type)
mat <- matrix(mat, ncol = ncol(mat), dimnames = dimnames(mat))
write.csv(mat, file = "results/niche_type_matrix.csv")
mat <- t(apply(mat,1, function(x) x/sum(x)))

col_fun = circlize::colorRamp2(breaks = c(0,0.4), hcl_palette = "Blue-Yellow")
row_ha <- rowAnnotation(Niche = rownames(mat),
                        col = list(Niche = c("Cortical layer 1" = "#4E79A7", "Cortical upper layer in development" = "#F28E2B", "Cortical layer 2 and 3" = "#E15759", "Cortical layer 4" = "#76B7B2", "Cortical layer 5" = "#59A14F", "Cortical layer 6 and subplate" = "#EDC948", "Intermediate zone" = "#B07AA1", "White matter and Meninge" = "#FF9DA7", "Ventricular zone and subventricular zone" = "#9C755F", "Dorsal lateral ganglionic eminence" = "#BAB0AC")),
                        gp = gpar(col = "white"),
                        simple_anno_size = unit(4, "mm"),
                        show_legend = F)
column_ha <-  HeatmapAnnotation(Type = colnames(mat),
                                col = list(Type = c("RG-tRG" = "#F6222EFF", "RG-oRG" = "#FE00FAFF", "IPC-EN" = "#16FF32FF", "EN-Newborn" = "#3283FEFF", "EN-IT-Immature" = "#FEAF16FF",
                                                    "EN-L2_3-IT" = "#B00068FF", "EN-L4-IT" = "#1CFFCEFF", "EN-L5-IT" = "#90AD1CFF", "EN-L6-IT" = "#2ED9FFFF", "EN-L5-ET" = "#AA0DFEFF",
                                                    "EN-L5_6-NP" = "#F8A19FFF", "EN-L6-CT" = "#325A9BFF", "EN-L6b" = "#C4451CFF", "IN-dLGE-Immature" = "#1C8356FF", "IN-CGE-Immature" = "#85660DFF",
                                                    "IN-CGE-VIP" = "#B10DA1FF", "IN-CGE-SNCG" = "#FBE426FF", "IN-CGE-LAMP5" =  "#1CBE4FFF", "IN-MGE-Immature" = "#FA0087FF", "IN-MGE-SST" = "#FC1CBFFF",
                                                    "IN-MGE-PV" = "#F7E1A0FF", "IPC-Glia" = "#C075A6FF", "Astrocyte-Protoplasmic" = "#AAF400FF", "Astrocyte-Fibrous" = "#BDCDFFFF", "OPC" = "#B5EFB5FF",
                                                    "Oligodendrocyte-Immature" = "#7ED7D1FF", "Oligodendrocyte" = "#1C7F93FF", "Microglia" = "#683B79FF", "Vascular" = "#3B00FBFF")),
                                gp = gpar(col = "white"),
                                simple_anno_size = unit(4, "mm"),
                                show_legend = F)

ht <- Heatmap(mat,
              name = "Proportion of the niche",
              cluster_rows = F,
              cluster_columns = F,
              col = col_fun,
              row_names_side = "left",
              column_names_side = "top",
              column_names_rot = 45,
              left_annotation = row_ha,
              top_annotation = column_ha,
              width = ncol(mat)*unit(6, "mm"), 
              height = nrow(mat)*unit(6, "mm"))
pdf(file = "results/niche_composition_heatmap.pdf", width = 14, height = 5)
draw(ht, heatmap_legend_side = "right")
dev.off()

#Immature interneuron proportions in 2nd trimester MZ vs VZ/SVZ
metadata <- combined.filt[[]]
metadata_filt <- metadata %>% filter(Group == "Second_trimester")
mat <- table(metadata_filt$niche_name, metadata_filt$type)
mat <- matrix(mat, ncol = ncol(mat), dimnames = dimnames(mat))
write.csv(mat, file = "results/niche_type_matrix_2nd_trimester.csv")

#EN-Newborn vs IN-dLGE-Immature in infancy VZ/SVZ 
metadata <- combined.filt[[]]
metadata_filt <- metadata %>% filter(Sample_ID == "UCSF2018-003-MFG")
mat <- table(metadata_filt$niche_name, metadata_filt$type)
mat <- matrix(mat, ncol = ncol(mat), dimnames = dimnames(mat))
write.csv(mat, file = "results/niche_type_matrix_UCSF2018-003-MFG.csv")

##########################
#dimplot for niches in individual samples
ARKFrozen62PFC.niche.plot.concensus <- ImageDimPlot(combined.filt, fov = "ARKFrozen62PFC",group.by = "niche_name", size = 0.75, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"))
NIH4365BA10.niche.plot.concensus <- ImageDimPlot(combined.filt, fov = "NIH4365BA10",group.by = "niche_name", size = 0.75, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"))
UCSF2018003MFG.niche.plot.concensus <- ImageDimPlot(combined.filt, fov = "UCSF2018003MFG",group.by = "niche_name", size = 0.75, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#FF9DA7", "#9C755F"))
ARKFrozen65V1.niche.plot.concensus <- ImageDimPlot(combined.filt, fov = "ARKFrozen65V1",group.by = "niche_name", size = 0.75, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"))
NIH5900BA17.niche.plot.concensus <- ImageDimPlot(combined.filt, fov = "NIH5900BA17",group.by = "niche_name", size = 0.75, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F"))
NIH4392BA17.niche.plot.concensus <- ImageDimPlot(combined.filt, fov = "NIH4392BA17",group.by = "niche_name", size = 0.75, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7"))

#remove legend for plotting
ARKFrozen62PFC_type_no_legend <- ARKFrozen62PFC_type + NoLegend()
ARKFrozen62PFC.niche.plot.concensus_no_legend <- ARKFrozen62PFC.niche.plot.concensus + NoLegend()
NIH4365BA10_type_no_legend <- NIH4365BA10_type + NoLegend()
NIH4365BA10.niche.plot.concensus_no_legend <- NIH4365BA10.niche.plot.concensus + NoLegend()
UCSF2018003MFG_type_no_legend <- UCSF2018003MFG_type + NoLegend()
UCSF2018003MFG.niche.plot.concensus_no_legend <- UCSF2018003MFG.niche.plot.concensus + NoLegend()
ARKFrozen65V1_type_no_legend <- ARKFrozen65V1_type + NoLegend()
ARKFrozen65V1.niche.plot.concensus_no_legend <- ARKFrozen65V1.niche.plot.concensus + NoLegend()
NIH5900BA17_type_no_legend <- NIH5900BA17_type + NoLegend()
NIH5900BA17.niche.plot.concensus_no_legend <- NIH5900BA17.niche.plot.concensus + NoLegend()
NIH4392BA17_type_no_legend <- NIH4392BA17_type + NoLegend()
NIH4392BA17.niche.plot.concensus_no_legend <- NIH4392BA17.niche.plot.concensus + NoLegend()

#plot all cell types and niches in one plot
ARKFrozen62PFC_type_no_legend + ARKFrozen62PFC.niche.plot.concensus_no_legend +
  NIH4365BA10_type_no_legend + NIH4365BA10.niche.plot.concensus_no_legend +
  UCSF2018003MFG_type_no_legend + UCSF2018003MFG.niche.plot.concensus_no_legend +
  ARKFrozen65V1_type_no_legend + ARKFrozen65V1.niche.plot.concensus_no_legend +
  NIH5900BA17_type_no_legend + NIH5900BA17.niche.plot.concensus_no_legend +
  NIH4392BA17_type_no_legend + NIH4392BA17.niche.plot.concensus_no_legend + patchwork::plot_layout(ncol = 6)
ggsave("results/final_type_niche_combined_plot.png", width = 30, height = 12)
#plot legends
ImageDimPlot(combined.filt, fov = "ARKFrozen62PFC", group.by = "type", cols = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                "#1C7F93FF", "#683B79FF", "#3B00FBFF"), size = 0.75, axes = T, border.size = 0, dark.background = F)
ggsave("results/final_type_legend_and_scale_bar.pdf", width = 8, height = 6)
ImageDimPlot(combined.filt, fov = "ARKFrozen62PFC",group.by = "niche_name", size = 0.75, axes = T, dark.background = F, cols = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"))
ggsave("results/final_niche_name_legend_and_scale_bar.pdf", width = 8, height = 6)

#plot all split plots in one plot for Extended Data Fig. 4
ARKFrozen62PFC_type_split <- ARKFrozen62PFC_type_split + xlim(0,13500) + ylim(0,13500) + theme(strip.text.x = element_text(size = 8), plot.margin = margin(1, 1, 1, 1, "cm"))
NIH4365BA10_type_split <- NIH4365BA10_type_split + xlim(0,13500) + ylim(0,13500)  + theme(strip.text.x = element_text(size = 8), plot.margin = margin(1, 1, 1, 1, "cm"),)
UCSF2018003MFG_type_split <- UCSF2018003MFG_type_split + xlim(0,13500) + ylim(0,13500) + theme(strip.text.x = element_text(size = 8), plot.margin = margin(1, 1, 1, 1, "cm"))
ARKFrozen65V1_type_split <- ARKFrozen65V1_type_split + xlim(0,13500) + ylim(0,13500) + theme(strip.text.x = element_text(size = 8), plot.margin = margin(1, 1, 1, 1, "cm"))
NIH5900BA17_type_split <- NIH5900BA17_type_split + xlim(0,13500) + ylim(0,13500) + theme(strip.text.x = element_text(size = 8), plot.margin = margin(1, 1, 1, 1, "cm"))
NIH4392BA17_type_split <- NIH4392BA17_type_split + xlim(0,13500) + ylim(0,13500) + theme(strip.text.x = element_text(size = 8), plot.margin = margin(1, 1, 1, 1, "cm"))

ARKFrozen62PFC_type_split + NIH4365BA10_type_split + UCSF2018003MFG_type_split + ARKFrozen65V1_type_split + NIH5900BA17_type_split + NIH4392BA17_type_split + patchwork::plot_layout(ncol = 2)
ggsave("results/final_type_split_plot.png", width = 20, height = 30)

######################
#plot QC results in Extended Data Fig. 4
combined.filt$Sample_ID <- factor(combined.filt$Sample_ID, levels = c("ARKFrozen-62-PFC", "ARKFrozen-65-V1", "NIH-4365-BA10", "NIH-5900-BA17", "UCSF2018-003-MFG", "NIH-4392-BA17"))
#UMAP figures in Extended Data Fig. 4
sample_id_plot <- DimPlot_scCustom(combined.filt, group.by = "Sample_ID", reduction = "umap", colors_use = c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"), pt.size = 0.1, raster = F) + NoLegend() + NoAxes() + labs(title = NULL)
DimPlot_scCustom(combined.filt, group.by = "Sample_ID", reduction = "umap", colors_use = c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C"), pt.size = 0.1, raster = T)
ggsave(file="results/sample_id_plot_legend.pdf", width = 10, height = 8, units = "in")

group_plot <- DimPlot_scCustom(combined.filt, group.by = "Group", reduction = "umap", colors_use = c("#f0f921", "#e16462", "#6a00a8"), pt.size = 0.1, raster = F) + NoLegend() + NoAxes() + labs(title = NULL)
DimPlot_scCustom(combined.filt, group.by = "Group", reduction = "umap", colors_use = c("#f0f921", "#e16462", "#6a00a8"), pt.size = 0.1, raster = T)
ggsave(file="results/group_plot_legend.pdf", width = 10, height = 8, units = "in")

region_plot <- DimPlot_scCustom(combined.filt, group.by = "Region", reduction = "umap", colors_use = c("#ffa500", "#0072b2"), pt.size = 0.1, raster = F) + NoLegend() + NoAxes() + labs(title = NULL)
DimPlot_scCustom(combined.filt, group.by = "Region", reduction = "umap", colors_use = c("#ffa500", "#0072b2"), pt.size = 0.1, raster = T)
ggsave(file="results/region_plot_legend.pdf", width = 10, height = 8, units = "in")

type_plot <- DimPlot_scCustom(combined.filt, group.by = "type", reduction = "umap", colors_use = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                   "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                   "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                   "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                   "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                   "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                   "#1C7F93FF", "#683B79FF", "#3B00FBFF"), pt.size = 0.1, raster = F) + NoLegend() + NoAxes() + labs(title = NULL)
DimPlot_scCustom(combined.filt, group.by = "type", reduction = "umap", colors_use = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                   "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                   "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                   "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                   "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                   "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                   "#1C7F93FF", "#683B79FF", "#3B00FBFF"), pt.size = 0.1, raster = T)
ggsave(file="results/type_plot_legend.pdf", width = 10, height = 8, units = "in")

niche_plot <- DimPlot_scCustom(combined.filt, group.by = "niche_name", reduction = "umap", colors_use = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"), pt.size = 0.1, raster = F) + NoLegend() + NoAxes() + labs(title = NULL)
DimPlot_scCustom(combined.filt, group.by = "niche_name", reduction = "umap", colors_use = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC"), pt.size = 0.1, raster = T)
ggsave(file="results/niche_plot_legend.pdf", width = 10, height = 8, units = "in")

sample_id_plot + type_plot + group_plot + niche_plot + region_plot + plot_layout(ncol = 2)
ggsave(file="results/extended_data_fig_umaps.png",width = 8.33, height = 11.25, units = "in",)

#########
#plot QC metrics violin plots in Extended Data Fig. 4
metadata <- combined.filt[[]]
nCount_RNA_plot <- ggplot(metadata, aes(x = Sample_ID, y = nCount_Vizgen)) + 
  geom_violin(aes(fill = Sample_ID), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab("# transcripts") +
  scale_fill_manual(values = c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")) +
  coord_cartesian(ylim = c(0, 1000)) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/nCount_RNA_per_dataset.png", width = 4, height = 2, units = "in")

nFeature_RNA_plot <- ggplot(metadata, aes(x = Sample_ID, y = nFeature_Vizgen)) + 
  geom_violin(aes(fill = Sample_ID), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab("# genes") +
  scale_fill_manual(values = c("#EF476F", "#F78C6B", "#FFD166", "#06D6A0", "#118AB2", "#073B4C")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/nFeature_RNA_per_dataset.png", width = 4, height = 2, units = "in")

#plot cell composition/proportion in individual samples
pt_sample <- table(metadata$type, metadata$Sample_ID)
pt_sample <- data.frame(pt_sample)
type_proportion_per_sample <-  metadata %>% group_by(Sample_ID) %>% mutate(total_cell_number = n()) %>% ungroup() %>% group_by(Sample_ID, total_cell_number, type) %>% summarise(type_number = n()) %>% mutate(type_proportion = type_number/total_cell_number)
write.csv(type_proportion_per_sample, "results/type_proportion_per_sample.csv", row.names = T)

type_porportion_plot <- ggplot(pt_sample, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.5) +
  xlab("Dataset") +
  ylab("Type proportion") +
  theme_classic() +
  scale_fill_manual(values = c("#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                               "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                               "#2ED9FFFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                               "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                               "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                               "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                               "#1C7F93FF", "#683B79FF", "#3B00FBFF")) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(ncol= 3))
ggsave("results/type_porportion_per_dataset.pdf", width = 10, height = 4, units = "in")
type_porportion_plot_no_legend <- type_porportion_plot + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank())

#plot niche proportion in individual samples
pt_niche <- table(metadata$niche_name, metadata$Sample_ID)
pt_niche <- data.frame(pt_niche)
niche_proportion_per_sample <-  metadata %>% group_by(Sample_ID) %>% mutate(total_cell_number = n()) %>% ungroup() %>% group_by(Sample_ID, total_cell_number, niche_name) %>% summarise(niche_number = n()) %>% mutate(niche_proportion = niche_number/total_cell_number)
write.csv(niche_proportion_per_sample, "results/niche_proportion_per_sample.csv", row.names = T)

niche_porportion_plot <- ggplot(pt_niche, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.5) +
  xlab("Dataset") +
  ylab("Niche proportion") +
  theme_classic() +
  scale_fill_manual(values = c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC")) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(ncol= 3))
ggsave("results/niche_porportion_per_dataset.pdf", width = 10, height = 4, units = "in")
niche_porportion_plot_no_legend <- niche_porportion_plot + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank())

#add annotation bar
annotation_bar_info <- distinct(metadata[,c("Sample_ID", "Region", "Group")])
age_group_bar <- ggplot(annotation_bar_info) +
  geom_bar(mapping = aes(x = Sample_ID, y = 1, fill = Group), 
           stat = "identity", 
           width = 0.85) +
  scale_fill_manual(values = c("#f0f921", "#e16462", "#6a00a8")) +
  theme_void() +
  theme(legend.direction="horizontal", panel.spacing.x = unit(1, "mm"), plot.margin = unit(c(0, 0, 2, 0), "mm"))

region_bar <- ggplot(annotation_bar_info) +
  geom_bar(mapping = aes(x = Sample_ID, y = 1, fill = Region), 
           stat = "identity", 
           width = 0.85) +
  scale_fill_manual(values = c("#ffa500", "#0072b2")) +
  theme_void() +
  theme(legend.direction="horizontal", panel.spacing.x = unit(1, "mm"), plot.margin = unit(c(0, 0, 2, 0), "mm"))

#add all legend together
legend <- plot_grid(get_legend(age_group_bar), get_legend(region_bar), ncol = 1)
age_group_bar <- age_group_bar + theme(legend.position = "none")
region_bar <- region_bar + theme(legend.position = "none")
#combine all graph together
plot <- plot_grid(age_group_bar, region_bar, nCount_RNA_plot, nFeature_RNA_plot, type_porportion_plot_no_legend, niche_porportion_plot_no_legend, align = "v", ncol = 1, axis = "b", rel_heights = c(1,1,7,7,7,7))
plot_grid(legend, plot, ncol = 1, rel_heights = c(1, 12))
ggsave("results/QC_plot_all_samples.pdf", width = 2.5, height = 7.5, units = "in")

##############
#plot cell type markers
DefaultAssay(combined.filt) <- "SCT"
type_markers <- c("TFAP2C", "CRYAB", "TNC", "EOMES", "NRP1", "CUX2", "RORB", "IL1RAPL2", "ZNF804B", "NWD2", "TSHZ2", "FAM160A1", "GAD1", "PBX3", "ADARB2", "VIP", "CCK", "CXCL14", "NXPH1", "SST", "ZNF385D", "EGFR", "OLIG1", "GJA1", "PCDH15", "PLP1" ,"SPP1", "IGFBP7")
type_marker_list <- list()
for (i in type_markers) {
  p <- FeaturePlot_scCustom(combined.filt, reduction = "umap", features = i, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  type_marker_list[[i]] <- p
}
patchwork::wrap_plots(type_marker_list, nrow = 4)
ggsave(file="results/all_marker_list.png", width = 25, height = 15, units = "in")



