library(Seurat)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(future)

plan("multicore", workers = 12)
options(future.globals.maxSize = 128000 * 1024^2)


setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/20221216_GPC_differentiation")

mat <- fread(file = "GSE135827_GE_mat_raw_count_with_week_info.txt")
genes <- mat[,1][[1]]
mat = data.frame(mat[,-1], row.names=genes)
meta <- read.csv(file = "meta.csv")
rownames(meta) <- colnames(mat)
Shi2021 <- CreateSeuratObject(counts = mat, project = "Shi2021", meta.data=meta, names.field = 2, names.delim = "\\.")
head(Shi2021)
UMAPoriginal <- cbind(Shi2021$UMAP.X, Shi2021$UMAP.Y, Shi2021$UMAP.Z)
colnames(UMAPoriginal) <- paste0("UMAPoriginal_", 1:3)
Shi2021[["UMAPoriginal"]] <- CreateDimReducObject(embeddings = UMAPoriginal, key = "UMAPoriginal_", assay = DefaultAssay(Shi2021))
DimPlot(Shi2021, dims = c(1,2), reduction = "UMAPoriginal", pt.size = 0.5, group.by = "Major.types")
ggsave(file="Major.types_dimplot.pdf",width = 8, height = 6, units = "in")
#save data
saveRDS(Shi2021, file = "Shi2021.rds")

#RPCA integration
# split the dataset into a list
Shi2021.list <- SplitObject(Shi2021, split.by = "orig.ident")
# normalize and identify variable features for each dataset independently
Shi2021.list <- lapply(X = Shi2021.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = Shi2021.list)
Shi2021.list <- lapply(X = Shi2021.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
#integration
anchors <- FindIntegrationAnchors(object.list = Shi2021.list, anchor.features = features, reduction = "rpca")
Shi2021.combined <- IntegrateData(anchorset = anchors)
#save data
saveRDS(Shi2021.combined, "Shi2021/Shi2021_integrated.rds")

#remove non-GE cells and then re-run UMAP to find clusters
Shi2021_GE <- subset(Shi2021.combined, subset = Major.types %in% c("CGE", "GE progenitor", "LGE", "MGE", "OPC"))
# Run the standard workflow for visualization and clustering
DefaultAssay(Shi2021.combined) <- "integrated"
Shi2021_GE <- ScaleData(Shi2021_GE, verbose = FALSE)
Shi2021_GE <- RunPCA(Shi2021_GE, npcs = 30, verbose = FALSE)
Shi2021_GE <- RunUMAP(Shi2021_GE, reduction = "pca", dims = 1:30)

#plot key markers: c("DCX", "BEST3", "CRABP1","LHX6", "MAF", "NPY", "VIP", "FOXP1", "ISL1", "PENK", "FOXP2", "TSHZ1", "MEIS2", "PAX6", "NR2F2", "PROX1")
IN.markers.plot.list <- lapply(X = c("DCX", "BEST3", "LHX6", "CRABP1", "MAF", "NPY", "VIP", "FOXP1", "ISL1", "PENK", "FOXP2", "TSHZ1", "MEIS2", "PAX6", "NR2F2", "PROX1"), FUN = function(x) {
  FeaturePlot_scCustom(Shi2021_GE, reduction = "umap", features = x, raster = F) + NoLegend() + NoAxes() + labs(title = NULL) + theme(plot.margin = margin(1, 0, 0, 0, "in"))
})
patchwork::wrap_plots(IN.markers.plot.list, nrow = 2)
ggsave(filename = "Shi2021/IN.markers.plot.list.png", width = 28, height = 9)

#find cluster
DefaultAssay(Shi2021_GE) <- "integrated"
use.pcs = 1:30
Shi2021_GE <- FindNeighbors(Shi2021_GE, reduction="pca", dims = use.pcs)
Shi2021_GE <- FindClusters(
  object = Shi2021_GE, 
  resolution = seq(3,8), 
  verbose = TRUE
)

for (i in c(3:8)) {
  DimPlot_scCustom(Shi2021_GE, group.by = paste0("integrated_snn_res.",i), reduction = "umap", label = TRUE, label.size = 6) + theme(legend.text = element_text(size = 18))
  ggsave(file= paste0("Shi2021/umap_cluster_integrated_snn_res.",i,".png"),width = 12, height = 10, units = "in")
}

#save data
saveRDS(Shi2021_GE, "Shi2021/Shi2021_GE.rds")

#find cluster markers
Idents(Shi2021_GE) <- "integrated_snn_res.6"
Shi2021_GE_markers_integrated_snn_res.6 <- FindAllMarkers(Shi2021_GE, assay = "RNA")
write.csv(Shi2021_GE_markers_integrated_snn_res.6, "Shi2021/Shi2021_GE_markers_integrated_snn_res.6.csv")

#add annotation
annotation <- read.csv("Shi2021/annotation.csv")
annotation$integrated_snn_res.6 <-  as.factor(annotation$integrated_snn_res.6)
metadata <- Shi2021_GE[[]]
metadata_annotated <- left_join(metadata, annotation, by = "integrated_snn_res.6")
rownames(metadata_annotated) <- rownames(metadata)
Shi2021_GE <- AddMetaData(Shi2021_GE, metadata_annotated)

Shi2021_GE$origin <- case_when(
  Shi2021_GE$Major.types %in% c("GE progenitor", "OPC") ~ "GE",
  .default = as.character(Shi2021_GE$Major.types)
)

#re-order metadata
Shi2021_GE$origin <- factor(Shi2021_GE$origin, levels = c("GE", "MGE", "LGE", "CGE"))
Shi2021_GE$type <- factor(Shi2021_GE$type, levels = c("Dividing_MKI67/TOP2A", "RG_VIM/HES1", "IPC-GE_OLIG2/EGFR", "IPC-IN_BEST3/ASCL1", "IN-Newborn_DCX/BEST3", "MGE_CRABP1/MAF",
                                                      "MGE_LHX6/MAF", "MGE_LHX6/NPY", "IN-Newborn_DCX/VIP", "LGE_FOXP1/ISL1", "LGE_FOXP1/PENK", "LGE_FOXP2/TSHZ1", "LGE/CGE_MEIS2/PAX6",
                                                      "CGE_NR2F2/PROX1", "VMF_ZIC1/ZIC2", "OPC_PDGFRA/SOX10", "Cortex_EMX1/TFAP2C", "Hippocampus_WNT7B/ZIC1"))
Shi2021_GE$inferred_terminal_type_abbr. <- factor(Shi2021_GE$inferred_terminal_type_abbr., levels = c("Transient", "STR-IN/CTX-IN", "CTX-IN", "CTX-IN/BN", "STR-dSPN", "STR-iSPN", "OB-IN/BN-eSPN", "OB-IN/DWM-IN", "SP"))
Shi2021_GE$inferred_terminal_type <- factor(Shi2021_GE$inferred_terminal_type, levels = c("Transient", "Striatal interneuron/cortical interneuron", "Cortical interneuron", "Cortical interneuron/basal nuclei inhibitory neuron", "D1 strialtal projection neuron", "D2 strialtal projection neuron", "Olfactory bulb interneuron/eccentric spiny projection neurons", "Olfactory bulb interneuron/deep white matter interneuron", "Subpallial neuron"))

#plot new annotation
p_subclass <- DimPlot_scCustom(Shi2021_GE, group.by = "subclass", reduction = "umap")
ggsave(file="Shi2021/umap_subclass.png",width = 7, height = 5, units = "in")

p_type <- DimPlot_scCustom(Shi2021_GE, group.by = "type", reduction = "umap", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", "#3283FEFF", 
                                                                                             "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", "#2ED9FFFF", 
                                                                                             "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", "#C4451CFF", 
                                                                                             "#1C8356FF", "#85660DFF", "#B10DA1FF"))
ggsave(file="Shi2021/umap_type.png",width = 8, height = 5, units = "in")

p_terminal_type <- DimPlot_scCustom(Shi2021_GE, group.by = "inferred_terminal_type_abbr.", reduction = "umap", colors_use = c("#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#2ED9FFFF", 
                                                                                                                              "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#C4451CFF")) + ggtitle("terminal type")
ggsave(file="Shi2021/umap_terminal_type.png",width = 7, height = 5, units = "in")

p_terminal_type_full_name <- DimPlot_scCustom(Shi2021_GE, group.by = "inferred_terminal_type", reduction = "umap", colors_use = c("#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#2ED9FFFF", 
                                                                                                                                  "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#C4451CFF")) + ggtitle("terminal type")
ggsave(file="Shi2021/umap_terminal_type_full_name.png",width = 11, height = 5, units = "in")

p_origin <- DimPlot_scCustom(Shi2021_GE, group.by = "origin", reduction = "umap")
ggsave(file="Shi2021/umap_origin.png",width = 6, height = 5, units = "in")

p_origin + p_subclass + p_type + p_terminal_type
ggsave(file="Shi2021/umap_origin_subclass_type_terminal_type.png",width = 14, height = 10, units = "in")
#save data
saveRDS(Shi2021_GE, "Shi2021/Shi2021_annotated.rds")





