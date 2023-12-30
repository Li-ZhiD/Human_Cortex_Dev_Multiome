library(Seurat)
library(sctransform)
library(tidyverse)
library(scCustomize)
library(future)
library(patchwork)

plan("multicore", workers = 4)
options(future.globals.maxSize = 200000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/20221216_GPC_differentiation")

data <- Read10X("aggr/GPC2/outs/count/filtered_feature_bc_matrix")
experiment.aggregate <- CreateSeuratObject(
  data,
  project = "GPC",
  min.cells = 20,
  min.features = 250,
  names.field = 2,
  names.delim = "\\-")

#add meta_data
metadata <- read.csv("meta_data.csv")
metadata$orig.ident <- as.factor(metadata$orig.ident)
old_metadata <- experiment.aggregate[[]]
new_metadata <- left_join(old_metadata, metadata)
rownames(new_metadata) <- rownames(old_metadata)
experiment.aggregate <- AddMetaData(experiment.aggregate, new_metadata[,4:6])

#add mito fraction as metadata
experiment.aggregate <- PercentageFeatureSet(experiment.aggregate, pattern = "^MT-", col.name = "percent.mito")
head(experiment.aggregate@meta.data, n = 20)

#plot some quality control metrics
Idents(experiment.aggregate) <- "Sample_ID"
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, pt.size = 0)
ggsave("results/qc.png", width = 8, height = 6)

#filter low quality cells
experiment.aggregate <- subset(experiment.aggregate, percent.mito <= 10 & nFeature_RNA >= 1000 & nFeature_RNA < 10000)
table(experiment.aggregate$Sample_ID)

#save original seurat object
saveRDS(experiment.aggregate, file="results/filtered_GPC.rds")

#pre-processing
experiment.aggregate <- NormalizeData(experiment.aggregate)
experiment.aggregate <- FindVariableFeatures(object = experiment.aggregate, selection.method = "vst")
all.genes <- rownames(experiment.aggregate)
experiment.aggregate <- ScaleData(object = experiment.aggregate, features = all.genes)
experiment.aggregate <- RunPCA(object = experiment.aggregate, features = VariableFeatures(object = experiment.aggregate))
ElbowPlot(experiment.aggregate, ndims = 50)
ggsave("results/ElbowPlot.png", width = 8, height = 6)

#use PC 1:30 for UMAP
use.pcs = 1:30
experiment.aggregate <- RunUMAP(experiment.aggregate, dims = use.pcs, n.neighbors = 30L)

#save results
saveRDS(experiment.aggregate, file="results/processed_GPC.rds")

#find clusters
use.pcs = 1:30
experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction="pca", dims = use.pcs)
experiment.aggregate <- FindClusters(
  object = experiment.aggregate, 
  resolution = seq(4,8,1), 
  verbose = TRUE
)

#annotation
#import annotation file
annotation <- read.csv("results/annotation.csv")
annotation$RNA_snn_res.6 <-  as.factor(annotation$RNA_snn_res.6)
metadata <- experiment.aggregate[[]]
metadata_annotated <- left_join(metadata, annotation, by = "RNA_snn_res.6")
rownames(metadata_annotated) <- rownames(metadata)
experiment.aggregate <- AddMetaData(experiment.aggregate, metadata_annotated[,14:17])

#save data
saveRDS(experiment.aggregate, "results/annotated_GPC.rds")

#remove low quality cells and re-run PCA and UMAP
experiment.aggregate.final <- subset(experiment.aggregate, subset = type != "Low quality")
experiment.aggregate.final <- RunPCA(experiment.aggregate.final)
experiment.aggregate.final <- RunUMAP(experiment.aggregate.final, reduction = "pca", dims = 1:30)

#re-order metadata
experiment.aggregate.final$Sample_ID <- factor(experiment.aggregate.final$Sample_ID, levels = c("GW22_tRG_DIV0", "GW22_tRG_DIV7", "GW22_tRG_DIV14", "GW22_oRG_DIV0", "GW22_oRG_DIV7", "GW22_oRG_DIV14", "GW22_IPC-Glia_DIV0", "GW22_IPC-Glia_DIV7", "GW22_IPC-Glia_DIV14"))
experiment.aggregate.final$Stage <- factor(experiment.aggregate.final$Stage, levels = c("DIV0", "DIV7", "DIV14"))
experiment.aggregate.final$type <- factor(experiment.aggregate.final$type, levels = c("Dividing", "RG", "Ependymal cell", "IPC-EN", "EN", "IPC-Glia", "Astrocyte", "OPC","IPC-IN", "IN"))

#plot metadata
Sample_ID <- DimPlot_scCustom(experiment.aggregate.final, group.by = "Sample_ID", reduction = "umap", colors_use = RColorBrewer::brewer.pal(9,"Set1")) + NoLegend() + NoAxes() + labs(title = NULL)
Stage <- DimPlot_scCustom(experiment.aggregate.final, group.by = "Stage", reduction = "umap", colors_use = c("#FFD166", "#EF476F", "#26547C")) + NoLegend() + NoAxes() + labs(title = NULL)
Initial_cell_type <-DimPlot_scCustom(experiment.aggregate.final, group.by = "Initial_cell_type", reduction = "umap", colors_use = c("#66a182", "#edae49", "#d1495b")) + NoLegend() + NoAxes() + labs(title = NULL)
type <- DimPlot_scCustom(experiment.aggregate.final, group.by = "type", reduction = "umap", colors_use = c("#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A")) + NoLegend() + NoAxes() + labs(title = NULL)
Sample_ID + Stage + Initial_cell_type + type + plot_layout(ncol = 4)
ggsave("results/umap_all_final.png", width = 40, height = 10, units = "in")

DimPlot_scCustom(experiment.aggregate.final, label = F, group.by = "type", split.by = "Sample_ID", split_seurat = TRUE, num_columns = 3, colors_use = c("#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A")) + theme(strip.text.x = element_blank()) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave("results/umap_type_by_sample.png", width = 15, height = 15, units = "in")

#plot marker genes
type_markers <- c("TFAP2C", "CRYAB", "HOPX", "FOXJ1", "EOMES", "RBFOX3", "OLIG2", "EGFR", "SPARCL1", "SOX10", "BEST3" ,"DLX5")
for (i in type_markers) {
  p <- FeaturePlot_scCustom(experiment.aggregate.final, reduction = "umap", features = i, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  type_marker_list[[i]] <- p
}
TFAP2C + CRYAB + HOPX + FOXJ1 + EOMES + RBFOX3 + OLIG2 + EGFR + SPARCL1 + SOX10 + BEST3 + DLX5 + plot_layout(ncol = 4)
ggsave(file= "results/umap_all_markers_final.png", width = 20, height = 15, units = "in")

#proportion of each cell type in each sample
metadata <- experiment.aggregate.final[[]]
pt <- table(metadata$type, metadata$Sample_ID)
pt <- as.data.frame(pt)
pt <- pt %>% separate(Var2, into = c("Age", "Initial_cell_type", "Stage"), sep = "_", remove = FALSE)
pt$Stage <- factor(pt$Stage, levels = c("DIV0", "DIV7", "DIV14"))

pt_tRG <-  pt %>% filter(Initial_cell_type == "tRG")
tRG_pt_plot <- ggplot(pt_tRG, aes(x = Stage, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  ylab("Proportion") +
  xlab("tRG") +
  scale_fill_manual(values = c("#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A")) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

pt_oRG <-  pt %>% filter(Initial_cell_type == "oRG")
oRG_pt_plot <- ggplot(pt_oRG, aes(x = Stage, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = c("#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A")) +
  xlab("oRG") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_blank(), legend.position = "none")

pt_IPC_Glia <-  pt %>% filter(Initial_cell_type == "IPC-Glia")
IPC_Glia_pt_plot <- ggplot(pt_IPC_Glia, aes(x = Stage, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  scale_fill_manual(values = c("#B15928", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FF7F00", "#CAB2D6", "#6A3D9A")) +
  xlab("IPC-Glia") +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.title.y = element_blank())
pt_plot <- tRG_pt_plot + oRG_pt_plot + IPC_Glia_pt_plot
ggsave("results/pt.pdf", width = 8, height = 6, units = "in")


#rename sample ID
experiment.aggregate.final$Sample_ID <- case_when(
  experiment.aggregate.final$Sample_ID == "GW22_tRG_DIV0" ~ "tRG_DIV0",
  experiment.aggregate.final$Sample_ID == "GW22_tRG_DIV7" ~ "tRG_DIV7",
  experiment.aggregate.final$Sample_ID == "GW22_tRG_DIV14" ~ "tRG_DIV14",
  experiment.aggregate.final$Sample_ID == "GW22_oRG_DIV0" ~ "oRG_DIV0",
  experiment.aggregate.final$Sample_ID == "GW22_oRG_DIV7" ~ "oRG_DIV7",
  experiment.aggregate.final$Sample_ID == "GW22_oRG_DIV14" ~ "oRG_DIV14",
  experiment.aggregate.final$Sample_ID == "GW22_IPC-Glia_DIV0" ~ "IPC-Glia_DIV0",
  experiment.aggregate.final$Sample_ID == "GW22_IPC-Glia_DIV7" ~ "IPC-Glia_DIV7",
  experiment.aggregate.final$Sample_ID == "GW22_IPC-Glia_DIV14" ~ "IPC-Glia_DIV14"
)

experiment.aggregate.final$Seeding_cell_type <- experiment.aggregate.final$Initial_cell_type
experiment.aggregate.final$RNA_snn_res.4 <- NULL
experiment.aggregate.final$RNA_snn_res.5 <- NULL
experiment.aggregate.final$RNA_snn_res.7 <- NULL
experiment.aggregate.final$RNA_snn_res.8 <- NULL
experiment.aggregate.final$seurat_clusters <- NULL

#save data
saveRDS(experiment.aggregate.final, "results/final_GPC.rds")

#save metadata
metadata <- experiment.aggregate.final[[]]
metadata$Cell_ID <- rownames(metadata)
metadata_filt <- metadata[,c("Cell_ID", "Sample_ID", "Seeding_cell_type", "Stage", "nCount_RNA", "nFeature_RNA", "percent.mito", "class", "subclass", "type", "RNA_snn_res.6")]
colnames(metadata_filt) <- c("Cell_ID", "Sample_ID", "Seeding_cell_type", "Stage", "nCount_RNA", "nFeature_RNA", "Percentage_in_mitochondrial_genes", "Class", "Subclass", "Type", "Cluster")
write.csv(metadata_filt, "results/cell_level_metadata.csv", row.names = F)


