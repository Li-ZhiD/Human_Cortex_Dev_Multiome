library(Seurat)
library(Signac)
library(tidyverse)
library(scCustomize)
library(patchwork)
library(future)

plan("multicore", workers = 10)
options(future.globals.maxSize = 200000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")
dat <- readRDS("multiome_20230730_filtered_annotated.rds")

#focus on glial cells
dat_glia <- subset(dat, subset = type %in% c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC", "Oligodendrocyte-Immature", "Oligodendrocyte"))

#WNN and UMAP
dat_glia <- FindMultiModalNeighbors(dat_glia, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:40), knn.graph.name = "wknn.Glia", snn.graph.name = "wsnn.Glia", weighted.nn.name = "weighted.nn.Glia")
dat_glia <- RunUMAP(dat_glia, nn.name = "weighted.nn.Glia", reduction.name = "wnn.Glia.umap2", reduction.key = "wnngliaumap2_", min.dist = 0.3, n.neighbors = 30, n.components = 2)

#plot metadata

type <- DimPlot_scCustom(dat_glia, group.by = "type", reduction = "wnn.Glia.umap2", colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#7ED7D1", "#AAF400", "#BDCDFF", "#782AB6", "#1C7F93", "#D85FF7", "#822E1C"), raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
type_plot_details <- ggplot_build(type)
panel$layout$panel_scales_x
#x limits are -14.7 to 10.6
panel$layout$panel_scales_y
#y limits are -11.1 to 13.6
region <- DimPlot_scCustom(dat_glia, group.by = "region_summary", reduction = "wnn.Glia.umap2", colors_use = c("#009e73", "#ffa500", "#0072b2"), raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
age <- DimPlot_scCustom(dat_glia, group.by = "Group", reduction = "wnn.Glia.umap2", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)

glial_markers <- c("TFAP2C", "ITGA2", "CD38", "CRYAB", "TNC", "HOPX", "F3", "OLIG1", "OLIG2", "EGFR", "PDGFRA", "SOX10", "MBP" ,"GJA1", "SPARCL1")
for (i in glial_markers) {
  p <- FeaturePlot_scCustom(dat_glia, reduction = "wnn.Glia.umap2", features = i, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  assign(i, p)
}
age + region + type +TFAP2C + CRYAB + HOPX + OLIG2 + EGFR + SPARCL1 + GJA1 + SOX10 + MBP + plot_layout(ncol = 3)
ggsave(file= "results/macroglial_lineages/all_final_umap.png", width = 20, height = 25, units = "in")

#plot in pdf for legend
DimPlot_scCustom(dat_glia, group.by = "type", reduction = "wnn.Glia.umap2", colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#7ED7D1", "#AAF400", "#BDCDFF", "#782AB6", "#1C7F93", "#D85FF7", "#822E1C"), raster = T)
ggsave(file= "results/macroglial_lineages/all_final_umap.pdf", width = 6, height = 4, units = "in")

#################
#late 2nd trimester (ARKFrozen-1-PFC, ARKFrozen-32-V1, ARKFrozen-43-PFC, ARKFrozen-8-V1)
dat_glia_late_2tri <- subset(dat_glia, subset = dataset %in% c("ARKFrozen-1-PFC", "ARKFrozen-32-V1", "ARKFrozen-43-PFC", "ARKFrozen-8-V1"))

#plot metadata
type <- DimPlot_scCustom(dat_glia_late_2tri, group.by = "type", reduction = "wnn.Glia.umap2", colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#7ED7D1", "#AAF400", "#BDCDFF", "#782AB6", "#1C7F93", "#D85FF7", "#822E1C"), raster = FALSE) + lims(x= c(-14.7, 10.6), y = c(-11.1, 13.6)) + NoLegend() + NoAxes() + labs(title = NULL)
region <- DimPlot_scCustom(dat_glia_late_2tri, group.by = "region_summary", reduction = "wnn.Glia.umap2", colors_use = c("#ffa500", "#0072b2"), raster = FALSE) + lims(x= c(-14.7, 10.6), y = c(-11.1, 13.6)) + NoLegend() + NoAxes() + labs(title = NULL)
age <- DimPlot_scCustom(dat_glia_late_2tri, group.by = "Group", reduction = "wnn.Glia.umap2", colors_use = c("#fca636"), raster = FALSE) + lims(x= c(-14.7, 10.6), y = c(-11.1, 13.6)) + NoLegend() + NoAxes() + labs(title = NULL)

glial_markers <- c("TFAP2C", "ITGA2", "CD38", "CRYAB", "TNC", "HOPX", "F3", "OLIG1", "OLIG2", "EGFR", "PDGFRA", "SOX10", "MBP" ,"GJA1", "SPARCL1")
for (i in glial_markers) {
  p <- FeaturePlot_scCustom(dat_glia_late_2tri, reduction = "wnn.Glia.umap2", features = i, raster = FALSE) + lims(x= c(-14.7, 10.6), y = c(-11.1, 13.6)) + NoLegend() + NoAxes() + labs(title = NULL)
  assign(i, p)
}
age + region + type +TFAP2C + CRYAB + HOPX + OLIG2 + EGFR + SPARCL1 + GJA1 + SOX10 + MBP  + plot_layout(ncol = 3)
ggsave(file= "results/macroglial_lineages/all_final_umap_late_2tri.png", width = 10, height = 12.5, units = "in")

#plot surface marker genes
glial_markers <- c("ITGA2", "CD38", "F3", "EGFR", "PDGFRA")
for (i in glial_markers) {
  p <- FeaturePlot_scCustom(dat_glia_late_2tri, reduction = "wnn.Glia.umap2", features = i, raster = FALSE) + lims(x= c(-14.7, 10.6), y = c(-11.1, 13.6)) + NoLegend() + NoAxes() + labs(title = NULL)
  assign(i, p)
}
ITGA2 + CD38 + F3 + EGFR + PDGFRA + plot_layout(ncol = 5)
ggsave(file= "results/macroglial_lineages/surface_markers_umap_late_2tri.png", width = 16.667, height = 3.125, units = "in")

# Create the stacked violin plot 
Idents(dat_glia_late_2tri) <- "type"
Stacked_VlnPlot(seurat_object = dat_glia_late_2tri, features = glial_markers, idents = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-Glia", "Astrocyte-Immature", "OPC", "Oligodendrocyte-Immature"), x_lab_rotate = TRUE, pt.size = 0.1, colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#7ED7D1", "#AAF400", "#BDCDFF", "#782AB6", "#1C7F93", "#D85FF7", "#822E1C"))
ggsave(file="results/macroglial_lineages/Stacked_VlnPlot_late_2tri.pdf",width = 6, height = 7, units = "in")

#save data
saveRDS(dat_glia, "results/macroglial_lineages/dat_glia.rds")
saveRDS(dat_glia_late_2tri, "results/macroglial_lineages/dat_glia_late_2tri.rds")

