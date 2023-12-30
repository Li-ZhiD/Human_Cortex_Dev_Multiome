library(Seurat)
library(Signac)
library(slingshot)
library(tidyverse)
library(scCustomize)
library(mclust)
library(future)
library(patchwork)

plan("multicore", workers = 12)
options(future.globals.maxSize = 2000000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")
dat <- readRDS("multiome_20230730_filtered_annotated.rds")

#focus on EN lineage
dat_EN_lineage <- subset(dat, subset = subclass %in% c("Radial glia", "IPC-EN", "Glutamatergic neuron"))
#WNN and UMAP
dat_EN_lineage <- FindMultiModalNeighbors(dat_EN_lineage, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:40), knn.graph.name = "wknn.EN", snn.graph.name = "wsnn.EN", weighted.nn.name = "weighted.nn.EN")
dat_EN_lineage <- RunUMAP(dat_EN_lineage, nn.name = "weighted.nn.EN", reduction.name = "wnn.EN.umap", reduction.key = "wnnenumap_", min.dist = 0.1, n.neighbors = 50L, n.components = 8)
dat_EN_lineage <- RunUMAP(dat_EN_lineage, nn.name = "weighted.nn.EN", reduction.name = "wnn.EN.umap2", reduction.key = "wnnenumap2_", min.dist = 0.1, n.neighbors = 50L, n.components = 2)

#plot metadata
DimPlot_scCustom(dat_EN_lineage, group.by = "subclass", reduction = "wnn.EN.umap2",raster = FALSE)
ggsave(file="results/slingshot/EN/subclass_umap.png",width = 8, height = 5, units = "in")
DimPlot_scCustom(dat_EN_lineage, group.by = "type", reduction = "wnn.EN.umap2",
                 colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C"), raster = FALSE)
ggsave(file="results/slingshot/EN/type_umap.png",width = 8, height = 5, units = "in")
DimPlot_scCustom(dat_EN_lineage, group.by = "dataset", reduction = "wnn.EN.umap2", raster = FALSE)
ggsave(file="results/slingshot/EN/dataset_umap_filtered.png",width = 15, height = 9, units = "in")
DimPlot_scCustom(dat_EN_lineage, group.by = "region_summary", reduction = "wnn.EN.umap2", colors_use = c("#009e73", "#ffa500", "#0072b2"), raster = FALSE)
ggsave(file="results/slingshot/EN/region_summary.png",width = 7, height = 5, units = "in")
DimPlot_scCustom(dat_EN_lineage, group.by = "Group", reduction = "wnn.EN.umap2", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), raster = FALSE)
ggsave(file="results/slingshot/EN/group_umap_filtered.png",width = 7, height = 5, units = "in")
FeaturePlot_scCustom(dat_EN_lineage, features = "log2_age", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file="results/slingshot/EN/log2_age_umap_filtered.png",width = 7, height = 5, units = "in")

DimPlot_scCustom(dat_EN_lineage, group.by = "wsnn_res.5", reduction = "wnn.EN.umap2", raster = FALSE)
ggsave(file="results/slingshot/EN/wsnn_res.5_umap_filtered.png",width = 7, height = 5, units = "in")


#save data
saveRDS(dat_EN_lineage, "results/slingshot/EN/dat_EN_lineage.rds")



#get different levels of clusters using mclust clusters
rd = Embeddings(dat_EN_lineage, reduction="wnn.EN.umap")
dat_EN_lineage_mclust <- dat_EN_lineage
for (i in c(23:28)) {
  cl <- Mclust(rd, G = i)$classification
  dat_EN_lineage_mclust <-  AddMetaData(dat_EN_lineage_mclust, metadata = cl, col.name = paste0("mclust",i))
  DimPlot_scCustom(dat_EN_lineage_mclust, group.by = paste0("mclust",i), reduction = "wnn.EN.umap2", label = TRUE,raster = FALSE)
  ggsave(file=paste0("results/slingshot/EN/mclust/mclust",i, "_umap.png"), width = 8, height = 5, units = "in")
}

#save data with mclust information
saveRDS(dat_EN_lineage_mclust, "results/slingshot/EN/dat_EN_lineage_mclust.rds")


#use mclust24 as clusters
#remove cluster 24 which is an outlier
dat_EN_lineage_mclust_filt <- subset(dat_EN_lineage_mclust, subset = mclust24 != "24")
rd = Embeddings(dat_EN_lineage_mclust_filt, reduction="wnn.EN.umap")
cl <- dat_EN_lineage_mclust_filt$mclust24
times <- dat_EN_lineage_mclust$log2_age
lin <- getLineages(rd, cl, start.clus = "9",
                   end.clus = c("6", "7", "2", "14", "11", "13", "23", "19", "8"),
                   times = times, use.median = FALSE, dist.method = "slingshot")

slingLineages(lin)
crv <- getCurves(lin, approx_points = 200)

#save crv
saveRDS(crv, "results/slingshot/EN/curve.rds")

#extract slingshot results
pseudotime_mtx <- slingPseudotime(crv)
pseudotime_avg <- slingAvgPseudotime(crv)
cell_weights <- slingCurveWeights(crv)
branchID <- slingBranchID(crv)


#embed curve in 2 dimensional umap
rd_2d <-  Embeddings(dat_EN_lineage_mclust_filt, reduction="wnn.EN.umap2")
crv_2d <- embedCurves(crv, rd_2d)

#save crv embeded in 2d umap
saveRDS(crv_2d, "results/slingshot/EN/curve_2d.rds")

#save results to seurat obj.
dat_EN_lineage_mclust_filt$average_pseudotime <- pseudotime_avg
colnames(pseudotime_mtx) <- c("EN.L6b_pseudotime", "EN.L6.CT_pseudotime", "EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L5.IT_pseudotime", "EN.L4.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.IT_pseudotime", "RG_pseudotime")
dat_EN_lineage_mclust_filt <- AddMetaData(dat_EN_lineage_mclust_filt, metadata = pseudotime_mtx, col.name = colnames(pseudotime_mtx))
colnames(cell_weights) <- c("EN.L6b_lineageWeight", "EN.L6.CT_lineageWeight", "EN.L2_3.IT_lineageWeight", "EN.L4.IT_V1_lineageWeight", "EN.L5.IT_lineageWeight", "EN.L4.IT_lineageWeight", "EN.L5_6.NP_EN.L5.ET_lineageWeight", "EN.L6.IT_lineageWeight", "RG_lineageWeight")
dat_EN_lineage_mclust_filt <- AddMetaData(dat_EN_lineage_mclust_filt, metadata = cell_weights, col.name = colnames(cell_weights))

#plot results
DimPlot_scCustom(dat_EN_lineage_mclust_filt, group.by = "mclust24", colors_use = DiscretePalette_scCustomize(num_colors = 23, palette = "alphabet2"), reduction = "wnn.EN.umap2", label = TRUE, raster = FALSE)
ggsave(file=paste0("results/slingshot/EN/mclust/mclust24_umap.png"), width = 8, height = 5, units = "in")

#lineage with mclust
png(filename="results/slingshot/EN/lineage_mclust24_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(rd_2d, col = DiscretePalette_scCustomize(num_colors = 23, palette = "alphabet2")[cl], asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()
#curve with mclust
png(filename="results/slingshot/EN/curve_mclust24_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(rd_2d, col = DiscretePalette_scCustomize(num_colors = 23, palette = "alphabet2")[cl], asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with type
p_type <-  DimPlot_scCustom(dat_EN_lineage_mclust_filt, group.by = "type",
                            colors_use = c("#5A5156", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
                                           "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C"),
                            reduction = "wnn.EN.umap2", raster = FALSE)
p_type <- ggplot_build(p_type)
col <- unlist(p_type$data[[1]]["colour"])
png(filename="results/slingshot/EN/lineage_type_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_type$data[[1]]["x"], p_type$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()
#curve with type
png(filename="results/slingshot/EN/curve_type_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_type$data[[1]]["x"], p_type$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()
#curve with type more transparent
png(filename="results/slingshot/EN/curve_type_umap_transparent.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_type$data[[1]]["x"], p_type$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = scales::alpha(col,0.01), asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with age group
p_age_group <-  DimPlot_scCustom(dat_EN_lineage_mclust_filt, group.by = "Group",
                                 colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"),
                                 reduction = "wnn.EN.umap2", raster = FALSE)
p_age_group <- ggplot_build(p_age_group)
col <- unlist(p_age_group$data[[1]]["colour"])
png(filename="results/slingshot/EN/lineage_age_group_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_age_group$data[[1]]["x"], p_age_group$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()
#curve with age group
png(filename="results/slingshot/EN/curve_age_group_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_age_group$data[[1]]["x"], p_age_group$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with region
p_region <-  DimPlot_scCustom(dat_EN_lineage_mclust_filt, group.by = "region_summary",
                                 colors_use = c("#009e73", "#ffa500", "#0072b2"),
                                 reduction = "wnn.EN.umap2", raster = FALSE)
p_region <- ggplot_build(p_region)
col <- unlist(p_region$data[[1]]["colour"])
png(filename="results/slingshot/EN/lineage_region_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_region$data[[1]]["x"], p_region$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()
#curve with age group
png(filename="results/slingshot/EN/curve_region_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_region$data[[1]]["x"], p_region$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#lineage with pseudotime
p_psueodtime <- FeaturePlot_scCustom(dat_EN_lineage_mclust_filt, features = "average_pseudotime", reduction = "wnn.EN.umap2", colors_use = viridis_light_high, na_cutoff = NULL, raster = T)
ggsave(file = "results/slingshot/EN/average_pseudotime.pdf",width = 7, height = 5, units = "in")
p_psueodtime <- ggplot_build(p_psueodtime)
col <- unlist(p_psueodtime$data[[1]]["colour"])
png(filename="results/slingshot/EN/lineage_pseudotime_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_psueodtime$data[[1]]["x"], p_psueodtime$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'lineages', show.constraints = TRUE)
dev.off()
#curve with pseudotime
png(file="results/slingshot/EN/curve_pseudotime_umap.png", width = 8, height = 6, units = 'in', res = 300)
plot(c(p_psueodtime$data[[1]]["x"], p_psueodtime$data[[1]]["y"]), xlab="wnnenumap2_1", ylab="wnnenumap2_2", col = col, asp = 1, pch = 16, cex = 0.3)
lines(SlingshotDataSet(crv_2d), lwd = 2, col = 'black', type = 'curve', show.constraints = TRUE)
dev.off()

#plot pseudotime in individual lineages
for (i in c("EN.L6b_pseudotime", "EN.L6.CT_pseudotime", "EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L5.IT_pseudotime", "EN.L4.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.IT_pseudotime", "RG_pseudotime")) {
  p <- FeaturePlot_scCustom(dat_EN_lineage_mclust_filt, features = i, reduction = "wnn.EN.umap2", colors_use = viridis_light_high, na_cutoff = NULL, raster = FALSE)
  ggsave(file=paste0("results/slingshot/EN/individual_lineages/", i, ".png"),width = 7, height = 5, units = "in")
}
#psuedotime of all lineages in the same plot
pseudotime_list <- list()
for (i in c("EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L4.IT_pseudotime", "EN.L5.IT_pseudotime", "EN.L6.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.CT_pseudotime", "EN.L6b_pseudotime", "RG_pseudotime")) {
  p <- FeaturePlot_scCustom(dat_EN_lineage_mclust_filt, features = i, reduction = "wnn.EN.umap2", colors_use = viridis_light_high, na_cutoff = NULL, raster = FALSE) + labs(title = NULL) + NoAxes() + NoLegend() + scale_colour_viridis_c(limits = c(0,26.48495), oob=scales::squish, na.value = "lightgrey")
  pseudotime_list[[i]] <- p
}
patchwork::wrap_plots(pseudotime_list, nrow = 2)
ggsave(file="results/slingshot/EN/individual_lineages/all_pseudotime.png",width = 20, height = 7, units = "in")

#highlight cells in individual lineages
colors <- c('#00429d', '#3e67ae', '#618fbf', '#85b7ce', '#b1dfdb', '#ffcab9', '#fd9291', '#e75d6f', '#8bc53f')
lineages <- c("EN.L2_3.IT", "EN.L4.IT_V1", "EN.L4.IT", "EN.L5.IT", "EN.L6.IT", "EN.L5_6.NP_EN.L5.ET", "EN.L6.CT", "EN.L6b", "RG")
for (i in 1:9) {
  Idents(dat_EN_lineage_mclust_filt) <- ifelse(is.na(dat_EN_lineage_mclust_filt[[paste0(lineages[i], "_pseudotime")]]), "ident.remove", "ident.keep")
  cells <- WhichCells(object = dat_EN_lineage_mclust_filt, ident = "ident.keep")
  cells <- list(cells)
  names(cells) <- lineages[i]
  p <- Cell_Highlight_Plot(dat_EN_lineage_mclust_filt, cells_highlight = cells, reduction = "wnn.EN.umap2", highlight_color = colors[i], raster = FALSE)
  assign(lineages[i], p)
  ggsave(file=paste0("results/slingshot/EN/individual_lineages/", lineages[i], "_assignment.png"),width = 7, height = 5, units = "in")
}

#save plot regarding "EN.L2_3.IT", "EN.L4.IT_V1", "EN.L4.IT", "EN.L5.IT"
EN.L2_3.IT + EN.L4.IT_V1 + EN.L4.IT + EN.L5.IT + plot_layout(ncol = 2)
ggsave(file = "results/slingshot/EN/individual_lineages/BP2.png", width = 11.5, height = 8, units = "in")

##########################
# plot the expression of markers from adult Human V1 study 
for (i in c("TXK", "FAT4", "CXCL13", "SV2C", "CCDC68", "SNCAIP", "ABI3BP", "HS3ST2")) {
  p <- FeaturePlot_scCustom(dat_EN_lineage_mclust_filt, features = i, reduction = "wnn.EN.umap2", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file=paste0("results/slingshot/EN/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}
TXK + FAT4 + CXCL13 + SV2C + CCDC68 + SNCAIP + ABI3BP + HS3ST2 + plot_layout(ncol = 4)
ggsave(file="results/slingshot/EN/markers/L4_IT_combined.png", width = 22, height = 10, units = "in")


#save data with slingshot results
saveRDS(dat_EN_lineage_mclust_filt, "results/slingshot/EN/dat_EN_lineage_mclust_slingshot.rds")

#save single cell level metadata
metadata <- data.frame(dat_EN_lineage_mclust_filt[[]])
write.csv(metadata, file = "results/slingshot/EN/EN_lineage_slingshot_results.csv",row.names = TRUE)

meta_filt <- metadata[,c("ID", "dataset", "Estimated_postconceptional_age_in_days", "Group", "sex", "region_summary", "class", "subclass", "type", "mclust24",  "EN.L2_3.IT_lineageWeight",
                         "EN.L4.IT_V1_lineageWeight", "EN.L4.IT_lineageWeight", "EN.L5.IT_lineageWeight", "EN.L6.IT_lineageWeight", "EN.L5_6.NP_EN.L5.ET_lineageWeight", "EN.L6.CT_lineageWeight", "EN.L6b_lineageWeight", "RG_lineageWeight",
                         "EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L4.IT_pseudotime", "EN.L5.IT_pseudotime", "EN.L6.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.CT_pseudotime", "EN.L6b_pseudotime", "RG_pseudotime", "average_pseudotime")]
colnames(meta_filt) <- c("Cell_ID", "Sample_ID", "Estimated_postconceptional_age_in_days", "Group", "Sex", "Region", "Class", "Subclass", "Type", "Mclust_cluster", "EN.L2_3.IT_lineageWeight",
                         "EN.L4.IT_V1_lineageWeight", "EN.L4.IT_lineageWeight", "EN.L5.IT_lineageWeight", "EN.L6.IT_lineageWeight", "EN.L5_6.NP_EN.L5.ET_lineageWeight", "EN.L6.CT_lineageWeight", "EN.L6b_lineageWeight", "RG_lineageWeight",
                         "EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L4.IT_pseudotime", "EN.L5.IT_pseudotime", "EN.L6.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.CT_pseudotime", "EN.L6b_pseudotime", "RG_pseudotime", "Average_pseudotime")
write.csv(meta_filt, file = "results/slingshot/EN/EN_lineage_slingshot_results_filt.csv", row.names = F)


