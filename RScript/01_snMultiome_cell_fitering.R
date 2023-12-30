library(Seurat)
library(Signac)
library(tidyverse)
library(future)
library(scCustomize)
library(patchwork)

plan("multisession", workers = 4)
options(future.globals.maxSize = 2000000 * 1024^2)
setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")

dat <- readRDS("multiome_020323_joint_cc_filtered.rds")

#add some additional meta data
dat$log2_age <- log(dat$Estimated_postconceptional_age_in_days, base = 2)
region <- dat$Region
dat$region_summary <- case_when(
  region %in% c("Telencephalon", "Cortex", "Forebrain") ~ "General",
  region %in% c("PFC", "BA10", "BA9") ~ "PFC",
  region %in% c("V1", "BA17") ~ "V1",
)

#re-order metadata levels
dat$Group <- factor(dat$Group, levels = c("First_trimester", "Second_trimester", "Third_trimester", "Infancy", "Adolescence"))
dat$region_summary <- factor(dat$region_summary, levels = c("General", "PFC", "V1"))
dat$dataset <- factor(dat$dataset, levels = c("HDBR-15020-FB", "HDBR-14584-CTX", "HDBR-14834-TC", "HDBR-14831-CTX", "ARKFrozen-45-CTX", "ARKFrozen-20-CTX",
                                              "ARKFrozen-37-CTX", "ARKFrozen-18-PFC", "ARKFrozen-19-V1", "ARKFrozen-41-PFC-2", "ARKFrozen-41-V1", "ARKFrozen-1-PFC",
                                              "ARKFrozen-32-V1", "ARKFrozen-43-PFC", "ARKFrozen-8-V1", "GW27-2-7-18-PFC", "NIH-5900-BA17", "NIH-4267-BA10-2",
                                              "NIH-M1154-BA10-2", "NIH-M2837-BA10-2", "NIH-4373-BA9-3", "NIH-4373-BA17", "NIH-1325-BA10-2", "NIH-1325-BA17",
                                              "NIH-4458-BA9-3", "NIH-4458-BA17", "NIH-4392-BA9", "NIH-4392-BA17", "NIH-1671-BA10-2", "NIH-1671-BA17",
                                              "NIH-5554-BA9", "NIH-5554-BA17", "NIH-5162-BA9", "NIH-5162-BA17", "NIH-5376-BA9", "NIH-5376-BA17",
                                              "NIH-4341-BA9", "NIH-4341-BA17"))

#RunUMAP
dat <- FindMultiModalNeighbors(dat, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:40))
dat <- RunUMAP(dat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", n.neighbors = 30L)

#plot metadata
DimPlot_scCustom(dat, group.by = "dataset", reduction = "wnn.umap",
                 colors_use = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                                "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                                "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                                "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                                "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                                "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                                "#FF739F", "#CC3D54"),
                 pt.size = 0.1, raster = FALSE)
ggsave(file="results/pre_filtering/dataset_umap.png",width = 10, height = 6, units = "in")
DimPlot_scCustom(dat, group.by = "Group", reduction = "wnn.umap", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), pt.size = 0.1, raster = F)
ggsave(file="results/pre_filtering/age_group_umap.png",width = 7, height = 5, units = "in")
DimPlot_scCustom(dat, group.by = "region_summary", reduction = "wnn.umap", colors_use = c("#009e73", "#ffa500", "#0072b2"), pt.size = 1, raster = T)
ggsave(file="results/pre_filtering/region_summary_umap.png",width = 6, height = 5, units = "in",)
FeaturePlot_scCustom(dat, features = "S.Score", reduction = "wnn.umap", na_cutoff = NULL, raster = FALSE)
ggsave(file="results/pre_filtering/S.Score_umap.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat, features = "G2M.Score", reduction = "wnn.umap", na_cutoff = NULL, raster = FALSE)
ggsave(file="results/pre_filtering/G2M.score_umap.png",width = 6, height = 5, units = "in")

#plot all markers
DefaultAssay(dat) <- "SCT"
dat <- PrepSCTFindMarkers(dat)
human_markers <- c('FOXG1', 'OTX2', 'EN1', 'IRX3', 'HOXB4', 'ADCYAP1', 'NHLH1','EOMES', 'PPP1R17', 'RYR2', 'TNC', 'DCX', 'SYT1', 'RBFOX3', 'RELN', 'GAD1', 'GAD2', 'SST', 'SATB2', 'BCL11B', 'CUX2', 'NRGN', 'SLC17A7', 'SLC17A6', 'TLE4', 'FOXP2', 'PCP4', 'RORB', 'ADRA1A',
                   'CCN2', 'HCRTR2', 'SEMA3E', 'SYT6', 'GFRA1', 'LRATD2', 'TRPC7', 'TLL1', 'CCDC168', 'PCDH8', 'PIK3C2G', 'NPNT', 'FAM160A1', 'OLFML2B', 'POU3F1', 'SYP', 'ADAM33', 'THEMIS', 'TSHZ2', 'PRSS12', 'PDGFRA', 'PLP1', 'OLIG2', 'BCAS1', 'MBP', 'EGFR', 'C1QA', 'AIF1', 'P2RY12', 'GFAP', 'S100B', 'AQP4', 'SPARCL1', 'ID1', 'ID3',
                   'CLDN5', 'PECAM1', 'ASPN', 'CD93', 'CSPG4', 'RGS5', 'PDGFRB', 'MKI67', 'ASCL1', 'DLX1', 'DLX2', 'DLX6', 'SOX2', 'PAX6', 'HES5', 'VIM', "NES", 'HOPX', 'CRYAB', 'PDGFD', 'FEZF2', 'ALDH1L1', 'GSX2', 'PROX1', 'NR2F2', 'ADARB2', 'LHX6', 'MAF', 'NKX2-1', 'ISL1', 'MEIS1', 'MEIS2',
                   'ADARB2', 'LAMP5', 'PAX6', 'VIP', 'WIF1', 'SV2C', 'NDNF', 'CALB2', 'SCGN', 'PVALB', 'TAC1', 'ETS1', 'BTBD11', 'ADAMTS5', 'TBR1', 'HMGB2', 'LUM', 'OGN', 'DCN', 'APOE', 'CD38', 'CD74', 'CCR2', 'MS4A7', 'MRC1',
                   'CUX1', 'SOX5', 'NR2F1', 'GRIN2B', 'GRIN1', 'NR4A2', 'TTR', 'FOLR1','SP8', 'ARX', 'ZEB2', 'NXPH1', 'ACKR3', 'ERBB4', 'EPHB3', 'LHX8', 'ETV1', 'PBX3', 'TSHZ1', 'GBX1', 'ZIC1', 'EBF1', 'BEST3', 'WNT8B', 'ZIC2', 'ZIC3', 'GBX2', 'SIX3', 'EMX1', 'FGF17', 'SHH',
                   'BARHL2', 'LHX2', 'DBX1', 'WNT3A', 'TP73', 'IRF7', 'IRF8', 'RSPO3', 'RSPO2', 'LMX1A', 'NEUROD1', 'EPHB1', 'UNC5D', 'OLIG1', 'GJA1', 'ALDOC', 'SLC1A2', 'NWD2', 'HS3ST4', 'GREM2', 'ADAMTSL1', 'MEGF11', 'ADAMTSL1', 'MEGF11', 'SULF1', 'EGR2', 'NFATC2', 'AGT', 'CA2', 'HEPACAM', 'POU6F2',
                   'TFAP2C', 'MYT1L', 'POSTN', 'DLX6-AS1', 'ITGA2', 'FBXO32', 'DACH1', "TFAP2C", "MOXD1", "ITGA2", "ID3", "GRM3", "UNC5C", "NFIX", "SOX10", "FOXJ1")
human_markers_to_plot <- human_markers[human_markers %in% rownames(dat)]
for (i in human_markers_to_plot) {
  FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, na_cutoff = NA, raster = FALSE)
  ggsave(file=paste0("results/pre_filtering/markers/", i,".png"),width = 6, height = 5, units = "in")
}


#find clusters
dat <- FindClusters(
  object = dat, 
  resolution = seq(3,7,1), 
  graph.name = "wsnn",
  algorithm = 3,
  verbose = TRUE
)

dat$wsnn_res.3 <- factor(dat$wsnn_res.3, levels = c(0:138))
dat$wsnn_res.4 <- factor(dat$wsnn_res.4, levels = c(0:161))
dat$wsnn_res.5 <- factor(dat$wsnn_res.5, levels = c(0:178))
dat$wsnn_res.6 <- factor(dat$wsnn_res.6, levels = c(0:188))
dat$wsnn_res.7 <- factor(dat$wsnn_res.7, levels = c(0:218))

#plot individual clusters
options(max.overlaps = Inf)
DimPlot_scCustom(dat, group.by = "wsnn_res.7", reduction = "wnn.umap", pt.size = 0.1, raster = F, label = T, label.size = 3) + theme(legend.text = element_text(size = 12))
ggsave(file= "results/pre_filtering/wsnn_res.7_umap.png",width = 15, height = 8, units = "in")
#check the distribution of individual clusters
for (i in 0:218) {
  Cluster_Highlight_Plot(dat, cluster_name = i, reduction = "wnn.umap", raster = FALSE)
  ggsave(file= paste0("results/pre_filtering/clusters/",i,".png"), width = 10, height = 9, units = "in")
}

#plot QC metrics
dat_tmp <- dat
dat_tmp$log10_nCount_RNA <- log10(dat_tmp$nCount_RNA)
dat_tmp$log10_nFeature_RNA <- log10(dat_tmp$nFeature_RNA)
dat_tmp$log10_atac_peak_region_fragments <- log10(dat_tmp$atac_peak_region_fragments)
dat_tmp$log10_nCount_ATAC <- log10(dat_tmp$nCount_ATAC)

log10_nCount_RNA_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_nCount_RNA", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/pre_filtering/log10_nCount_RNA_umap.png",width = 6, height = 5, units = "in")

log10_nFeature_RNA_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_nFeature_RNA", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/pre_filtering/log10_nFeature_RNA_umap.png",width = 6, height = 5, units = "in")

log10_atac_peak_region_fragments_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_atac_peak_region_fragments", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/pre_filtering/log10_atac_peak_region_fragments_umap.png",width = 6, height = 5, units = "in")

log10_atac_peak_region_fragments_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_nCount_ATAC", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/pre_filtering/log10_nCount_ATAC_umap.png",width = 6, height = 5, units = "in")

TSS_enrichment_plot <- FeaturePlot_scCustom(dat_tmp, features = "TSS.enrichment", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F) + scale_color_gradientn(colors = viridis::viridis(n = 10, direction = 1), limits = c(1, 8))
ggsave(file="results/pre_filtering/TSS_enrichment_umap.png",width = 6, height = 5, units = "in")

nucleosome_signal_plot <- FeaturePlot_scCustom(dat_tmp, features = "pct_reads_in_peaks", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/pre_filtering/pct_reads_in_peaks_umap.png",width = 6, height = 5, units = "in")

nucleosome_signal_plot <- FeaturePlot_scCustom(dat_tmp, features = "nucleosome_signal", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/pre_filtering/nucleosome_signal_umap.png",width = 6, height = 5, units = "in")


#save data
saveRDS(dat, "multiome_031023_with_clusters.rds")

#############
#remove nuclei from diencephalon and striatum: clusters positive for ISL1, SIX3, OTX2, GBX2 but negative for FOXG1 and EMX1, and obvious additional doublet cells both NRGN and MBP positive
DefaultAssay(dat) <- "SCT"
for (i in c("FOXG1", "EMX1", "ISL1", "SIX3", "OTX2", "GBX2", "NRGN", "MBP", "PAX2", "PAX6")) {
  p <- FeaturePlot_scCustom(dat, features = i, reduction = "wnn.umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file=paste0("results/pre_filtering/markers_for_removing_clusters/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}
#plot annotated cell types from later analysis
meta_filt = read.csv( "results/cell_level_metadata.csv")
rownames(meta_filt) = meta_filt$Cell_ID
dat_tmp <- AddMetaData(dat, metadata = meta_filt)
dat_tmp$subclass_new_name <- ifelse(is.na(dat_tmp$Subclass), "Non-telencephalon, striatum, or doublets", dat_tmp$Subclass)
dat_tmp$subclass_new_name <- factor(dat_tmp$subclass_new_name, levels = c("Radial glia", "IPC-EN", "Glutamatergic neuron",
                                                                          "GABAergic neuron", "IPC-Glia", "Astrocyte",
                                                                          "OPC", "Oligodendrocyte", "Cajal-Retzius cell",
                                                                          "Microglia", "Vascular", "Unknown",
                                                                          "Non-telencephalon, striatum, or doublets"))
dat_tmp$Group <- factor(dat_tmp$Group, levels = c("First_trimester", "Second_trimester", "Third_trimester", "Infancy", "Adolescence"))


all_subtype <- DimPlot_scCustom(dat_tmp, group.by = "subclass_new_name", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                        "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                        "#2ED9FFFF", "#AA0DFEFF", "#E4E1E3FF", "#F8A19FFF"), reduction = "wnn.umap", pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = "results/pre_filtering/subclass.png", width = 7, height = 5, units = "in")
all_subtype_legend <- DimPlot_scCustom(dat_tmp, group.by = "subclass_new_name", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                        "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                        "#2ED9FFFF", "#AA0DFEFF", "#E4E1E3FF", "#F8A19FFF"), reduction = "wnn.umap", pt.size = 0.1, raster = T)
ggsave(file = "results/pre_filtering/subclass_legend.pdf", width = 7, height = 5, units = "in")

age_group <- DimPlot_scCustom(dat_tmp, group.by = "Group", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), reduction = "wnn.umap", pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)

#plot nuclei to remove
clusters_to_remove <- Cluster_Highlight_Plot(dat, cluster_name = c(17,191,181,81,161,107,155,204,85,176,35,44), highlight_color = "black", reduction = "wnn.umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ggsave(file = "results/pre_filtering/nuclei_to_remove.png", width = 10, height = 9, units = "in")
clusters_to_keep <- c(0:218)[!c(0:218) %in% c(17,191,181,81,161,107,155,204,85,176,35,44)]
dat_TC <- subset(dat, subset = wsnn_res.7 %in% clusters_to_keep)

all_subtype + clusters_to_remove + age_group + EMX1 + ISL1 + SIX3 + OTX2+ GBX2 + NRGN + MBP + plot_layout(ncol = 4)
ggsave(file = "results/pre_filtering/combined.png", width = 16, height = 12, units = "in")

#RunUMAP again on filtered data
dat_TC <- FindMultiModalNeighbors(dat_TC, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:50, 2:40))
dat_TC <- RunUMAP(dat_TC, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", min.dist = 0.27, n.neighbors = 50)

#plot metadata
DimPlot_scCustom(dat_TC, group.by = "dataset", reduction = "wnn.umap",
                 colors_use = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                                "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                                "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                                "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                                "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                                "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                                "#FF739F", "#CC3D54"),
                 pt.size = 0.1, raster = FALSE)
ggsave(file="results/annotation/dataset_umap_filtered.png",width = 10, height = 6, units = "in")
DimPlot_scCustom(dat_TC, group.by = "Group", reduction = "wnn.umap", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), pt.size = 0.1, raster = F)
ggsave(file="results/annotation/age_group_umap_filtered.png",width = 7, height = 5, units = "in")
DimPlot_scCustom(dat_TC, group.by = "region_summary", reduction = "wnn.umap", colors_use = c("#009e73", "#ffa500", "#0072b2"), pt.size = 0.1, raster = F)
ggsave(file="results/annotation/region_summary_umap_filtered.png",width = 6, height = 5, units = "in",)
FeaturePlot_scCustom(dat_TC, features = "S.Score", reduction = "wnn.umap", na_cutoff = NULL, raster = FALSE)
ggsave(file="results/annotation/S.Score_umap_filtered.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat_TC, features = "G2M.Score", reduction = "wnn.umap", na_cutoff = NULL, raster = FALSE)
ggsave(file="results/annotation/G2M.score_umap_filtered.png",width = 6, height = 5, units = "in")

#save data
saveRDS(dat_TC, "multiome_031123_filtered.rds")

