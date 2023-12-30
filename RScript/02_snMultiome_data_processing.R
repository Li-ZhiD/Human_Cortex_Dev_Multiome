library(Seurat)
library(Signac)
library(sctransform)
library(tidyverse)
library(future)
library(scCustomize)
library(cowplot)
library(patchwork)
library(viridis)


plan("multicore", workers = 12)
options(future.globals.maxSize = 200000 * 1024^2)
setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")

dat = readRDS("multiome_031123_filtered.rds")


#find clusters
dat <- FindClusters(
  object = dat, 
  resolution = seq(3,7,1), 
  graph.name = "wsnn",
  algorithm = 3,
  verbose = TRUE
)

dat$wsnn_res.3 <- factor(dat$wsnn_res.3, levels = 0:137)
dat$wsnn_res.4 <- factor(dat$wsnn_res.4, levels = 0:153)
dat$wsnn_res.5 <- factor(dat$wsnn_res.5, levels = 0:175)
dat$wsnn_res.6 <- factor(dat$wsnn_res.6, levels = 0:185)
dat$wsnn_res.7 <- factor(dat$wsnn_res.7, levels = 0:209)

#plot clusters
for (i in seq(3,7,1)) {
  set.seed(2)
  M = nrow(unique(dat[[]][paste0("wsnn_res.",i)]))
  color = grDevices::colors()[intersect(grep('gr(e|a)y', grDevices::colors(), invert = T), grep('light', grDevices::colors(), invert = T))]
  color = color[!(color %in% c("azure", "azure1", "white", "ivory", "ivory1", "ivory2", "snow", "snow1", "snow2"))]
  myColours_type = sample(color,M)
  DimPlot(object = dat, group.by=paste0("wsnn_res.",i), reduction = "wnn.umap", label = T, label.size = 3, cols = myColours_type, raster = F) + theme(legend.text = element_text(size = 12))
  ggsave(file= paste0("results/annotation/wsnn_res.",i,"_umap_filtered.png"),width = 20, height = 12, units = "in")
}

#plot individual clusters
Idents(dat) <- "wsnn_res.5"
for (i in 0:175) {
  Cluster_Highlight_Plot(dat, cluster_name = i, reduction = "wnn.umap", raster = FALSE)
  ggsave(file= paste0("results/annotation/clusters/",i,".png"), width = 10, height = 9, units = "in")
}

#plot known markers
DefaultAssay(dat) <- "SCT"
human_markers <- c('FOXG1', 'OTX2', 'EN1', 'IRX3', 'HOXB4', 'ADCYAP1', 'NHLH1','EOMES', 'PPP1R17', 'RYR2', 'TNC', 'DCX', 'SYT1', 'RBFOX3', 'RELN', 'GAD1', 'GAD2', 'SST', 'SATB2', 'BCL11B', 'CUX2', 'NRGN', 'SLC17A7', 'SLC17A6', 'TLE4', 'FOXP2', 'PCP4', 'RORB', 'ADRA1A',
                   'CCN2', 'HCRTR2', 'SEMA3E', 'SYT6', 'GFRA1', 'LRATD2', 'TRPC7', 'TLL1', 'CCDC168', 'PCDH8', 'PIK3C2G', 'NPNT', 'FAM160A1', 'OLFML2B', 'POU3F1', 'SYP', 'ADAM33', 'THEMIS', 'TSHZ2', 'PRSS12', 'PDGFRA', 'PLP1', 'OLIG2', 'BCAS1', 'MBP', 'EGFR', 'C1QA', 'AIF1', 'P2RY12', 'GFAP', 'S100B', 'AQP4', 'SPARCL1', 'ID1', 'ID3',
                   'CLDN5', 'PECAM1', 'ASPN', 'CD93', 'CSPG4', 'RGS5', 'PDGFRB', 'MKI67', 'ASCL1', 'DLX1', 'DLX2', 'DLX6', 'SOX2', 'PAX6', 'HES5', 'VIM', "NES", 'HOPX', 'CRYAB', 'PDGFD', 'FEZF2', 'ALDH1L1', 'GSX2', 'PROX1', 'NR2F2', 'ADARB2', 'LHX6', 'MAF', 'NKX2-1', 'ISL1', 'MEIS1', 'MEIS2',
                   'ADARB2', 'LAMP5', 'PAX6', 'VIP', 'WIF1', 'SV2C', 'NDNF', 'CALB2', 'SCGN', 'PVALB', 'TAC1', 'ETS1', 'BTBD11', 'ADAMTS5', 'TBR1', 'HMGB2', 'LUM', 'OGN', 'DCN', 'APOE', 'CD38', 'CD74', 'CCR2', 'MS4A7', 'MRC1', 'CUX1', 'SOX5', 'NR2F1', 'GRIN2B', 'GRIN1', 'NR4A2', 'TTR', 'FOLR1',
                   'SP8', 'ARX', 'ZEB2', 'NXPH1', 'ACKR3', 'ERBB4', 'EPHB3', 'LHX8', 'ETV1', 'PBX3', 'TSHZ1', 'GBX1', 'ZIC1', 'EBF1', 'BEST3', 'WNT8B', 'ZIC2', 'ZIC3', 'GBX2', 'SIX3', 'EMX1', 'FGF17', 'SHH', 'BARHL2', 'LHX2', 'DBX1', 'WNT3A', 'TP73', 'IRF7', 'IRF8', 'RSPO3',
                   'RSPO2', 'LMX1A', 'NEUROD1', 'EPHB1', 'UNC5D', 'OLIG1', 'GJA1', 'ALDOC', 'SLC1A2', 'NWD2', 'HS3ST4', 'GREM2', 'MEGF11', 'ADAMTSL1', 'MEGF11', 'SULF1', 'EGR2', 'NFATC2', 'AGT', 'CA2', 'HEPACAM', 'POU6F2', 'TFAP2C', 'MYT1L', 'POSTN', 'DLX6-AS1', 'ITGA2', 'FBXO32', 'DACH1', "TFAP2C", 'MOXD1', 'ITGA2', 'ID3', 'GRM3',
                   'UNC5C', 'NFIX', 'SOX10', 'FOXJ1', 'CHODL', 'NRP1', 'NEUROD6', 'TRHDE', 'SNCG', 'CNR1', 'ST8SIA5', 'LIFR', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'LINC00343', 'SMYD1', 'IL1RAPL2', 'BMPER', 'CCK', 'ZNF385D', 'MOB3B', 'GRIK1')
human_markers_to_plot <- human_markers[human_markers %in% rownames(dat)]
for (i in human_markers_to_plot) {
  FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, pt.size = 0.1, raster = FALSE)
  ggsave(file=paste0("results/annotation/markers_final/", i,".png"),width = 15, height = 12, units = "in")
}

#find all markers in each cluster
Idents(dat) <- "wsnn_res.5"
wsnn_res.5 <- FindAllMarkers(dat, assay = "SCT", recorrect_umi = FALSE)
write.csv(wsnn_res.5, "results/annotation/wsnn_res.5_markers.csv")
wsnn_res.5 <- read.csv("results/annotation/wsnn_res.5_markers.csv", row.names = 1)
#retain top 50 markers for each cluster
wsnn_res.5_filt <- wsnn_res.5 %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% slice_head(n=50)
write.csv(wsnn_res.5_filt, "results/annotation/wsnn_res.5_markers_filt.csv", row.names = F)

#plot cell type specific markers
#plot RG markers
DefaultAssay(dat) <- "SCT"
RG_markers_list <- list()
for (i in c("HES1", "TFAP2C", "CRYAB", "FBXO32", "TNC", "HOPX")) {
  p <- FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL) + theme(plot.margin = margin(1, 0, 0, 0, "in"))
  RG_markers_list[[i]] <- p
}

#plot EN markers
EN_markers_list <- list()
for (i in c("SLC17A6", "EOMES", "NRP1", "GLIS3", "CUX2", "RORB", "IL1RAPL2", "ZNF804B", "FEZF2", "NWD2", "TSHZ2", "SYT6", "FAM160A1")) {
  p <- FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL) + theme(plot.margin = margin(1, 0, 0, 0, "in"))
  EN_markers_list[[i]] <- p
}

#plot IN markers
IN_markers_list <- list()
for (i in c("GAD1", "MEIS2", "ADARB2", "VIP", "CCK", "LAMP5", "NXPH1", "SST", "PVALB")) {
  p <- FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL) + theme(plot.margin = margin(1, 0, 0, 0, "in"))
  IN_markers_list[[i]] <- p
}

#plot Glia markers
Glia_markers_list <- list()
for (i in c("EGFR", "MAP3K1", "AQP4", "MOXD1", "GRM3", "GFAP", "OLIG1", "SOX10", "BCAS1", "MBP")) {
  p <- FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL) + theme(plot.margin = margin(1, 0, 0, 0, "in"))
  Glia_markers_list[[i]] <- p
}

#plot other markers
Other_markers_list <- list()
for (i in c("RELN", "IRF8", "IGFBP7")) {
  p <- FeaturePlot_scCustom(dat, reduction = "wnn.umap", features = i, pt.size = 0.1, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL) + theme(plot.margin = margin(1, 0, 0, 0, "in"))
  Other_markers_list[[i]] <- p
}

all_marker_list <- c(RG_markers_list, EN_markers_list, IN_markers_list, Glia_markers_list, Other_markers_list)
patchwork::wrap_plots(all_marker_list, nrow = 7)
ggsave(file="results/annotation/markers_examples/all_marker_list.png", width = 24, height = 36, units = "in")



#save data
saveRDS(dat, "multiome_031123_filtered_with_clusters.rds")

#manual annotation of clusters
#import annotation file
annotation <- read.csv("annotation.csv")
annotation$wsnn_res.5 <-  as.factor(annotation$wsnn_res.5)
metadata <- dat[[]]
metadata_annotated <- left_join(metadata, annotation, by = "wsnn_res.5")
rownames(metadata_annotated) <- rownames(metadata)
dat <- AddMetaData(dat, metadata_annotated)


############
#re-level metadata
dat$Group <- factor(dat$Group, levels = c("First_trimester", "Second_trimester", "Third_trimester", "Infancy", "Adolescence"))
dat$region_summary <- factor(dat$region_summary, levels = c("General", "PFC", "V1"))
dat$type <- factor(dat$type, levels = c("RG-vRG", "RG-tRG", "RG-oRG", "IPC-EN", "EN-Newborn", "EN-IT-Immature", "EN-L2_3-IT",
                                        "EN-L4-IT", "EN-L5-IT", "EN-L6-IT", "EN-Non-IT-Immature", "EN-L5-ET", "EN-L5_6-NP",
                                        "EN-L6-CT", "EN-L6b", "IN-dLGE-Immature", "IN-CGE-Immature", "IN-CGE-VIP", "IN-CGE-SNCG", "IN-CGE-LAMP5", "IN-MGE-Immature", "IN-MGE-SST",
                                        "IN-MGE-PV", "IPC-Glia", "Astrocyte-Immature", "Astrocyte-Protoplasmic", "Astrocyte-Fibrous", "OPC",
                                        "Oligodendrocyte-Immature", "Oligodendrocyte", "Cajal-Retzius cell", "Microglia", "Vascular", "Unknown"))
dat$dataset <- factor(dat$dataset, levels = c("HDBR-15020-FB", "HDBR-14584-CTX", "HDBR-14834-TC", "HDBR-14831-CTX", "ARKFrozen-45-CTX", "ARKFrozen-20-CTX",
                                              "ARKFrozen-37-CTX", "ARKFrozen-18-PFC", "ARKFrozen-19-V1", "ARKFrozen-41-PFC-2", "ARKFrozen-41-V1", "ARKFrozen-1-PFC",
                                              "ARKFrozen-32-V1", "ARKFrozen-43-PFC", "ARKFrozen-8-V1", "GW27-2-7-18-PFC", "NIH-5900-BA17", "NIH-4267-BA10-2",
                                              "NIH-M1154-BA10-2", "NIH-M2837-BA10-2", "NIH-4373-BA9-3", "NIH-4373-BA17", "NIH-1325-BA10-2", "NIH-1325-BA17",
                                              "NIH-4458-BA9-3", "NIH-4458-BA17", "NIH-4392-BA9", "NIH-4392-BA17", "NIH-1671-BA10-2", "NIH-1671-BA17",
                                              "NIH-5554-BA9", "NIH-5554-BA17", "NIH-5162-BA9", "NIH-5162-BA17", "NIH-5376-BA9", "NIH-5376-BA17",
                                              "NIH-4341-BA9", "NIH-4341-BA17"))


#find markers in each cell type
Idents(dat) <- "type"
type_markers <- FindAllMarkers(dat, assay = "SCT", recorrect_umi = FALSE, only.pos = TRUE)
write.csv(type_markers, "results/annotation/type_markers.csv")
#retain top 50 markers for each cell type
type_markers_filt <- type_markers %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% slice_head(n=50)
write.csv(type_markers_filt, "results/annotation/type_markers_filt.csv", row.names = F)

#####################
#plot results
#UMAP figures in Fig. 1
group_plot <- DimPlot_scCustom(dat, group.by = "Group", reduction = "wnn.umap", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), pt.size = 0.1, raster = F)
ggsave(file="results/group_umap_filtered.png",width = 11, height = 9, units = "in")
DimPlot_scCustom(dat, group.by = "Group", reduction = "wnn.umap", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), pt.size = 1, raster = T)
ggsave(file="results/group_umap_filtered.pdf",width = 6, height = 5, units = "in")

region_plot <- DimPlot_scCustom(dat, group.by = "region_summary", reduction = "wnn.umap", colors_use = c("#009e73", "#ffa500", "#0072b2"), pt.size = 0.1, raster = FALSE)
ggsave(file="results/region_summary_umap_filtered.png",width = 10, height = 9, units = "in")
DimPlot_scCustom(dat, group.by = "region_summary", reduction = "wnn.umap", colors_use = c("#009e73", "#ffa500", "#0072b2"), pt.size = 1, raster = T)
ggsave(file="results/region_summary_umap_filtered.pdf",width = 6, height = 5, units = "in")

type_plot <- DimPlot_scCustom(dat, group.by = "type", reduction = "wnn.umap", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                             "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                             "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                             "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                             "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                             "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                             "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF"), pt.size = 0.1, raster = FALSE)
ggsave(file="results/type_umap_filtered.png",width = 15, height = 9, units = "in")
DimPlot_scCustom(dat, group.by = "type", reduction = "wnn.umap", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF"), pt.size = 1, raster = T)
ggsave(file="results/type_umap_filtered.pdf",width = 15, height = 9, units = "in")

type_plot + group_plot + region_plot
ggsave(file="results/type_group_region.png",width = 34, height = 9, units = "in")

###########
#UMAP figures in Extended Data Fig. 1
dat_tmp <- dat
dat_tmp$log10_nCount_RNA <- log10(dat_tmp$nCount_RNA)
dat_tmp$log10_nFeature_RNA <- log10(dat_tmp$nFeature_RNA)
dat_tmp$log10_atac_peak_region_fragments <- log10(dat_tmp$atac_peak_region_fragments)

dput(DiscretePalette_scCustomize(num_colors = 38, palette = "varibow"))
dataset_plot <- DimPlot_scCustom(dat_tmp, group.by = "dataset", reduction = "wnn.umap",
                                 colors_use = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                                                "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                                                "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                                                "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                                                "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                                                "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                                                "#FF739F", "#CC3D54"),
                                 pt.size = 0.1, raster = FALSE)
ggsave(file="results/dataset_umap_filtered.png",width = 9, height = 5, units = "in")

log10_nCount_RNA_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_nCount_RNA", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/log10_nCount_RNA_umap_filtered.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat_tmp, features = "log10_nCount_RNA", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 1, raster = T, colors_use = viridis_light_high, order = F)
ggsave(file="results/log10_nCount_RNA_umap_filtered.pdf",width = 6, height = 5, units = "in")

log10_nFeature_RNA_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_nFeature_RNA", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/log10_nFeature_RNA_umap_filtered.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat_tmp, features = "log10_nFeature_RNA", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 1, raster = T, colors_use = viridis_light_high, order = F)
ggsave(file="results/log10_nFeature_RNA_umap_filtered.pdf",width = 6, height = 5, units = "in")

log10_atac_peak_region_fragments_plot <- FeaturePlot_scCustom(dat_tmp, features = "log10_atac_peak_region_fragments", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/log10_atac_peak_region_fragments_umap_filtered.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat_tmp, features = "log10_atac_peak_region_fragments", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 1, raster = T, colors_use = viridis_light_high, order = F)
ggsave(file="results/log10_atac_peak_region_fragments_umap_filtered.pdf",width = 6, height = 5, units = "in")

TSS_enrichment_plot <- FeaturePlot_scCustom(dat_tmp, features = "TSS.enrichment", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F) + scale_color_gradientn(colors = viridis(n = 10, direction = 1), limits = c(1, 8))
ggsave(file="results/TSS_enrichment_umap_filtered.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat_tmp, features = "TSS.enrichment", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 1, raster = T, colors_use = viridis_light_high, order = F) + scale_color_gradientn(colors = viridis(n = 10, direction = 1), limits = c(1, 8))
ggsave(file="results/TSS_enrichment_umap_filtered.pdf",width = 6, height = 5, units = "in")

nucleosome_signal_plot <- FeaturePlot_scCustom(dat_tmp, features = "nucleosome_signal", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 0.1, raster = F, colors_use = viridis_light_high, order = F)
ggsave(file="results/nucleosome_signal_umap_filtered.png",width = 6, height = 5, units = "in")
FeaturePlot_scCustom(dat_tmp, features = "nucleosome_signal", reduction = "wnn.umap", na_cutoff = NULL, pt.size = 1, raster = T, colors_use = viridis_light_high, order = F)
ggsave(file="results/nucleosome_signal_umap_filtered.pdf",width = 6, height = 5, units = "in")

dataset_plot + log10_nCount_RNA_plot + log10_nFeature_RNA_plot + log10_atac_peak_region_fragments_plot + TSS_enrichment_plot + nucleosome_signal_plot + plot_layout(ncol = 2)
ggsave(file="results/QC_umap_plots.png",width = 18, height = 20, units = "in")

###########
#cell type proportions in Fig. 1
#cell types per stage
metadata <- dat[[]]

pt_all <- table(metadata$type, metadata$Group)
pt_all <- as.data.frame(pt_all)
ggplot(pt_all, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.5) +
  xlab("Age Group") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                               "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                               "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                               "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                               "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                               "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                               "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF")) +
  scale_x_discrete(limits = rev(levels(pt_all$Var2))) +
  theme(legend.title = element_blank()) +
  coord_flip()
ggsave("results/percentage_per_group.pdf", width = 10, height = 4, units = "in")

#cell types per stage in PFC
metadata_PFC <- metadata %>% filter(region_summary == "PFC")
pt_PFC <- table(metadata_PFC$type, metadata_PFC$Group)
pt_PFC <- as.data.frame(pt_PFC)
ggplot(pt_PFC, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.5) +
  xlab("Age Group") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                               "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                               "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                               "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                               "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                               "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                               "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF")) +
  scale_x_discrete(limits = rev(levels(pt_PFC$Var2))) +
  theme(legend.title = element_blank()) +
  coord_flip()
ggsave("results/percentage_per_group_PFC.pdf", width = 10, height = 4, units = "in")

#cell types per stage in V1
metadata_V1 <- metadata %>% filter(region_summary == "V1")
pt_V1 <- table(metadata_V1$type, metadata_V1$Group)
pt_V1 <- as.data.frame(pt_V1)
ggplot(pt_V1, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.5) +
  xlab("Age Group") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                               "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                               "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                               "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                               "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                               "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                               "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF")) +
  scale_x_discrete(limits = rev(levels(pt_V1$Var2))) +
  theme(legend.title = element_blank()) +
  coord_flip()
ggsave("results/percentage_per_group_V1.pdf", width = 10, height = 4, units = "in")


#########
#plot QC metrics violin plots in Extended Data Fig. 1
nCount_RNA_plot <- ggplot(metadata, aes(x = dataset, y = nCount_RNA)) + 
  geom_violin(aes(fill = dataset), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab("# UMIs") +
  scale_fill_manual(values = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                               "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                               "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                               "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                               "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                               "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                               "#FF739F", "#CC3D54")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/nCount_RNA_per_dataset.png", width = 10, height = 2, units = "in")

nFeature_RNA_plot <- ggplot(metadata, aes(x = dataset, y = nFeature_RNA)) + 
  geom_violin(aes(fill = dataset), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab("# genes") +
  scale_fill_manual(values = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                               "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                               "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                               "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                               "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                               "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                               "#FF739F", "#CC3D54")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/nFeature_RNA_per_dataset.png", width = 10, height = 2, units = "in")

nFragments_in_peaks_plot <- ggplot(metadata, aes(x = dataset, y = log10(atac_peak_region_fragments))) + 
  geom_violin(aes(fill = dataset), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab(expression("Log10(# fragments in ATAC peaks)")) +
  scale_fill_manual(values = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                               "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                               "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                               "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                               "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                               "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                               "#FF739F", "#CC3D54")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/nCount_ATAC_per_dataset.png", width = 10, height = 2, units = "in")

TSS_enrichment_plot <- ggplot(metadata, aes(x = dataset, y = TSS.enrichment)) + 
  geom_violin(aes(fill = dataset), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  coord_cartesian(ylim=c(0, 10)) +
  ylab("TSS enrichment") +
  scale_fill_manual(values = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                               "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                               "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                               "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                               "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                               "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                               "#FF739F", "#CC3D54")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/TSS_enrichment_per_dataset.png", width = 10, height = 2, units = "in")

nucleosome_signal_plot <- ggplot(metadata, aes(x = dataset, y = nucleosome_signal)) + 
  geom_violin(aes(fill = dataset), trim = T) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  ylab("Nucleosome signal") +
  scale_fill_manual(values = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                               "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                               "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                               "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                               "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                               "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                               "#FF739F", "#CC3D54")) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x=element_blank())
ggsave("results/neucleosome_signal_per_dataset.png", width = 10, height = 2, units = "in")

#plot cell composition/proportion in individual samples
pt_sample <- table(metadata$type, metadata$dataset)
pt_sample <- data.frame(pt_sample)
type_proportion_per_sample <-  metadata %>% group_by(dataset) %>% mutate(total_cell_number = n()) %>% ungroup() %>% group_by(region_summary, Group, dataset, total_cell_number, log2_age, type) %>% summarise(type_number = n()) %>% mutate(type_proportion = type_number/total_cell_number)
write.csv(type_proportion_per_sample, "results/type_proportion_per_sample.csv", row.names = T)

porportion_plot <- ggplot(pt_sample, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_col(position = position_fill(reverse = TRUE), width = 0.5) +
  xlab("Dataset") +
  ylab("Proportion") +
  theme_classic() +
  scale_fill_manual(values = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                               "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                               "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                               "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                               "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                               "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                               "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF")) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  guides(fill=guide_legend(ncol= 3))
ggsave("results/percentage_per_dataset.pdf", width = 10, height = 4, units = "in")
porportion_plot_no_legend <- porportion_plot + theme(legend.position = "none")


#add annotation bar
annotation_bar_info <- distinct(metadata[,c("dataset", "region_summary", "Group")])
annotation_bar_info$dataset <- factor(annotation_bar_info$dataset, levels = c("HDBR-15020-FB", "HDBR-14584-CTX", "HDBR-14834-TC", "HDBR-14831-CTX", "ARKFrozen-45-CTX", "ARKFrozen-20-CTX",
                                              "ARKFrozen-37-CTX", "ARKFrozen-18-PFC", "ARKFrozen-19-V1", "ARKFrozen-41-PFC-2", "ARKFrozen-41-V1", "ARKFrozen-1-PFC",
                                              "ARKFrozen-32-V1", "ARKFrozen-43-PFC", "ARKFrozen-8-V1", "GW27-2-7-18-PFC", "NIH-5900-BA17", "NIH-4267-BA10-2",
                                              "NIH-M1154-BA10-2", "NIH-M2837-BA10-2", "NIH-4373-BA9-3", "NIH-4373-BA17", "NIH-1325-BA10-2", "NIH-1325-BA17",
                                              "NIH-4458-BA9-3", "NIH-4458-BA17", "NIH-4392-BA9", "NIH-4392-BA17", "NIH-1671-BA10-2", "NIH-1671-BA17",
                                              "NIH-5554-BA9", "NIH-5554-BA17", "NIH-5162-BA9", "NIH-5162-BA17", "NIH-5376-BA9", "NIH-5376-BA17",
                                              "NIH-4341-BA9", "NIH-4341-BA17"))
annotation_bar_info$Group <- factor(annotation_bar_info$Group, levels = c("First_trimester", "Second_trimester", "Third_trimester", "Infancy", "Adolescence"))
annotation_bar_info$region_summary <- factor(annotation_bar_info$region_summary, levels = c("General", "PFC", "V1"))
age_group_bar <- ggplot(annotation_bar_info) +
  geom_bar(mapping = aes(x = dataset, y = 1, fill = Group), 
           stat = "identity", 
           width = 0.85) +
  scale_fill_manual(values = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8")) +
  theme_void() +
  theme(legend.direction="horizontal", panel.spacing.x = unit(1, "mm"), plot.margin = unit(c(0, 0, 2, 0), "mm"))

region_bar <- ggplot(annotation_bar_info) +
  geom_bar(mapping = aes(x = dataset, y = 1, fill = region_summary), 
           stat = "identity", 
           width = 0.85) +
  scale_fill_manual(values = c("#009e73", "#ffa500", "#0072b2")) +
  theme_void() +
  theme(legend.direction="horizontal", panel.spacing.x = unit(1, "mm"), plot.margin = unit(c(0, 0, 2, 0), "mm"))

#add all legend together
legend <- plot_grid(get_legend(age_group_bar), get_legend(region_bar), ncol = 1)
age_group_bar <- age_group_bar + theme(legend.position = "none")
region_bar <- region_bar + theme(legend.position = "none")
#combine all graph together
plot <- plot_grid(age_group_bar, region_bar, nCount_RNA_plot, nFeature_RNA_plot, nFragments_in_peaks_plot, TSS_enrichment_plot, nucleosome_signal_plot, porportion_plot_no_legend, align = "v", ncol = 1, axis = "b", rel_heights = c(1,1,7,7,7,7,7,14))
plot_grid(legend, plot, ncol = 1, rel_heights = c(1, 12))
ggsave("results/QC_plot_all_samples.pdf", width = 10, height = 13, units = "in")

##############
#recalculate dimension reductions based on RNA and ATAC only for Extended Data Fig. 1
dat[["umap"]] <- NULL
dat[["umap_gex_cc"]] <- NULL
dat[["umap_atac"]] <- NULL
#RNA only
dat <-  RunUMAP(dat, reduction = "pca", dims = 1:50, reduction.name = "gex.umap", reduction.key = "gexUMAP_", min.dist = 0.27, n.neighbors = 50)

dataset_umap_filtered_RNA_only <- DimPlot_scCustom(dat, group.by = "dataset", reduction = "gex.umap",
                                                   colors_use = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                                                                  "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                                                                  "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                                                                  "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                                                                  "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                                                                  "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                                                                  "#FF739F", "#CC3D54"),
                                                   pt.size = 0.1, raster = FALSE)
ggsave(file="results/dataset_umap_filtered_RNA_only.png",width = 9, height = 5, units = "in")

age_group_umap_filtered_RNA_only <- DimPlot_scCustom(dat, group.by = "Group", reduction = "gex.umap", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), pt.size = 0.1, raster = F)
ggsave(file="results/age_group_umap_filtered_RNA_only.png",width = 11, height = 9, units = "in")

region_umap_filtered_RNA_only <- DimPlot_scCustom(dat, group.by = "region_summary", reduction = "gex.umap", colors_use = c("#009e73", "#ffa500", "#0072b2"), pt.size = 0.1, raster = FALSE)
ggsave(file="results/region_umap_filtered_RNA_only.png",width = 10, height = 9, units = "in")

type_umap_filtered_RNA_only <- DimPlot_scCustom(dat, group.by = "type", reduction = "gex.umap", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                               "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                               "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                               "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                               "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                               "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                               "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF"), pt.size = 0.1, raster = FALSE)
ggsave(file="results/type_umap_filtered_RNA_only.png",width = 15, height = 9, units = "in")

dataset_umap_filtered_RNA_only <- dataset_umap_filtered_RNA_only + NoLegend() + NoAxes() + labs(title = NULL)
age_group_umap_filtered_RNA_only <- age_group_umap_filtered_RNA_only + NoLegend() + NoAxes() + labs(title = NULL)
region_umap_filtered_RNA_only <- region_umap_filtered_RNA_only + NoLegend() + NoAxes() + labs(title = NULL)
type_umap_filtered_RNA_only <- type_umap_filtered_RNA_only + NoLegend() + NoAxes() + labs(title = NULL)
dataset_umap_filtered_RNA_only + age_group_umap_filtered_RNA_only + region_umap_filtered_RNA_only + type_umap_filtered_RNA_only + plot_layout(nrow = 4)
ggsave(file="results/dataset_type_group_region_umap_filtered_RNA_only.png",width = 5, height = 20, units = "in")

#ATAC only
dat <-  RunUMAP(dat, reduction = "integrated_lsi", dims = 2:40, reduction.name = "atac.umap", reduction.key = "atacUMAP_", min.dist = 0.27, n.neighbors = 50)

dataset_umap_filtered_ATAC_only <- DimPlot_scCustom(dat, group.by = "dataset", reduction = "atac.umap",
                                                   colors_use = c("#FF7373", "#CC543D", "#994017", "#FF7900", "#CCA35C", "#99822E", 
                                                                  "#FFF426", "#B7CC00", "#839945", "#B4FF4D", "#68CC1F", "#289900", 
                                                                  "#82FF73", "#3DCC45", "#179932", "#00FF5E", "#5CCC97", "#2E9977", 
                                                                  "#26FFDD", "#00CCCC", "#458C99", "#4DC7FF", "#1F7ACC", "#003899", 
                                                                  "#7390FF", "#3D45CC", "#251799", "#4300FF", "#8B5CCC", "#6C2E99", 
                                                                  "#C626FF", "#B700CC", "#994595", "#FF4DD9", "#CC1F8C", "#990048", 
                                                                  "#FF739F", "#CC3D54"),
                                                   pt.size = 0.1, raster = FALSE)
ggsave(file="results/dataset_umap_filtered_ATAC_only.png",width = 9, height = 5, units = "in")

age_group_umap_filtered_ATAC_only <- DimPlot_scCustom(dat, group.by = "Group", reduction = "atac.umap", colors_use = c("#f0f921", "#fca636", "#e16462", "#b12a90", "#6a00a8"), pt.size = 0.1, raster = F)
ggsave(file="results/age_group_umap_filtered_ATAC_only.png",width = 11, height = 9, units = "in")

region_umap_filtered_ATAC_only <- DimPlot_scCustom(dat, group.by = "region_summary", reduction = "atac.umap", colors_use = c("#009e73", "#ffa500", "#0072b2"), pt.size = 0.1, raster = FALSE)
ggsave(file="results/region_umap_filtered_ATAC_only.png",width = 10, height = 9, units = "in")

type_umap_filtered_ATAC_only <- DimPlot_scCustom(dat, group.by = "type", reduction = "atac.umap", colors_use = c("#5A5156FF", "#F6222EFF", "#FE00FAFF", "#16FF32FF", 
                                                                                                                 "#3283FEFF", "#FEAF16FF", "#B00068FF", "#1CFFCEFF", "#90AD1CFF", 
                                                                                                                 "#2ED9FFFF", "#DEA0FDFF", "#AA0DFEFF", "#F8A19FFF", "#325A9BFF", 
                                                                                                                 "#C4451CFF", "#1C8356FF", "#85660DFF", "#B10DA1FF", "#FBE426FF", 
                                                                                                                 "#1CBE4FFF", "#FA0087FF", "#FC1CBFFF", "#F7E1A0FF", "#C075A6FF", 
                                                                                                                 "#782AB6FF", "#AAF400FF", "#BDCDFFFF", "#B5EFB5FF", "#7ED7D1FF",
                                                                                                                 "#1C7F93FF", "#D85FF7FF", "#683B79FF", "#3B00FBFF", "#E4E1E3FF"), pt.size = 0.1, raster = FALSE)
ggsave(file="results/type_umap_filtered_ATAC_only.png",width = 15, height = 9, units = "in")

dataset_umap_filtered_ATAC_only <- dataset_umap_filtered_ATAC_only + NoLegend() + NoAxes() + labs(title = NULL)
age_group_umap_filtered_ATAC_only <- age_group_umap_filtered_ATAC_only + NoLegend() + NoAxes() + labs(title = NULL)
region_umap_filtered_ATAC_only <- region_umap_filtered_ATAC_only + NoLegend() + NoAxes() + labs(title = NULL)
type_umap_filtered_ATAC_only <- type_umap_filtered_ATAC_only + NoLegend() + NoAxes() + labs(title = NULL)
dataset_umap_filtered_ATAC_only + age_group_umap_filtered_ATAC_only + region_umap_filtered_ATAC_only + type_umap_filtered_ATAC_only+ plot_layout(nrow = 4)
ggsave(file="results/dataset_type_group_region_umap_filtered_ATAC_only.png",width = 5, height = 20, units = "in")

#save data
saveRDS(dat, "multiome_20230730_filtered_annotated.rds")

#save cell-level metadata
meta <- dat[[]]
meta_filt <- meta[,c("ID", "dataset", "Estimated_postconceptional_age_in_days", "Group", "sex", "region_summary", "gex_barcode", "atac_barcode", "nCount_RNA", "nFeature_RNA", "atac_peak_region_fragments", "pct_reads_in_peaks", "TSS.enrichment", "nucleosome_signal", "SCR_score",  "S.Score", "G2M.Score", "class", "subclass", "type", "wsnn_res.5")]
colnames(meta_filt) <- c("Cell_ID", "Sample_ID", "Estimated_postconceptional_age_in_days", "Group", "Sex", "Region", "GEX_barcode", "ATAC_barcode", "nCount_RNA", "nFeature_RNA", "ATAC_fragments_in_peaks", "Percentage_reads_in_peaks", "TSS.enrichment", "Nucleosome_signal", "Scrublet_doublet_score", "S.Score", "G2M.Score", "Class", "Subclass", "Type", "Cluster")
write.csv(meta_filt, "results/cell_level_metadata.csv", row.names = F)



