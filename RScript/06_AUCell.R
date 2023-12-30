library(AUCell)
library(GSEABase)
library(Seurat)
library(Signac)
library(scCustomize)
library(tidyverse)
library(BiocParallel)
library(future)
plan("multicore", workers = 12)
options(future.globals.maxSize = 320000 * 1024^2)

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")

#import expression matrix
dat <- readRDS("multiome_20230730_filtered_annotated.rds")
DefaultAssay(object = dat) <- "RNA"
exprMatrix <- GetAssayData(object = dat, slot = "counts")
#import scenicplus eGRN gene sets
geneSets <- getGmt("results/AUCell/582.sceplus.gmt")
nGenes(geneSets)
#Run AUCell
cells_AUC <- AUCell_run(exprMatrix, geneSets)
saveRDS(cells_AUC, file = "results/AUCell/cells_AUC.rds")

#generate new Seurat object with AUC scores as features, replace "_" in gene set names by "." because Seurat does not allow "_" in feature name
AUC_mtx <- getAUC(cells_AUC)
eGRN_name <- rownames(AUC_mtx)
eGRN_name_new <- str_replace_all(eGRN_name, pattern = "_", ",")
rownames(AUC_mtx) <- eGRN_name_new
dat[["scenicplus_AUC"]] <- CreateAssayObject(counts = AUC_mtx)
DefaultAssay(dat) <- "scenicplus_AUC"
for (i in eGRN_name_new) {
  FeaturePlot_scCustom(dat, features = i, reduction = "wnn.umap", na_cutoff = NULL, raster = FALSE)
  ggsave(file=paste0("results/AUCell/plots/", i,".png"),width = 6.5, height = 5, units = "in")
}
#save data
saveRDS(dat, file = "results/AUCell/seurat_obj_with_AUC_scores.rds")

