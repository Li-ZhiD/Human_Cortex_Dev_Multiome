library(Seurat)
library(Signac)
library(future)
library(tidyverse)
library(scater)
library(scran)
library(edgeR)
library(glmmSeq)
library(scCustomize)
library(patchwork)
library(ggrepel)

plan("multisession", workers = 12)
options(future.globals.maxSize = 2000000 * 1024^2)
setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")

setwd("/kriegsteinlab/data3/LiWang/analysis/MERFISH/analysis")

#import primary data
primary_data <-  readRDS("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/multiome_20230730_filtered_annotated.rds")
#subset to focus on L4-IT
primary_L4_IT <- subset(primary_data, subset = type == "EN-L4-IT")
DefaultAssay(primary_L4_IT) <- "integrated"

#re-annotate primary data
primary_L4_IT <- RunUMAP(primary_L4_IT, reduction = "pca", dims = 1:50)
primary_L4_IT <- FindNeighbors(primary_L4_IT, dims = 1:50)
primary_L4_IT <- FindClusters(primary_L4_IT, resolution = 0.5)
DimPlot_scCustom(primary_L4_IT, reduction = "umap", group.by = "integrated_snn_res.0.5", raster = FALSE, label = T)
ggsave(filename = "results/L4_IT/primary_L4_IT_integrated_snn_res.0.5.png", width = 8, height = 6)

#annotate the lineages of individual clusters
annotation <- read.csv("results/L4_IT/annotation.csv")
annotation$integrated_snn_res.0.5 <-  as.factor(annotation$integrated_snn_res.0.5)
metadata <- primary_L4_IT[[]]
metadata_annotated <- left_join(metadata, annotation, by = "integrated_snn_res.0.5")
rownames(metadata_annotated) <- rownames(metadata)
primary_L4_IT <- AddMetaData(primary_L4_IT, metadata_annotated)

#save data
saveRDS(primary_L4_IT, "results/L4_IT/primary_L4_IT.rds")

#pseudobulk analysis
DefaultAssay(primary_L4_IT) <- "RNA"
primary_L4_IT[["SCT"]] <- NULL
primary_L4_IT[["ATAC"]] <- NULL
primary_L4_IT[["integrated"]] <- NULL
primary_L4_IT.sce <- as.SingleCellExperiment(primary_L4_IT)
primary_L4_IT.sce$ident <- NULL

#aggregate across cells
summed <- aggregateAcrossCells(primary_L4_IT.sce, ids = colData(primary_L4_IT.sce)[,c("subtype", "dataset")])

# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(summed), samples=colData(summed))

#remove pseudobulk cells consisting of less than 50 cells
discarded <- summed$ncells < 50
y <- y[,!discarded]
# Another typical step in bulk RNA-seq analyses is to remove genes that are lowly expressed. This reduces computational work, improves the accuracy of mean-variance trend modelling and decreases the severity of the multiple testing correction. Here, we use the filterByExpr() function from edgeR to remove genes that are not expressed above a log-CPM threshold in a minimum number of samples (determined from the size of the smallest treatment group in the experimental design).
keep <- filterByExpr(y, group=summed$subtype)
y <- y[keep,]
summary(keep)

#estiamate dispersion
sizeFactors <- calcNormFactors(y$counts, method="TMM")
disp <- setNames(edgeR::estimateDisp(y$counts)$tagwise.dispersion, rownames(y$counts))

#perform DE test using a generalized linear mixed model using glmmSeq package, LRT test was used here with a reduced model without the "subtype" covariate
results <- glmmSeq(modelFormula = ~ subtype + log2_age + (1 | dataset),
                   reduced = ~ log2_age + (1 | dataset),
                   countdata = y$counts,
                   metadata = y$samples,
                   dispersion = disp,
                   sizeFactors = sizeFactors,
                   cores = 12,
                   progress = TRUE)
results <- glmmQvals(results)
stats <- data.frame(summary(results))
colnames(stats)[c(14)] <- c("Q_LRT")
# add log2 fold change information from predicted values
predData <- results@predict[, 1:nrow(modelData)]
head(predData)
stats$l2fc <- log2(predData[,2])-log2(predData[,1])
stats$l2fc <- l2fc
write.csv(stats, "results/L4_IT/glmmseq_V1_vs_Common_res.csv")


#plot top markers in umap
DefaultAssay(primary_L4_IT) <- "SCT"
for (i in c("LINC02232", "CNR1", "KCNAB1", "CUX1", "KCNIP1", "LINC02055")) {
  p <- FeaturePlot_scCustom(primary_L4_IT, features = i, reduction = "umap", raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
  ggsave(file=paste0("results/L4_IT/markers/", i, ".png"), width = 7, height = 5, units = "in")
  assign(i, p)
}
primary_L4_IT_subtype_plot <- DimPlot_scCustom(primary_L4_IT, reduction = "umap", group.by = "subtype", raster = F, colors_use = c("#e63946", "#1d3557")) + NoLegend() + NoAxes() + labs(title = NULL)
primary_L4_IT_region_plot <- DimPlot_scCustom(primary_L4_IT, reduction = "umap", group.by = "region_summary", raster = F, colors_use = c("#009e73", "#ffa500", "#0072b2")) + NoLegend() + NoAxes() + labs(title = NULL)
primary_L4_IT_region_plot + primary_L4_IT_subtype_plot + CNR1 + CUX1 + LINC02232 + LINC02055 + KCNAB1 + KCNIP1 + plot_layout(ncol = 2)

ggsave(file="results/L4_IT/subtype_markers_combined.png", width = 10, height = 20, units = "in")

#volcono plot
stats$p.adj <- p.adjust(stats$P_LRT, method = 'BH')

stats$significant <- ifelse(stats$p.adj < 0.05, TRUE, FALSE)
stats <- stats %>% mutate(gene_type = case_when(l2fc > 1 & p.adj < 0.05 ~ "V1 biased",
                                                l2fc < -1 & p.adj < 0.05 ~ "Common biased",
                                                TRUE ~ "Other"))
stats$gene_type <- factor(stats$gene_type, levels = c("V1 biased", "Common biased", "Other"))
stats$symbol <- rownames(stats)
write.csv(stats, "results/L4_IT/glmmseq_V1_vs_Common_res_final.csv", row.names = F)

#select genes to highlight in volcano plot
genes_to_highlight <- stats %>%
  filter((l2fc > 2.3 & p.adj < 0.005) | (l2fc < -2.3 & p.adj < 0.005))
genes_to_highlight_specific <- genes_to_highlight[c("CUX1", "VAV3", "KCNIP1", "LINC02055", "CNTN5", "PLD5", "SLC35F4","CNR1", "LINC02232", "CBLN2", "TRPC5", "KCNAB1", "HS3ST2", "ARHGAP29"),]

# Add colour, size and alpha (transparency) to volcano plot --------------------
cols <- c("V1 biased" = "#1d3557", "Common biased" = "#e63946", "Other" = "grey")

ggplot(stats, aes(x = l2fc, y = -log10(p.adj))) +
  geom_point(mapping = aes(color = gene_type), alpha = 0.2, size = 0.5) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed") +
  scale_size_manual(values = 1) + # Modify point size
  scale_x_continuous(breaks = c(seq(-6, 6, 2)), # Modify x-axis tick intervals    
                     limits = c(-5.5, 5.5)) +
  geom_point(data = genes_to_highlight_specific,
             mapping = aes(color = gene_type),
             size = 2,
             alpha = 1) +
  scale_color_manual(values = cols) + # Modify point colour
  geom_text_repel(data = genes_to_highlight_specific, # Add labels last to appear as the top layer  
                  aes(label = symbol),
                  size = 6,
                  nudge_x = 0.2,
                  nudge_y = 0.4,
                  max.overlaps = 10,
                  min.segment.length = 0.21
  ) +
  ylim(0,4) +
  xlim(-8,8) +
  labs(x = "log2(fold change)",
       y = "-log10(adjusted P-value)",
       color = "Expression \ndifference") +
  theme_bw() + # Select theme with a white background  
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth= 0.5),    
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        #legend.text=element_text(size=15),
        #legend.title=element_text(size=18)
  )
ggsave(filename = "results/L4_IT/volcano.pdf", width = 7, height = 5, units = "in")
