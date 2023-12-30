library(tradeSeq)
library(SingleCellExperiment)
library(slingshot)
library(Seurat)
library(Signac)
library(tidyverse)
library(scCustomize)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)


setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome")
#load data
crv_2d <- readRDS("results/slingshot/EN/curve_2d.rds")
dat_EN_lineage_mclust_slingshot <- readRDS("results/slingshot/EN/dat_EN_lineage_mclust_slingshot.rds")
dat_AUC <- readRDS("results/AUCell/seurat_obj_with_AUC_scores.rds")


pseudotime <- slingPseudotime(crv_2d, na = FALSE)
colnames(pseudotime) <- c("EN.L6b", "EN.L6.CT", "EN.L2_3.IT", "EN.L4.IT_V1", "EN.L5.IT", "EN.L4.IT", "EN.L5_6.NP_EN.L5.ET", "EN.L6.IT", "RG")
cellWeights <- slingCurveWeights(crv_2d)
colnames(cellWeights) <- c("EN.L6b", "EN.L6.CT", "EN.L2_3.IT", "EN.L4.IT_V1", "EN.L5.IT", "EN.L4.IT", "EN.L5_6.NP_EN.L5.ET", "EN.L6.IT", "RG")

AUC <- GetAssayData(dat_AUC, slot = "data")
AUC <- as.matrix(AUC)
AUC_EN <- AUC[,rownames(pseudotime)]

#add AUC socres to dat_EN_lineage_mclust_slingshot object
dat_EN_lineage_mclust_slingshot[["scenicplus_AUC"]] <- CreateAssayObject(counts = AUC_EN)
DefaultAssay(dat_EN_lineage_mclust_slingshot) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC <- dat_EN_lineage_mclust_slingshot
#save EN lineage data with AUC scores
saveRDS(dat_EN_lineage_mclust_slingshot_AUC, file = "results/tradeSeq/dat_EN_lineage_mclust_slingshot_AUC.rds")
rm(list = c("dat_AUC", "AUC", "AUC_EN", "dat_EN_lineage_mclust_slingshot"))

##############################
#find the location of knots
set.seed(1)
assignCells <- function(cellWeights) {
  if (is.null(dim(cellWeights))) {
    if (any(cellWeights == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      return(matrix(1, nrow = length(cellWeights), ncol = 1))
    }
  } else {
    if (any(rowSums(cellWeights) == 0)) {
      stop("Some cells have no positive cell weights.")
    } else {
      # normalize weights
      normWeights <- sweep(cellWeights, 1,
                           FUN = "/",
                           STATS = apply(cellWeights, 1, sum)
      )
      # sample weights
      wSamp <- apply(normWeights, 1, function(prob) {
        stats::rmultinom(n = 1, prob = prob, size = 1)
      })
      # If there is only one lineage, wSamp is a vector so we need to adjust for that
      if (is.null(dim(wSamp))) {
        wSamp <- matrix(wSamp, ncol = 1)
      } else {
        wSamp <- t(wSamp)
      }
      return(wSamp)
    }
  }
}
findKnots <- function(nknots, pseudotime, wSamp) {
  # Easier to recreate them all here than to pass them on
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("t",ii), pseudotime[,ii])
  }
  for (ii in seq_len(ncol(pseudotime))) {
    assign(paste0("l",ii),1*(wSamp[,ii] == 1))
  }

  # Get the times for the knots
  tAll <- c()
  for (ii in seq_len(nrow(pseudotime))) {
    tAll[ii] <- pseudotime[ii, which(as.logical(wSamp[ii,]))]
  }

  knotLocs <- stats::quantile(tAll, probs = (0:(nknots - 1)) / (nknots - 1))
  if (any(duplicated(knotLocs))) {
    # fix pathological case where cells can be squeezed on one pseudotime value.
    # take knots solely based on longest lineage
    knotLocs <- stats::quantile(t1[l1 == 1],
                                probs = (0:(nknots - 1)) / (nknots - 1))
    # if duplication still occurs, get average btw 2 points for dups.
    if (any(duplicated(knotLocs))) {
      dupId <- duplicated(knotLocs)
      # if it's the last knot, get duplicates from end and replace by mean
      if (max(which(dupId)) == length(knotLocs)) {
        dupId <- duplicated(knotLocs, fromLast = TRUE)
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      } else {
        knotLocs[dupId] <- mean(c(knotLocs[which(dupId) - 1],
                                  knotLocs[which(dupId) + 1]))
      }
    }
    # if this doesn't fix it, get evenly spaced knots with warning
    if (any(duplicated(knotLocs))) {
      knotLocs <- seq(min(tAll), max(tAll), length = nknots)
    }
  }

  maxT <- max(pseudotime[,1])
  if (ncol(pseudotime) > 1) {
    maxT <- c()
    # note that first lineage should correspond to the longest, hence the
    # 100% quantile end point is captured.
    for (jj in 2:ncol(pseudotime)) {
      maxT[jj - 1] <- max(get(paste0("t", jj))[get(paste0("l",jj)) == 1])
    }
  }
  # if max is already a knot we can remove that
  if (all(maxT %in% knotLocs)) {
    knots <- knotLocs
  } else {
    maxT <- maxT[!maxT %in% knotLocs]
    replaceId <- vapply(maxT, function(ll){
      which.min(abs(ll - knotLocs))
    }, FUN.VALUE = 1)
    knotLocs[replaceId] <- maxT
    if (!all(maxT %in% knotLocs)) {
      # if not all end points are at knots, return a warning, but keep
      # quantile spaced knots.
      warning(paste0("Impossible to place a knot at all endpoints.",
                     "Increase the number of knots to avoid this issue."))
    }
    knots <- knotLocs
  }

  # guarantees that first knot is 0 and last knot is maximum pseudotime.
  knots[1] <- min(tAll)
  knots[nknots] <- max(tAll)

  knotList <- lapply(seq_len(ncol(pseudotime)), function(i){
    knots
  })
  names(knotList) <- paste0("t", seq_len(ncol(pseudotime)))

  return(knotList)
}
wSamp <- assignCells(cellWeights)
nknots=6
knotList <- findKnots(nknots, pseudotime, wSamp = wSamp)
saveRDS(knotList, "results/tradeSeq/EN/knotList.rds")


######################
# register BiocParallel
BPPARAM <- BiocParallel::MulticoreParam(workers = 20)
#run fitGAM
sce_EN_beta <- fitGAM(counts = AUC_EN, pseudotime = pseudotime, cellWeights = cellWeights, offset = rep(0, ncol(AUC_EN)),
                      nknots = 6, verbose = TRUE, parallel = TRUE, BPPARAM=BPPARAM,
                      family = "betar")
saveRDS(sce_EN_beta, "results/tradeSeq/EN/sce_EN_beta.rds")


##########################
#find eRegulon clusters
##########################
ysmooth_tidy <- predictSmooth(models = sce_EN_beta, rownames(sce_EN_beta), nPoints = 100)
saveRDS(ysmooth_tidy, "results/tradeSeq/EN/clustering/ysmooth_tidy.rds")
ysmooth_tidy_filt <- ysmooth_tidy %>% filter(lineage %in% 1:8)
saveRDS(ysmooth_tidy_filt, "results/tradeSeq/EN/clustering/ysmooth_tidy_filt.rds")

ysmooth <- predictSmooth(models = sce_EN_beta, rownames(sce_EN_beta), nPoints = 100, tidy = FALSE)
saveRDS(ysmooth, "results/tradeSeq/EN/clustering/ysmooth.rds")
#remove RG lineage
ysmooth_filt <- ysmooth[,1:800]
saveRDS(ysmooth_filt, "results/tradeSeq/EN/clustering/ysmooth_filt.rds")
ysmooth_filt_scaled <- t(scale(t(ysmooth_filt)))
saveRDS(ysmooth_filt_scaled, "results/tradeSeq/EN/clustering/ysmooth_scaled.rds")

#get scaled values in tidy format
df <- data.frame(ysmooth_filt_scaled)
df$gene <- rownames(df)
ysmooth_filt_scaled_tidy = df %>% pivot_longer(cols = !gene,
                                               cols_vary = "fastest",
                                               names_prefix = "lineage",
                                               names_sep = "_",
                                               names_to = c("lineage", "timepoint"),
                                               values_to = "yhat_scaled")
ysmooth_filt_scaled_tidy$time <- ysmooth_tidy_filt$time
ysmooth_filt_scaled_tidy$yhat <- ysmooth_tidy_filt$yhat
saveRDS(ysmooth_filt_scaled_tidy, "results/tradeSeq/EN/clustering/ysmooth_filt_scaled_tidy.rds")

# Decide how many clusters to look at for kmeans clustering
n_clusters <- 30
# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)
set.seed(123)
# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(ysmooth_filt_scaled, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}
# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot
ggsave(file="results/tradeSeq/EN/clustering/kmeans_scree.png",width = 8, height = 5, units = "in")

#we will used k = 6 based on scree plot results
set.seed(123)
km.out <- kmeans(ysmooth_filt_scaled, centers = 6, nstart = 100)
table(km.out$cluster)
cluster <-  data.frame(km.out$cluster)
colnames(cluster) <- "cluster"
cluster$gene <- rownames(cluster)
ysmooth_filt_scaled_tidy_with_cluster <- ysmooth_filt_scaled_tidy %>% left_join(cluster)
ysmooth_filt_scaled_tidy_with_cluster <- ysmooth_filt_scaled_tidy_with_cluster %>% mutate(
  lineage_name = case_when(
    lineage == 1 ~ "EN.L6b",
    lineage == 2 ~ "EN.L6.CT",
    lineage == 3 ~ "EN.L2_3.IT",
    lineage == 4 ~"EN.L4.IT_V1",
    lineage == 5 ~ "EN.L5.IT",
    lineage == 6 ~ "EN.L4.IT",
    lineage == 7 ~ "EN.L5_6.NP_EN.L5.ET",
    lineage == 8 ~"EN.L6.IT",
  ))
ysmooth_filt_scaled_tidy_with_cluster$lineage_name <- factor(ysmooth_filt_scaled_tidy_with_cluster$lineage_name, levels = c("EN.L2_3.IT", "EN.L4.IT_V1", "EN.L4.IT", "EN.L5.IT", "EN.L6.IT", "EN.L5_6.NP_EN.L5.ET", "EN.L6.CT", "EN.L6b"))
ysmooth_filt_scaled_tidy_with_cluster <- ysmooth_filt_scaled_tidy_with_cluster %>% mutate(
  module = case_when(
    cluster == 1 ~ "module6",
    cluster == 2 ~ "module2",
    cluster == 3 ~ "module5",
    cluster == 4 ~ "module3",
    cluster == 5 ~ "module1",
    cluster == 6 ~ "module4",
  ))

saveRDS(ysmooth_filt_scaled_tidy_with_cluster, "results/tradeSeq/EN/clustering/ysmooth_filt_scaled_tidy_with_cluster.rds")

ysmooth_filt_scaled_tidy_with_cluster_mean <- ysmooth_filt_scaled_tidy_with_cluster %>% group_by(module, time, lineage_name) %>% summarise(mean = mean(yhat_scaled))

ggplot(ysmooth_filt_scaled_tidy_with_cluster) +
  geom_line(aes(x = time, y = yhat_scaled, group = interaction(gene, lineage_name), color = lineage_name), alpha = 0.04, lwd = 0.2) +
  geom_line(aes(x = time, y = mean, group = lineage_name, color = lineage_name), data = ysmooth_filt_scaled_tidy_with_cluster_mean, lwd = 1) +
  facet_wrap(vars(module), scales = "free") +
  theme_bw() +
  ylab("Scaled gene-based AUC") +
  xlab("Pseudotime") +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(color = "beige", fill="beige")) +
  coord_cartesian(ylim = c(-3, 5)) +
  scale_color_manual(values = c('#00429d', '#3e67ae', '#618fbf', '#85b7ce', '#b1dfdb', '#ffcab9', '#fd9291', '#e75d6f'))
ggsave(file="results/tradeSeq/EN/clustering/cluster_expression_patterns.pdf",width = 8, height = 6, units = "in")

#export eRegulon with module information
eregulon_with_module <- ysmooth_filt_scaled_tidy_with_cluster %>% select(gene, module) %>% unique()
colnames(eregulon_with_module) <- c("eRegulon", "Module_membership")
eregulon_with_module$eRegulon <- str_replace_all(eregulon_with_module$eRegulon, pattern = ",", "_")
write.csv(eregulon_with_module, "results/tradeSeq/EN/clustering/module_membership.csv", row.names = F)

#GO enrichment analysis of eRegulon target genes in each module
#import eRegulon networks
geneSets <- getGmt("results/AUCell/582.sceplus.gmt")
#change gene name to be the same as geneSets, replacing "," with "_"
ysmooth_filt_scaled_tidy_with_cluster$gene <- str_replace_all(ysmooth_filt_scaled_tidy_with_cluster$gene, ",", "_")
#gene high-confident target genes in six modules (target genes have to be identified in at least 8% of all eRegulons in that module)
module_genes_filt <- c()
for (module in c("module1", "module2", "module3", "module4", "module5", "module6")) {
  module_TF <- ysmooth_filt_scaled_tidy_with_cluster %>% filter(module == !!module) %>% pull(gene) %>% unique()
  all_genes <- c()
  for (TF in module_TF) {
    genes <- geneIds(geneSets[[TF]])
    all_genes <- c(all_genes, genes)
  }
  tab <- table(all_genes)
  module_genes_filt[[module]] <- names(tab)[tab > length(module_TF)*0.08]
}

#change all gene symbols to ENTREZ ID for enrichGO analysis
module_ENTREZID_filt <- c()
for (module in names(module_genes_filt)) {
  IDs <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                               keys = module_genes_filt[[module]],
                               column = c("ENTREZID"),
                               keytype = "SYMBOL")
  IDs <- IDs[!is.na(IDs)]
  module_ENTREZID_filt[[module]] <- IDs
}

#get all genes from all eRegulons as background
universe_name <- unique(unlist(geneIds(geneSets)))
universe <- AnnotationDbi::mapIds(org.Hs.eg.db, 
                                  keys = universe_name,
                                  column = c("ENTREZID"),
                                  keytype = "SYMBOL")
universe <- universe[!is.na(universe)]

#go analysis
res_GO <- compareCluster(geneCluster = module_ENTREZID_filt,
                         fun = "enrichGO",
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         minGSSize = 10,
                         maxGSSize = 500,
                         universe = universe,
                         readable = TRUE)
res_GO_simplified <- simplify(res_GO)

write.csv(res_GO@compareClusterResult, file = "results/tradeSeq/EN/clustering/GO_enrichment/res_GO.csv")
write.csv(res_GO_simplified@compareClusterResult, file = "results/tradeSeq/EN/clustering/GO_enrichment/res_GO_simplified.csv")

#plot results in dotplot
result <- res_GO_simplified@compareClusterResult
gsize <- as.numeric(sub("/\\d+$", "", as.character(result$GeneRatio)))
gcsize <- as.numeric(sub("^\\d+/", "", as.character(result$GeneRatio)))
result$GeneRatio = gsize/gcsize

top_pathways <- c("mitotic cell cycle process", "regulation of cell cycle process", "cell division",
                     "ERK1 and ERK2 cascade", "response to epidermal growth factor", "positive regulation of cell differentiation",
                     "cytoplasmic translation", "regulation of neuron projection development", "cell projection morphogenesis",
                     "neuron projection morphogenesis", "regulation of cell projection organization", "neuron migration",
                     "cell junction assembly", "synapse assembly", "neuron recognition",
                     "trans-synaptic signaling", "regulation of synaptic plasticity", "inorganic cation transmembrane transport")

result_filt <- result %>% filter(Description %in% top_pathways)
result_filt <- result_filt %>% mutate("-log10(p.adj)" = -log10(p.adjust))
result_filt$Description <- factor(result_filt$Description, levels = rev(top_pathways))
#draw dotplot
p <- ggplot(data = result_filt, aes(x = Cluster, y = Description, size = GeneRatio, color = `-log10(p.adj)`)) +
  geom_point() + 
  scale_colour_viridis_c(limits = c(min(result_filt$`-log10(p.adj)`),8), oob=scales::squish)
p <- p + xlab("") + ylab("") + ggtitle("") +
  theme_dose(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.direction="horizontal")
p
ggsave(file="results/tradeSeq/EN/clustering/GO_enrichment/dotplot.pdf",width = 8.5, height = 5, units = "in")

#########################
#differential gene expression analysis
##########################
#global association test to find genes changes with pseudotime
assoRes <- associationTest(sce_EN_beta, global = TRUE, lineages = TRUE, l2fc = 1)
oassoRes <- assoRes[order(assoRes$waldStat, decreasing = TRUE, na.last = T),]
write.csv(oassoRes, "results/tradeSeq/EN/oassoRes.csv", row.names = T)

#Start vs End test to find genes different between RG and differentiated neurons
startvsendRes <- startVsEndTest(sce_EN_beta, global = TRUE, lineages = TRUE, l2fc = 1)
ostartvsendRes <- startvsendRes[order(startvsendRes$waldStat, decreasing = TRUE, na.last = T),]
write.csv(ostartvsendRes, "results/tradeSeq/EN/ostartvsendRes.csv", row.names = T)

#differentiated cell type differences
endRes <- diffEndTest(sce_EN_beta, global = TRUE, pairwise = TRUE, l2fc = 1)
oendRes <- endRes[order(endRes$waldStat, decreasing = TRUE, na.last = T),]
write.csv(oendRes, "results/tradeSeq/EN/oendRes.csv", row.names = T)

#identify cells between individual knots
#plot segments in all cells
metadata <- dat_EN_lineage_mclust_slingshot_AUC[[]]
metadata <- metadata %>%  mutate(
  knot_segment = case_when(
    average_pseudotime >= 23.732125 ~ 5,
    average_pseudotime >= 21.163459 ~ 4,
    average_pseudotime >= 18.988775 ~ 3,
    average_pseudotime >= 7.384903 ~ 2,
    .default = 1
  )
)
dat_EN_lineage_mclust_slingshot_AUC2 <- AddMetaData(dat_EN_lineage_mclust_slingshot_AUC, metadata)
DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC2, group.by = "knot_segment", reduction = "wnn.EN.umap2", raster = FALSE, colors_use =  c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"))
ggsave(file=paste0("results/tradeSeq/EN/plots/knots_segment.png"),width = 7, height = 5, units = "in")

plotGeneCount(curve = crv_2d, counts = counts(sce_EN_beta),
              clusters = dat_EN_lineage_mclust_slingshot_AUC2$knot_segment,
              models = sce_EN_beta) + scale_color_manual(values=DiscretePalette_scCustomize(num_colors = 5, palette = "varibow")[1:5])
ggsave(file="results/tradeSeq/EN/plots/curve_knots_segment.png",width = 8, height = 6, units = "in")

#plot segments in individual lineages
metadata <- dat_EN_lineage_mclust_slingshot_AUC[[]]
for (i in c("EN.L6b_pseudotime", "EN.L6.CT_pseudotime", "EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L5.IT_pseudotime", "EN.L4.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.IT_pseudotime", "RG_pseudotime")) {
  new_name <- paste0(i, "_knots")
  metadata <- metadata %>%  mutate(
    !!new_name := case_when(
      metadata[[i]] >= 23.732125 ~ 5,
      metadata[[i]] >= 21.163459 ~ 4,
      metadata[[i]] >= 18.988775 ~ 3,
      metadata[[i]] >= 7.384903 ~ 2,
      metadata[[i]] >= 0.000000 ~ 1,
      .default = NA
    )
  )
}
dat_EN_lineage_mclust_slingshot_AUC3 <- AddMetaData(dat_EN_lineage_mclust_slingshot_AUC2, metadata)
for (i in c("EN.L6b_pseudotime", "EN.L6.CT_pseudotime", "EN.L2_3.IT_pseudotime", "EN.L4.IT_V1_pseudotime", "EN.L5.IT_pseudotime", "EN.L4.IT_pseudotime", "EN.L5_6.NP_EN.L5.ET_pseudotime", "EN.L6.IT_pseudotime", "RG_pseudotime")) {
  knots <- paste0(paste0(i, "_knots"))
  DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = knots, reduction = "wnn.EN.umap2", raster = FALSE, colors_use =  DiscretePalette_scCustomize(num_colors = 5, palette = "varibow"))
  ggsave(file=paste0("results/tradeSeq/EN/plots/individual_lineages_knots/", i, ".png"), width = 7, height = 5, units = "in")
}


#use knots 2-4 (segment 2-3) for early differences between IT and non-IT neurons
earlyDERes_2_4 <- earlyDETest(sce_EN_beta, knots = c(2, 4), global = TRUE, pairwise  = TRUE, nPoints = 12, l2fc=log2(2))
write.csv(earlyDERes_2_4, "results/tradeSeq/EN/earlyDERes_2_4.csv", row.names = T)
#use knots 3-5 (segment 3-4) for early differences between (EN.L2_3.IT, EN.L4.IT_V1) vs (EN.L4.IT and EN.L5.IT) vs L6.IT, and EN.L6b vs EN.L6.CT vs EN.L5_6.NP_EN.L5.ET
earlyDERes_3_5 <- earlyDETest(sce_EN_beta, knots = c(3, 5), global = TRUE, pairwise  = TRUE, nPoints = 12, l2fc=log2(2))
write.csv(earlyDERes_3_5, "results/tradeSeq/EN/earlyDERes_3_5.csv", row.names = T)
#use diffEndTest results for differences between EN.L2_3.IT, EN.L4.IT_V1 and between EN.L4.IT and EN.L5.IT
#or use knots 4-6 (segment 4-5)
earlyDERes_4_6 <- earlyDETest(sce_EN_beta, knots = c(4, 6), global = TRUE, pairwise  = TRUE, nPoints = 12, l2fc=log2(2))
write.csv(earlyDERes_4_6, "results/tradeSeq/EN/earlyDERes_4_6.csv", row.names = T)

#save seurat obj with knots information
saveRDS(dat_EN_lineage_mclust_slingshot_AUC3, file = "results/tradeSeq/dat_EN_lineage_mclust_slingshot_AUC_knots.rds")


###################
#plot representative eGRN AUC values at branch points
###################
#summary 1 dotplot
#subset datasets into cells belonging to specific lineage branches and specific pseudotime segments
#branch point 1 IT vs non-IT
IT_2_4 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L2_3.IT_pseudotime_knots %in% c(2:3) |
                                                                   EN.L4.IT_V1_pseudotime_knots %in% c(2:3) |
                                                                   EN.L5.IT_pseudotime_knots %in% c(2:3) |
                                                                   EN.L4.IT_pseudotime_knots %in% c(2:3) |
                                                                   EN.L6.IT_pseudotime_knots %in% c(2:3)))
IT_2_4_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(IT_2_4)))
DefaultAssay(IT_2_4) <- "SCT"
IT_2_4_expression_avg <- rowMeans(as.matrix(GetAssayData(IT_2_4)))

non.IT_2_4 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L6b_pseudotime_knots %in% c(2:3) |
                                                                       EN.L6.CT_pseudotime_knots %in% c(2:3) |
                                                                       EN.L5_6.NP_EN.L5.ET_pseudotime_knots %in% c(2:3)))
non.IT_2_4_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(non.IT_2_4)))
DefaultAssay(non.IT_2_4) <- "SCT"
non.IT_2_4_expression_avg <- rowMeans(as.matrix(GetAssayData(non.IT_2_4)))

#branch point 2 EN.L2_4.IT vs EN.L4_5.IT vs EN.L6.IT
EN.L2_4.IT_3_5 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L2_3.IT_pseudotime_knots %in% c(3:4) |
                                                                   EN.L4.IT_V1_pseudotime_knots %in% c(3:4)))
EN.L2_4.IT_3_5_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L2_4.IT_3_5)))
DefaultAssay(EN.L2_4.IT_3_5) <- "SCT"
EN.L2_4.IT_3_5_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L2_4.IT_3_5)))

EN.L4_5.IT_3_5 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L5.IT_pseudotime_knots %in% c(3:4) |
                                                                           EN.L4.IT_pseudotime_knots %in% c(3:4)))
EN.L4_5.IT_3_5_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L4_5.IT_3_5)))
DefaultAssay(EN.L4_5.IT_3_5) <- "SCT"
EN.L4_5.IT_3_5_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L4_5.IT_3_5)))

EN.L6.IT_3_5 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L6.IT_pseudotime_knots %in% c(3:4)))
EN.L6.IT_3_5_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L6.IT_3_5)))
DefaultAssay(EN.L6.IT_3_5) <- "SCT"
EN.L6.IT_3_5_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L6.IT_3_5)))

#branch point 3 EN.L2_3.IT vs EN.L4.IT_V1
EN.L2_3.IT_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L2_3.IT_pseudotime_knots %in% c(4:5)))
EN.L2_3.IT_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L2_3.IT_4_6)))
DefaultAssay(EN.L2_3.IT_4_6) <- "SCT"
EN.L2_3.IT_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L2_3.IT_4_6)))

EN.L4.IT_V1_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L4.IT_V1_pseudotime_knots %in% c(4:5)))
EN.L4.IT_V1_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L4.IT_V1_4_6)))
DefaultAssay(EN.L4.IT_V1_4_6) <- "SCT"
EN.L4.IT_V1_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L4.IT_V1_4_6)))

#branch point 4 EN.L4.IT vs EN.L5.IT
EN.L4.IT_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L4.IT_pseudotime_knots %in% c(4:5)))
EN.L4.IT_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L4.IT_4_6)))
DefaultAssay(EN.L4.IT_4_6) <- "SCT"
EN.L4.IT_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L4.IT_4_6)))

EN.L5.IT_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L5.IT_pseudotime_knots %in% c(4:5)))
EN.L5.IT_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L5.IT_4_6)))
DefaultAssay(EN.L5.IT_4_6) <- "SCT"
EN.L5.IT_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L5.IT_4_6)))

#branch point 5 EN.L5_6.NP_EN.L5.ET vs EN.L6.CT vs EN.L6b
EN.L5_6.NP_EN.L5.ET_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L5_6.NP_EN.L5.ET_pseudotime_knots %in% c(4:5)))
EN.L5_6.NP_EN.L5.ET_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L5_6.NP_EN.L5.ET_4_6)))
DefaultAssay(EN.L5_6.NP_EN.L5.ET_4_6) <- "SCT"
EN.L5_6.NP_EN.L5.ET_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L5_6.NP_EN.L5.ET_4_6)))

EN.L6.CT_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L6.CT_pseudotime_knots %in% c(4:5)))
EN.L6.CT_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L6.CT_4_6)))
DefaultAssay(EN.L6.CT_4_6) <- "SCT"
EN.L6.CT_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L6.CT_4_6)))

EN.L6b_4_6 <- subset(dat_EN_lineage_mclust_slingshot_AUC3, subset = (EN.L6b_pseudotime_knots %in% c(4:5)))
EN.L6b_4_6_scenicplus_avg <- rowMeans(as.matrix(GetAssayData(EN.L6b_4_6)))
DefaultAssay(EN.L6b_4_6) <- "SCT"
EN.L6b_4_6_expression_avg <- rowMeans(as.matrix(GetAssayData(EN.L6b_4_6)))

#combine all scenicplus results
expression_avg <- cbind(IT_2_4_expression_avg, non.IT_2_4_expression_avg, EN.L2_4.IT_3_5_expression_avg, EN.L4_5.IT_3_5_expression_avg, EN.L6.IT_3_5_expression_avg,
                        EN.L2_3.IT_4_6_expression_avg, EN.L4.IT_V1_4_6_expression_avg, EN.L4.IT_4_6_expression_avg, EN.L5.IT_4_6_expression_avg,
                        EN.L5_6.NP_EN.L5.ET_4_6_expression_avg, EN.L6.CT_4_6_expression_avg, EN.L6b_4_6_expression_avg)
colnames(expression_avg) <- c("IT_2_4", "non.IT_2_4", "EN.L2_4.IT_3_5", "EN.L4_5.IT_3_5", "EN.L6.IT_3_5", "EN.L2_3.IT_4_6", "EN.L4.IT_V1_4_6", "EN.L4.IT_4_6", "EN.L5.IT_4_6", "EN.L5_6.NP_EN.L5.ET_4_6", "EN.L6.CT_4_6", "EN.L6b_4_6")
write.csv(expression_avg, "results/tradeSeq/EN/pseudobulk_expression_levels/expression_avg.csv",row.names = T)
#scale gene expression values


auc_avg <- cbind(IT_2_4_scenicplus_avg, non.IT_2_4_scenicplus_avg, EN.L2_4.IT_3_5_scenicplus_avg, EN.L4_5.IT_3_5_scenicplus_avg, EN.L6.IT_3_5_scenicplus_avg,
                        EN.L2_3.IT_4_6_scenicplus_avg, EN.L4.IT_V1_4_6_scenicplus_avg, EN.L4.IT_4_6_scenicplus_avg, EN.L5.IT_4_6_scenicplus_avg,
                        EN.L5_6.NP_EN.L5.ET_4_6_scenicplus_avg, EN.L6.CT_4_6_scenicplus_avg, EN.L6b_4_6_scenicplus_avg)
colnames(auc_avg) <- c("IT_2_4", "non.IT_2_4", "EN.L2_4.IT_3_5", "EN.L4_5.IT_3_5", "EN.L6.IT_3_5", "EN.L2_3.IT_4_6", "EN.L4.IT_V1_4_6", "EN.L4.IT_4_6", "EN.L5.IT_4_6", "EN.L5_6.NP_EN.L5.ET_4_6", "EN.L6.CT_4_6", "EN.L6b_4_6")
write.csv(auc_avg, "results/tradeSeq/EN/pseudobulk_expression_levels/auc_avg.csv",row.names = T)

#dotplot
TF_list <- c("POU3F1","FOXP1", "NFIB", "NFIA", "SOX4", "FEZF2", "SMAD3", "GLIS1", "CUX2", "BARX2", "NR1D2", "FOXP2","TCF4", "ZEB2",
             "EGR3", "ETV6", "POU6F2", "JDP2", "CUX1", "TWIST2", "HSF4","ETV1", "ETV5", "RXRG", "ERG", "BHLHE40","LCORL", "EBF3", "NR4A2")
expression_avg_filt <- expression_avg[TF_list,]
expression_avg_filt_scaled <- data.frame(t(scale(t(expression_avg_filt))))
expression_avg_filt_scaled$gene <- row.names(expression_avg_filt_scaled)
expression_avg_filt_scaled_long <- pivot_longer(expression_avg_filt_scaled, cols = 1:12, names_to = "lineage", values_to = "expression")

eGRN_list <- c("POU3F1,+,+","FOXP1,+,+", "NFIB,+,+", "NFIA,+,+", "SOX4,+,+", "FEZF2,extended,+,+", "SMAD3,+,+", "GLIS1,+,+", "CUX2,+,+", "BARX2,+,+", "NR1D2,extended,+,+", "FOXP2,+,+","TCF4,+,+", "ZEB2,+,+",
               "EGR3,+,+", "ETV6,+,+", "POU6F2,+,+", "JDP2,+,+", "CUX1,+,+", "TWIST2,extended,+,+", "HSF4,+,+","ETV1,+,+", "ETV5,+,+", "RXRG,+,+", "ERG,+,+", "BHLHE40,+,+","LCORL,+,+", "EBF3,+,+", "NR4A2,+,+")
auc_avg_filt <- auc_avg[eGRN_list,]
auc_avg_filt_scaled <- data.frame(t(scale(t(auc_avg_filt))))
auc_avg_filt_scaled$gene <- row.names(expression_avg_filt_scaled)
auc_avg_filt_scaled_long <- pivot_longer(auc_avg_filt_scaled, cols = 1:12, names_to = "lineage", values_to = "eGRN_auc")

combined <- left_join(expression_avg_filt_scaled_long, auc_avg_filt_scaled_long, by = c("gene", "lineage"))
combined$gene <- factor(combined$gene, levels = rev(TF_list))
combined$lineage <- factor(combined$lineage, levels = colnames(expression_avg))

ggplot(combined, aes(x=lineage, y=gene)) +
  geom_point(aes(size = eGRN_auc, fill = expression), color="black", shape=21) +
  scale_radius("Scaled eGRN AUC", range = c(-2,8), breaks = c(-1,0,1,2)) +
  scale_fill_viridis_c(option = "A", direction = 1, guide = guide_colorbar(ticks.colour = "black", frame.colour = "black"), name = "Sacled TF expression") +
  ylab("") + xlab("") +
  theme_classic() +
  theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"))
ggsave("results/tradeSeq/EN/plots/lineage_specific_TF_dotplot.pdf", width = 7, height = 8)

###################
#summary 2 representative eGRNs in UMAP plots
#branch point 1
branch_point1 <- ifelse((dat_EN_lineage_mclust_slingshot_AUC3$EN.L2_3.IT_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L4.IT_V1_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L5.IT_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L4.IT_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L6.IT_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L6b_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L6.CT_pseudotime_knots %in% c(2:3) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L5_6.NP_EN.L5.ET_pseudotime_knots %in% c(2:3)),
                        TRUE, FALSE)

#segments in branch point 1
dat_EN_lineage_mclust_slingshot_AUC3$branch_point1 <- dat_EN_lineage_mclust_slingshot_AUC3$knot_segment
dat_EN_lineage_mclust_slingshot_AUC3$branch_point1[!branch_point1] <- NA
DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point1", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF"), na.value = "lightgrey")
ggsave(file=paste0("results/tradeSeq/EN/plots/branch_point1.png"),width = 6, height = 5, units = "in")

#eGRN AUC
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC3$NFIB_AUC_branch_point1 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["NFIB,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$NFIB_AUC_branch_point1[!branch_point1] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NFIB_AUC_branch_point1", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/NFIB_AUC_branch_point1.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$POU3F1_AUC_branch_point1 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["POU3F1,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$POU3F1_AUC_branch_point1[!branch_point1] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "POU3F1_AUC_branch_point1", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/POU3F1_AUC_branch_point1.png"),width = 6, height = 5, units = "in")

#TF expression
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "SCT"
dat_EN_lineage_mclust_slingshot_AUC3$NFIB_expression_branch_point1 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["NFIB",]
dat_EN_lineage_mclust_slingshot_AUC3$NFIB_expression_branch_point1[!branch_point1] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NFIB_expression_branch_point1", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/NFIB_expression_branch_point1.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$POU3F1_expression_branch_point1 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["POU3F1",]
dat_EN_lineage_mclust_slingshot_AUC3$POU3F1_expression_branch_point1[!branch_point1] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "POU3F1_expression_branch_point1", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/POU3F1_expression_branch_point1.png"),width = 6, height = 5, units = "in")

###################
#branch point 2
branch_point2 <- ifelse((dat_EN_lineage_mclust_slingshot_AUC3$EN.L2_3.IT_pseudotime_knots %in% c(3:4) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L4.IT_V1_pseudotime_knots %in% c(3:4) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L5.IT_pseudotime_knots %in% c(3:4) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L4.IT_pseudotime_knots %in% c(3:4) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L6.IT_pseudotime_knots %in% c(3:4)),
                        TRUE, FALSE)

#segments in branch point 2
dat_EN_lineage_mclust_slingshot_AUC3$branch_point2 <- dat_EN_lineage_mclust_slingshot_AUC3$knot_segment
dat_EN_lineage_mclust_slingshot_AUC3$branch_point2[!branch_point2] <- NA
DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point2", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#AFCC3D", "#17994B", "#0066FF"), na.value = "lightgrey")
ggsave(file=paste0("results/tradeSeq/EN/plots/branch_point2.png"),width = 6, height = 5, units = "in")

#eGRN AUC
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC3$SMAD3_AUC_branch_point2 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["SMAD3,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$SMAD3_AUC_branch_point2[!branch_point2] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "SMAD3_AUC_branch_point2", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/SMAD3_AUC_branch_point2.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$TCF4_AUC_branch_point2 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["TCF4,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$TCF4_AUC_branch_point2[!branch_point2] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TCF4_AUC_branch_point2", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/TCF4_AUC_branch_point2.png"),width = 6, height = 5, units = "in")

#TF expression
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "SCT"
dat_EN_lineage_mclust_slingshot_AUC3$SMAD3_expression_branch_point2 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["SMAD3",]
dat_EN_lineage_mclust_slingshot_AUC3$SMAD3_expression_branch_point2[!branch_point2] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "SMAD3_expression_branch_point2", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/SMAD3_expression_branch_point2.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$TCF4_expression_branch_point2 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["TCF4",]
dat_EN_lineage_mclust_slingshot_AUC3$TCF4_expression_branch_point2[!branch_point2] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TCF4_expression_branch_point2", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/TCF4_expression_branch_point2.png"),width = 6, height = 5, units = "in")

###################
#branch point 3
branch_point3 <- ifelse((dat_EN_lineage_mclust_slingshot_AUC3$EN.L2_3.IT_pseudotime_knots %in% c(4:5) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L4.IT_V1_pseudotime_knots %in% c(4:5)),
                        TRUE, FALSE)

#segments in branch point 3
dat_EN_lineage_mclust_slingshot_AUC3$branch_point3 <- dat_EN_lineage_mclust_slingshot_AUC3$knot_segment
dat_EN_lineage_mclust_slingshot_AUC3$branch_point3[!branch_point3] <- NA
DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point3", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#17994B", "#0066FF", "#B65CCC"), na.value = "lightgrey")
ggsave(file=paste0("results/tradeSeq/EN/plots/branch_point3.png"),width = 6, height = 5, units = "in")

#eGRN AUC
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC3$ETV6_AUC_branch_point3 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["ETV6,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$ETV6_AUC_branch_point3[!branch_point3] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV6_AUC_branch_point3", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/ETV6_AUC_branch_point3.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$CUX1_AUC_branch_point3 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["CUX1,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$CUX1_AUC_branch_point3[!branch_point3] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "CUX1_AUC_branch_point3", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/CUX1_AUC_branch_point3.png"),width = 6, height = 5, units = "in")

#TF expression
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "SCT"
dat_EN_lineage_mclust_slingshot_AUC3$ETV6_expression_branch_point3 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["ETV6",]
dat_EN_lineage_mclust_slingshot_AUC3$ETV6_expression_branch_point3[!branch_point3] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV6_expression_branch_point3", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/ETV6_expression_branch_point3.png"),width = 6, height = 5, units = "in")

DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "SCT"
dat_EN_lineage_mclust_slingshot_AUC3$CUX1_expression_branch_point3 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["CUX1",]
dat_EN_lineage_mclust_slingshot_AUC3$CUX1_expression_branch_point3[!branch_point3] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "CUX1_expression_branch_point3", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/CUX1_expression_branch_point3.png"),width = 6, height = 5, units = "in")

###################
#branch point 4
branch_point4 <- ifelse((dat_EN_lineage_mclust_slingshot_AUC3$EN.L4.IT_pseudotime_knots %in% c(4:5) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L5.IT_pseudotime_knots %in% c(4:5)),
                        TRUE, FALSE)

#segments in branch point 4
dat_EN_lineage_mclust_slingshot_AUC3$branch_point4 <- dat_EN_lineage_mclust_slingshot_AUC3$knot_segment
dat_EN_lineage_mclust_slingshot_AUC3$branch_point4[!branch_point4] <- NA
DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point4", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#17994B", "#0066FF", "#B65CCC"), na.value = "lightgrey", drop = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/branch_point4.png"),width = 6, height = 5, units = "in")

#eGRN AUC
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC3$TWIST2_AUC_branch_point4 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["TWIST2,extended,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$TWIST2_AUC_branch_point4[!branch_point4] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TWIST2_AUC_branch_point4", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/TWIST2_AUC_branch_point4.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$ETV1_AUC_branch_point4 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["ETV1,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$ETV1_AUC_branch_point4[!branch_point4] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV1_AUC_branch_point4", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/ETV1_AUC_branch_point4.png"),width = 6, height = 5, units = "in")

#TF expression
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "SCT"
dat_EN_lineage_mclust_slingshot_AUC3$TWIST2_expression_branch_point4 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["TWIST2",]
dat_EN_lineage_mclust_slingshot_AUC3$TWIST2_expression_branch_point4[!branch_point4] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TWIST2_expression_branch_point4", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/TWIST2_expression_branch_point4.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$ETV1_expression_branch_point4 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["ETV1",]
dat_EN_lineage_mclust_slingshot_AUC3$ETV1_expression_branch_point4[!branch_point4] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV1_expression_branch_point4", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/ETV1_expression_branch_point4.png"),width = 6, height = 5, units = "in")

###################
#branch point 5
branch_point5 <- ifelse((dat_EN_lineage_mclust_slingshot_AUC3$EN.L5_6.NP_EN.L5.ET_pseudotime_knots %in% c(4:5) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L6.CT_pseudotime_knots %in% c(4:5) |
                           dat_EN_lineage_mclust_slingshot_AUC3$EN.L6b_pseudotime_knots %in% c(4:5)),
                        TRUE, FALSE)

#segments in branch point 5
dat_EN_lineage_mclust_slingshot_AUC3$branch_point5 <- dat_EN_lineage_mclust_slingshot_AUC3$knot_segment
dat_EN_lineage_mclust_slingshot_AUC3$branch_point5[!branch_point5] <- NA
DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point5", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#17994B", "#0066FF", "#B65CCC"), na.value = "lightgrey")
ggsave(file=paste0("results/tradeSeq/EN/plots/branch_point5.png"),width = 6, height = 5, units = "in")

#eGRN AUC
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC3$ERG_AUC_branch_point5 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["ERG,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$ERG_AUC_branch_point5[!branch_point5] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ERG_AUC_branch_point5", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/ERG_AUC_branch_point5.png"),width = 6, height = 5, units = "in")

DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "scenicplus_AUC"
dat_EN_lineage_mclust_slingshot_AUC3$NR4A2_AUC_branch_point5 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["NR4A2,+,+",]
dat_EN_lineage_mclust_slingshot_AUC3$NR4A2_AUC_branch_point5[!branch_point5] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NR4A2_AUC_branch_point5", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/NR4A2_AUC_branch_point5.png"),width = 6, height = 5, units = "in")

#TF expression
DefaultAssay(dat_EN_lineage_mclust_slingshot_AUC3) <- "SCT"
dat_EN_lineage_mclust_slingshot_AUC3$ERG_expression_branch_point5 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["ERG",]
dat_EN_lineage_mclust_slingshot_AUC3$ERG_expression_branch_point5[!branch_point5] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ERG_expression_branch_point5", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/ERG_expression_branch_point5.png"),width = 6, height = 5, units = "in")

dat_EN_lineage_mclust_slingshot_AUC3$NR4A2_expression_branch_point5 <- GetAssayData(dat_EN_lineage_mclust_slingshot_AUC3)["NR4A2",]
dat_EN_lineage_mclust_slingshot_AUC3$NR4A2_expression_branch_point5[!branch_point5] <- NA
FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NR4A2_expression_branch_point5", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE)
ggsave(file=paste0("results/tradeSeq/EN/plots/AUC_examples/NR4A2_expression_branch_point5.png"),width = 6, height = 5, units = "in")

#save EN lineage data with AUC scores and tradeSeq results
saveRDS(dat_EN_lineage_mclust_slingshot_AUC3, file = "results/tradeSeq/dat_EN_lineage_mclust_slingshot_AUC_tradeSeq.rds")


#plots in extended data fig related to Fig.3

library(patchwork)
branch_point1 <- DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point1", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF"), na.value = "lightgrey") + NoLegend() + NoAxes() + labs(title = NULL)
branch_point2 <- DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point2", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#AFCC3D", "#17994B", "#0066FF"), na.value = "lightgrey") + NoLegend() + NoAxes() + labs(title = NULL)
branch_point3 <- DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point3", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#17994B", "#0066FF", "#B65CCC"), na.value = "lightgrey") + NoLegend() + NoAxes() + labs(title = NULL)
branch_point4 <- DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point4", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#17994B", "#0066FF", "#B65CCC"), na.value = "lightgrey") + NoLegend() + NoAxes() + labs(title = NULL)
branch_point5 <- DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "branch_point5", reduction = "wnn.EN.umap2", colors_use = c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC"), raster = FALSE) + scale_colour_manual(values = c("#17994B", "#0066FF", "#B65CCC"), na.value = "lightgrey") + NoLegend() + NoAxes() + labs(title = NULL)

all <- DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "knot_segment", reduction = "wnn.EN.umap2", raster = FALSE, colors_use =  c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC")) + NoLegend() + NoAxes() + labs(title = NULL)

all + branch_point1 + branch_point2 + branch_point3 + branch_point4 + branch_point5 + plot_layout(nrow = 1)
ggsave(file=paste0("results/tradeSeq/EN/plots/all_branch_points.png"),width = 24, height = 3.5, units = "in")

DimPlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, group.by = "knot_segment", reduction = "wnn.EN.umap2", pt.size = 1, raster = TRUE, colors_use =  c("#FF7373", "#AFCC3D", "#17994B", "#0066FF", "#B65CCC")) + theme(legend.direction="horizontal")
ggsave(file="results/tradeSeq/EN/plots/legends_for_branch_points.pdf",width = 6, height = 5, units = "in")


NFIB_expression_branch_point1 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NFIB_expression_branch_point1", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
NFIB_AUC_branch_point1 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NFIB_AUC_branch_point1", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
TCF4_expression_branch_point2 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TCF4_expression_branch_point2", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
TCF4_AUC_branch_point2 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TCF4_AUC_branch_point2", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ETV6_expression_branch_point3 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV6_expression_branch_point3", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ETV6_AUC_branch_point3 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV6_AUC_branch_point3", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
TWIST2_expression_branch_point4 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TWIST2_expression_branch_point4", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
TWIST2_AUC_branch_point4 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "TWIST2_AUC_branch_point4", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ETV1_expression_branch_point4 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV1_expression_branch_point4", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ETV1_AUC_branch_point4 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ETV1_AUC_branch_point4", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ERG_expression_branch_point5 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ERG_expression_branch_point5", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
ERG_AUC_branch_point5 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "ERG_AUC_branch_point5", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
NR4A2_expression_branch_point5 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NR4A2_expression_branch_point5", reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)
NR4A2_AUC_branch_point5 <- FeaturePlot_scCustom(dat_EN_lineage_mclust_slingshot_AUC3, features = "NR4A2_AUC_branch_point5", colors_use = viridis_dark_high, reduction = "wnn.EN.umap2", na_cutoff = NULL, raster = FALSE) + NoLegend() + NoAxes() + labs(title = NULL)

NFIB_expression_branch_point1 + NFIB_AUC_branch_point1 + TCF4_expression_branch_point2 + TCF4_AUC_branch_point2 + ETV6_expression_branch_point3 + ETV6_AUC_branch_point3 +
  TWIST2_expression_branch_point4 + TWIST2_AUC_branch_point4+ETV1_expression_branch_point4 + ETV1_AUC_branch_point4 + plot_spacer() + plot_spacer() +
  ERG_expression_branch_point5+ERG_AUC_branch_point5+NR4A2_expression_branch_point5+NR4A2_AUC_branch_point5 + plot_layout(nrow = 3, ncol = 6)
ggsave(file="results/tradeSeq/EN/plots/AUC_examples/EDF_eGRN_examples.png",width = 22, height = 10, units = "in")

