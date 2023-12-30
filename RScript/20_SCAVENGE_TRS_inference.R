##03 TRS inference by gchromVAR & SCAVENGE

library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(dplyr)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg38)
library(igraph)

disease <- "DISEASE"
wdir <- paste0("/path/to/output_dir/",disease)
setwd(wdir)
setwd(wdir)

peaks <- readRDS("/path/to/rds/atac/peaks.rds")
seu <- readRDS(/path/to/rds/atac/integrated.rds)
set.seed(9527)
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
peaks <- keepSeqlevels(peaks,seqlevels(peaks)[1:22],pruning.model='coarse')

# quantify counts in each peak
counts <- FeatureMatrix(fragments = Fragments(seu),features = peaks,cells = colnames(seu))
assay <- CreateChromatinAssay(counts = counts,fragments = Fragments(seu),annotation = annotation)

SE <- SummarizedExperiment(assays = list(counts = assay),rowData = peaks, colData = seu@meta.data)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)
SE_bg <- getBackgroundPeaks(SE, niterations=200)
lsi_mat <- Embeddings(object = seu[["integrated_lsi"]])[,1:40]
ukbb <- importBedScore(rowRanges(SE), "PP.abf.hg38.bed", colidx = 5)
SE_DEV <- computeWeightedDeviations(object = SE, weights = ukbb, background_peaks = SE_bg)
z_score_mat <- data.frame(SummarizedExperiment::colData(SE), z_score=t(SummarizedExperiment::assays(SE_DEV)[["z"]]) |> c())
mutualknn30 <- getmutualknn(lsimat = lsi_mat2, num_k = 30)
seed_idx <- seedindex(z_score = z_score_mat$z_score,percent_cut = 0.05)
scale_factor <- cal_scalefactor(z_score = z_score_mat$z_score, percent_cut = 0.01)
np_score <- randomWalk_sparse(intM = mutualknn30, queryCells = rownames(mutualknn30)[seed_idx], gamma = 0.05)
omit_idx <- np_score==0
sum(omit_idx)
mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
np_score <- np_score[!omit_idx]
seed_idx <- seed_idx[!omit_idx]
TRS <- capOutlierQuantile(x = np_score, q_ceiling = 0.95) |> max_min_scale()
TRS <- TRS * scale_factor
mono_mat <- data.frame(z_score_mat[!omit_idx, ], seed_idx, np_score, TRS)
write.table(mono_mat,"TRS",sep="\t",row.names = T)

## permutation

mono_permu <- get_sigcell_simple(knn_sparse_mat=mutualknn30, seed_idx=mono_mat$seed_idx, topseed_npscore=mono_mat$np_score, 
                                 permutation_times=1000,true_cell_significance=0.001, rda_output=FALSE, mycores=8, rw_gamma=0.05)

mono_mat2 <- data.frame(mono_mat, mono_permu)
write.table(mono_mat2,"permutation.TRS'),sep="\t",row.names=T,quote=F)