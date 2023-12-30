library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)

setwd('/path/to/workdir/')
seu <- readRDS("/path/to/rds")
DefaultAssay(seu) <- "ATAC"
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- as.logical(seqnames(granges(seu)) %in% main.chroms)
seu <- seu[keep.peaks, ]
pfm <- getMatrixSet(x = JASPAR2020,opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
seu <- AddMotifs(object = seu,genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm)

##chromvar
seu <- RunChromVAR(object = seu,genome = BSgenome.Hsapiens.UCSC.hg38)
saveRDS(seu, "/path/to/output/rds.rds")