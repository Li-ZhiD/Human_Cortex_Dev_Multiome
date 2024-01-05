library(Seurat)
library(CellChat)
library(reticulate)
wdir="/path/to/cellchat/output/"
seu_file="/path/to/rds/integrated.rds"
setwd(wdir)
seu=readRDS(seu_file)
Idents(seu)="type"
data.input <- GetAssayData(seu, assay = "SCT", slot = "data") 
labels <- Idents(seu)
group_of_interset <- unique(labels)
meta <- data.frame(group = labels, row.names = names(labels)) 
cell.use <- rownames(meta)[meta$group %in% group_of_interest]
data.input <- data.input[, cell.use]
meta <- data.frame(labels = meta[cell.use,], row.names = colnames(data.input))
meta$labels <- droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

df.net <- subsetCommunication(cellchat)
celltypes <- levels(cellchat@idents)[1:33]
groupSize=length(celltypes)

mat <- cellchat@net$weight[celltypes,celltypes]
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
