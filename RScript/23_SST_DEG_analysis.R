library("Seurat")
library("tidyverse")
library("scran")
library("qvalue")
library("clusterProfiler")
library("org.Hs.eg.db")
library("ReactomePA")
library("enrichplot")
library("DOSE")

setwd("/kriegsteinlab/data3/LiWang/analysis/human_cortex_multiome/SST")

dat_filt <- readRDS("sst_annotated_filt_agonist.rds")

DefaultAssay(dat_filt) <- "RNA"
dat_filt[["SCT"]] <- NULL

#convert to sce
sce.sst <- as.SingleCellExperiment(dat_filt)
sce.sst$ident = NULL
#create pseudobulk
summed <- aggregateAcrossCells(sce.sst, id = colData(sce.sst)[, c("type", "individual", "treatment")])

# Removing all pseudo-bulk samples with 'insufficient' cells.
summed.filt <- summed[,summed$ncells >= 10]

design <- model.matrix(~ individual + treatment, data = colData(summed.filt))
head(design)
summed.filt$individual <- factor(summed.filt$individual)

#focus on agonist
colData(summed.filt.agonist) <- droplevels(colData(summed.filt.agonist))


L.054264.de.results <- pseudoBulkDGE(summed.filt, 
                                     label = summed.filt$type,
                                     design = ~individual + treatment,
                                     coef = c("treatmentL-054264"),
                                     condition = summed.filt$treatment,
                                     method = "voom"
)
write.csv(L.054264.de.results[[1]], "results/DEG/L.054264.EN-IT-Immature.deg.csv")
write.csv(L.054264.de.results[[2]], "results/DEG/L.054264.EN-L5_6-NP.deg.csv")
write.csv(L.054264.de.results[[3]], "results/DEG/L.054264.EN-L6-CT.deg.csv")
write.csv(L.054264.de.results[[4]], "results/DEG/L.054264.EN-L6-IT.deg.csv")
write.csv(L.054264.de.results[[5]], "results/DEG/L.054264.EN-L6b.deg.csv")
write.csv(L.054264.de.results[[6]], "results/DEG/L.054264.EN-Newborn.deg.csv")
write.csv(L.054264.de.results[[9]], "results/DEG/L.054264.IPC-EN.degcsv")

Ostreotide.de.results <- pseudoBulkDGE(summed.filt, 
                                       label = summed.filt$type,
                                       design = ~individual + treatment,
                                       coef = c("treatmentOstreotide"),
                                       condition = summed.filt$treatment,
                                       method = "voom"
)
write.csv(Ostreotide.de.results[[1]], "results/DEG/Ostreotide.EN-IT-Immature.deg.csv")
write.csv(Ostreotide.de.results[[2]], "results/DEG/Ostreotide.EN-L5_6-NP.deg.csv")
write.csv(Ostreotide.de.results[[3]], "results/DEG/Ostreotide.EN-L6-CT.deg.csv")
write.csv(Ostreotide.de.results[[4]], "results/DEG/Ostreotide.EN-L6-IT.deg.csv")
write.csv(Ostreotide.de.results[[5]], "results/DEG/Ostreotide.EN-L6b.deg.csv")
write.csv(Ostreotide.de.results[[5]], "results/DEG/Ostreotide.EN-Newborn.deg.csv")
write.csv(Ostreotide.de.results[[9]], "results/DEG/Ostreotide.IPC-EN.deg.csv")

#save all deg results
saveRDS(L.054264.de.results, "results/DEG/L.054264.de.results.rds")
saveRDS(Ostreotide.de.results, "results/DEG/Ostreotide.de.results.rds")





####################################
#GSEA
setwd("SST/results/DEG/")
#import deg results
L.054264.EN.IT.Immature <- read.csv("L.054264.EN-IT-Immature.deg.csv")
L.054264.EN.L5_6.NP <- read.csv("L.054264.EN-L5_6-NP.deg.csv")
L.054264.EN.L6.CT <- read.csv("L.054264.EN-L6-CT.deg.csv")
L.054264.EN.L6.IT <- read.csv("L.054264.EN-L6-IT.deg.csv")
L.054264.EN.L6b <- read.csv("L.054264.EN-L6b.deg.csv")
L.054264.EN.Newborn <- read.csv("L.054264.EN-Newborn.deg.csv")

Ostreotide.EN.IT.Immature <- read.csv("Ostreotide.EN-IT-Immature.deg.csv")
Ostreotide.EN.L5_6.NP <- read.csv("Ostreotide.EN-L5_6-NP.deg.csv")
Ostreotide.EN.L6.CT <- read.csv("Ostreotide.EN-L6-CT.deg.csv")
Ostreotide.EN.L6.IT <- read.csv("Ostreotide.EN-L6-IT.deg.csv")
Ostreotide.EN.L6b <- read.csv("Ostreotide.EN-L6b.deg.csv")
Ostreotide.EN.Newborn <- read.csv("Ostreotide.EN-Newborn.deg.csv")
#GSEA on GO
deg_names <- c("L.054264.EN.IT.Immature", "L.054264.EN.L5_6.NP", "L.054264.EN.L6.CT", "L.054264.EN.L6.IT", "L.054264.EN.L6b", "L.054264.EN.Newborn",
               "Ostreotide.EN.IT.Immature", "Ostreotide.EN.L5_6.NP", "Ostreotide.EN.L6.CT", "Ostreotide.EN.L6.IT", "Ostreotide.EN.L6b", "Ostreotide.EN.Newborn")
degs <- list(L.054264.EN.IT.Immature, L.054264.EN.L5_6.NP, L.054264.EN.L6.CT, L.054264.EN.L6.IT, L.054264.EN.L6b, L.054264.EN.Newborn,
             Ostreotide.EN.IT.Immature, Ostreotide.EN.L5_6.NP, Ostreotide.EN.L6.CT, Ostreotide.EN.L6.IT, Ostreotide.EN.L6b, Ostreotide.EN.Newborn)

gsea_go_res_list <- list()
for (i in seq_along(deg_names)) {
  deg_name <- deg_names[i]
  deg <- degs[[i]]
  table_filt <- deg[!is.na(deg[,"logFC"]),]
  gene_list <- deg[,'t']
  names(gene_list) <- deg[,'X']
  gene_list <- sort(gene_list, decreasing = T)
  res <- gseGO(gene_list,
               OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL",
               exponent = 1,
               minGSSize = 10,
               maxGSSize = 500,
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               verbose = TRUE,
               seed = 0)
  gsea_go_res_list[[i]] <- res
  names(gsea_go_res_list)[i] <- deg_name
}

GSEA_go_res<- merge_result(gsea_go_res_list)
write.csv(GSEA_go_res@compareClusterResult, "GSEA_go_res.csv")

GSEA_table <- GSEA_go_res@compareClusterResult



#################################
#import chosen GO terms for plotting
# Load the chosen_GO_terms file
chosen_GO_terms <- read.csv("chosen_GO_terms.csv", header = 1)
# Load the GSEA_table file
GSEA_table_to_plot <- GSEA_table[GSEA_table$Description %in% chosen_GO_terms$Description,]
GSEA_table_to_plot$Significance <- ifelse(GSEA_table_to_plot$p.adjust < 0.05, "yes", "no")
GSEA_table_to_plot$Direction <- ifelse(GSEA_table_to_plot$NES < 0, "down", "up")
GSEA_table_to_plot$Description <- factor(GSEA_table_to_plot$Description, levels = rev(chosen_GO_terms$Description))
#dotplot of representative GO terms
p <- ggplot(data = GSEA_table_to_plot, aes(x = Cluster, y = Description, size = abs(NES))) +
  geom_point(shape = 21 ,aes(fill = Direction, color = Significance), stroke = 1.5) +
  scale_color_manual(breaks=c("yes","no"),
                     values=c("red", rgb(0, 0, 0, alpha=0))) +
  scale_fill_manual(values = c("#94d769", "#b98dba")) +
  scale_size_continuous(range = c(-1, 5))
p <- p + xlab("") + ylab("") + ggtitle("") +
  theme_dose(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(
    color = guide_legend(override.aes = list(size = 4)),  # Increase color legend key size
    fill = guide_legend(override.aes = list(size = 4)),   # Increase fill legend key size
    #size = guide_legend()
  )
p <- p + 
  scale_x_discrete(limits = c("L.054264.EN.Newborn", "Ostreotide.EN.Newborn", " ", "L.054264.EN.IT.Immature", "Ostreotide.EN.IT.Immature", " ",
                              "L.054264.EN.L6.IT", "Ostreotide.EN.L6.IT", " ", "L.054264.EN.L5_6.NP", "Ostreotide.EN.L5_6.NP", " ",
                              "L.054264.EN.L6.CT", "Ostreotide.EN.L6.CT", " ", "L.054264.EN.L6b", "Ostreotide.EN.L6b")) +
  scale_y_discrete(limits = rev(c("synapse organization", "cell junction assembly", "regulation of axonogenesis", " ", 
                                  "cell morphogenesis involved in neuron differentiation", "axon guidance", 
                                  "regulation of synapse assembly", " ", "chloride transmembrane transport", 
                                  "synapse assembly", "gamma-aminobutyric acid signaling pathway", " ", 
                                  "regulation of membrane potential", "potassium ion transport", 
                                  "modulation of chemical synaptic transmission", " ", "regulation of synaptic plasticity", 
                                  "synaptic signaling", "chemical synaptic transmission", " ", "regulation of neuron projection development", 
                                  "synaptic vesicle cycle", "developmental cell growth", " ", "protein acylation", 
                                  "forebrain development", "nucleoside triphosphate metabolic process", " ", " ", " ", " ",
                                  "alcohol metabolic process", "sterol metabolic process", "protein targeting to mitochondrion", " ", 
                                  "positive regulation of protein polyubiquitination", "reactive oxygen species metabolic process", 
                                  "regulation of ncRNA transcription", " ", "ATP synthesis coupled electron transport", 
                                  "cellular respiration", "hydrogen peroxide metabolic process", " ", 
                                  "cellular lipid catabolic process", "response to glucocorticoid", 
                                  "cellular response to corticosteroid stimulus", " ", "cholesterol efflux", 
                                  "regulation of PERK-mediated unfolded protein response", "tRNA metabolic process", " ", 
                                  "response to fatty acid", "cellular response to fatty acid", 
                                  "mitochondrial gene expression", " ", "basic amino acid transport", 
                                  "collagen fibril organization", "amino acid transmembrane transport"
  )))
p

ggsave(filename = "../gseGO_dotplot.pdf", width = 10, height = 11, units = "in")