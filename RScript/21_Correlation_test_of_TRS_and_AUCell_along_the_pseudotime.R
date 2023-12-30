library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
library(tidyverse)
library(SummarizedExperiment)
library(Matrix)
library(BuenColors)
library(cowplot)
library(diffloop)
library(igraph)
library(stringr)
traitname <- "Trait"
wdir <- "/path/to/traits/TRS/"
setwd(wdir)
out <- readRDS("path/to/EN_trajectory.rds")
TRS <- fread("/path/to/traits/TRS")
TRS2 <- TRS[match(TRS$cell_name, rownames(out)),]

## Scatterplots and Boxplots

p1 <- ggplot(data=TRS2, aes(UMAP1, UMAP2, color=pseudotime)) + geom_point(size=0.8, na.rm = TRUE) +
	pretty_plot() + scale_color_gradientn(colors = jdb_palette("blue_cyan")) + scale_alpha() + geom_path(data=data.frame(out[[2]]), aes(x, y, color=NULL), 
        arrow = arrow(type = "open", angle = 20, length = unit(0.1, "inches")), size = 3)

p2 <- ggplot(data=TRS2,  aes(x=pseudotime, y=SCZ2_TRS, color=type))  +
  ggrastr::geom_point_rast(size=0.8, na.rm = TRUE) + geom_smooth(colour = "black") + 
  scale_color_manual(values=lev[levels(TRS2$type),'color']) +
  pretty_plot() + scale_color_gradientn(colors=jdb_palette("brewer_marine")) +
  ylab("Trait relevant score (TRS)") + xlab("Pseudo time") + ggtitle(traitname) +
  theme(axis.title.y=element_blank())
pdf(paste0(traitname,".pseudotime.pdf"),wi=15,he=10)
p1
p2
dev.off()

## correlation test with AUC data of all eRegulons
auc_data <- readRDS("/path/to/auc_data.RDS")
mycores=1
res_cor.test_list <- parallel::mclapply(1:nrow(auc_data), mc.cores = mycores, function(i){
  x <- cor.test(auc_data[i, ] %>% unlist, TRS2[,traitname])
  res_cor.test <- c(x$estimate, x$p.value)
  if (i %% 100 == 0) {print(i)}
  return(res_cor.test)
})
res_cor.test_df <- as.data.frame(res_cor.test_list) %>% t %>% as.data.frame
rownames(res_cor.test_df) <- rownames(auc_data)
colnames(res_cor.test_df) <- c("rho", "p.value")
res_cor.test_df$fdr <- p.adjust(res_cor.test_df[, 2], method = "bonferroni")
res_cor.test_df$log10fdr <- -log10(res_cor.test_df$fdr)
res_cor.test_df$posorneg <- res_cor.test_df[, 1]>0
res_cor.test_df <- res_cor.test_df[order(0-res_cor.test_df$rho), ]
res_cor.test_df$order <- 1:nrow(res_cor.test_df)
rownames(res_cor.test_df)=gsub(pattern = "\\,",replace="_",x=rownames(res_cor.test_df))
res_cor.test_df$category="active"
res_cor.test_df[grepl(pattern = "\\-_",rownames(res_cor.test_df)),'category']='repressive'

pos_tf <- rownames(res_cor.test_df[res_cor.test_df=='active',])[1:10]
neg_tf <- rownames(res_cor.test_df[res_cor.test_df=='repressive',])[1:10]
pos_tf_idx <- res_cor.test_df[pos_tf,'order'] %>% as.numeric
neg_tf_idx <- res_cor.test_df[neg_tf,'order'] %>% as.numeric
pos_tf_name <- str_split(pos_tf, "_", simplify=T)[, 1]
neg_tf_name <- str_split(neg_tf, "_", simplify=T)[, 1]
tf <- c(pos_tf, neg_tf)
tf_name <- c(pos_tf_name, neg_tf_name)

manual_selected_top <- res_cor.test_df[tf, ]
manual_selected_top <- data.frame(manual_selected_top, tf_name=tf_name)
my_col <- c("#007AB7", "#007AB7", "#208DC3", "#9ACDE7", "#E3E3E3", "#E3E3E3", "#FCEDAA", "#F8D32F", "#EDAD0B", "#EDAD0B")

p1 <- ggplot(res_cor.test_df, aes(x = order, y = rho, fill = rho)) + 
geom_col(position = "identity") + pretty_plot() + scale_fill_gradientn(colors = my_col)+
ggrepel::geom_label_repel(data=manual_selected_top, aes(order, rho, label = tf_name), color = 'black', fill="white",
size = 1.5, max.overlaps = Inf,point.padding = 0, # additional padding around each point
min.segment.length = 0, # draw all line segments
max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
box.padding = 0.3) # additional padding around each text label
+ ylab("Spearman correlation")+ pretty_plot(fontsize = 10)

tfs <- c(pos_tf,neg_tf)
tfs <- tfs[grepl(pattern = "\\+_+",x = tfs)]
rownames(auc_data)=gsub(pattern = "\\,",replace="_",x=rownames(auc_data))
colnames(auc_data)=colnames(auc)
auc_part_df <- as.data.frame(t(auc_data[tfs,rownames(TRS2)]))
colnames(auc_part_df) <- gsub(pattern = "\\_.*",replacement ="",x=colnames(auc_part_df))
auc_part_std <- as.data.frame(apply(auc_part_df,2,scale))[,tfs]
rownames(auc_part_std) <-rownames(auc_part_df)
auc_part_std$pseudotime <- TRS2$pseudotime
auc_part_long <- gather(auc_part_std,TF,auc,MEF2C：PAX6,factor_key = T) ## modify for each trait

p2 <- ggplot(data = auc_part_long, aes(x = pseudotime,fill=auc,y=TF)) + ggrastr::geom_tile_rast() +scale_fill_gradientn(colors  =jdb_palette("brewer_jamaica"),limits=c(-2.5,2.5))+pretty_plot()

pdf(paste0(traitname,"_AUC_time.pdf"))
plot_grid(p1, p2,labels = "AUTO", ncol = 1)
dev.off()