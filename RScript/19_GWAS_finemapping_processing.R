
## 01 processing GWAS posterior from GCTA-COJO

library(data.table)
library(dplyr)
raw_gwas_data <-"/path/to/raw_gwas_summary"
disease <-"DISEASE"
wdir <-paste0("/path/to/output_dir/",disease)
setwd(wdir)
gcta_file <- paste0("/path/to/gcta_cojo/",disease,"gcta_input.ma")
full_summarystats <- data.table::fread(gcta_file) # columns may include: SNP, CHR, BP, NEA, EA, EA_FRQ, BETA, SE, P, N
topsnp <- data.table::fread("leadsnp")
rownames(full_summarystats)=full_summarystats$SNP
raw_data<-data.table::fread(raw_gwas_data,fill=T)
raw_data$SNP <- raw_data$SNPID_UKB
full_summarystats <- merge(full_summarystats,raw_data[,list(SNP,CHR,BP)],by.x="SNP",by.y="SNP",all.x=T,all.y=F)

result <- list()
for (i in 1:nrow(topsnp)) {
  window_size <- (1e6 / 2)
  window_name <- topsnp[i,]$V2
  leadsnp_name <- topsnp[i,]$V2
  leadsnp_pos <- topsnp[i,]$V3
  leadsnp_chr <- topsnp[i,]$V1
  pos_min <- leadsnp_pos - window_size
  pos_max <- leadsnp_pos + window_size
  
  print(paste0("working on: ", window_name))
  windowss <- full_summarystats %>% 
    distinct(SNP, .keep_all=TRUE) %>%
    drop_na() %>%
    dplyr::filter(CHR==leadsnp_chr & (between(BP, pos_min, pos_max)) ) %>%
    mutate(freq = ifelse(freq>0.5, 1-freq, freq)) %>%
    dplyr::filter(freq>0.01)
  
  coloc_ds_input <- list(
    type="binary",			## "quant" if traits is a quantitative trait
    pvalues=windowss$p,
    N=windowss$N,
    MAF=windowss$freq,
    beta=windowss$b,
    varbeta=(windowss$se^2),
    #s=NULL,
    #sdy=NULL,
    snp=windowss$SNP)
  
  priornum <- 1e-04
  
  result[[paste0(window_name)]] <- coloc::finemap.abf(coloc_ds_input, p1 = priornum) %>% dplyr::filter(as.character(prior) == as.character(priornum)) %>% arrange(-SNP.PP) %>% 
	mutate(SNP.PP.cumsum = cumsum(SNP.PP)) %>% mutate(incl95cs = ifelse( round(SNP.PP.cumsum, 2) <= 0.95, TRUE, FALSE )) %>% mutate(LDsentinel = ifelse (snp == leadsnp_name, TRUE, FALSE))
}
result.df <- data.table::rbindlist(result)
write.table(result.df,"abf_pp.tsv",sep="\t",row.names = F,quote = F)
colnames(result.df)[5]='SNP'
result.df <- left_join(result.df,raw_data,by='SNP')
result.df$CHR <- paste0('chr',result.df$CHR)
result.df$start <- result.df$BP-1
result.df <- result.df[,c("CHR","start","BP",'SNP','SNP.PP')]
write.table(result.df,'PP.abf.hg19.bed'),sep="\t",row.names = F,quote = F,col.names = F)


## 02 Liftover by UCSC toolkits
liftOver PP.abf.hg19.bed /path/to/ucsc/chains/hg19ToHg38.over.chain PP.abf.hg38.bed unMapped