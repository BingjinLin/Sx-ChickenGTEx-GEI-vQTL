library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
tis = ARGS[1]

vdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore_update/")
task_file <- file.path(vdir, paste0(tis,"_cieQTL_intersection.txt"))
if(file.exists(task_file)){

task_df <- read.table(task_file,head=T,sep="\t")
cell_list <- unique(task_df$cell)

for(cie in cell_list){
name <- paste(tis, cie,sep="_")
ieDir<-file.path("/ossfs/ieQTL_mapping/DWLS_Update",tis,name)
gene_list <- task_df$ieGene[task_df$cell == cie]

genodir=paste0("/ossfs/ZHH_update/ZHH_update/Genotype/",tis)
sample=read.table(paste0(genodir,"/",tis,".autosomal.geno.fam"))
sample_size <- as.numeric(dim(sample)[1])
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/coloc/")
#dir.create(Otdir)
outfile = paste0(Otdir, name, ".vQTL_cieQTL.pph4")

if(length(gene_list) >0){
chr_list <- unique(task_df$Chr[task_df$ieGene %in% gene_list & task_df$cell == cie])
for(chr in chr_list){
file1 = paste0(vdir, tis, ".vqtl_chr",chr, ".all.query.gz")
file2 = paste0(ieDir, "/", name, ".cis_qtl_pairs.",chr, ".txt.gz") # nominal file of cieQTL
if(file.exists(file1) & file.exists(file2)){
nominal1 = fread(file1, header = T)
nominal2 = fread(file2, header = T)
Genes <- task_df$ieGene[task_df$cell == cie & task_df$Chr==chr & task_df$ieGene %in% gene_list]
for(g in Genes){
print(paste0("Processing on ", name,"'s chr",chr, ".", g, " (total is ",length(Genes), ")"))
nominal1.1 <- nominal1[nominal1$Gene==g,]
nominal2.1 <- nominal2[nominal2$pheno_id==g,]
merged_data <- merge(nominal1.1, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)

if(nrow(merged_data)>0){
  
  # ignore genes with missing data
  df_test <- as.data.frame(merged_data)
  if (any(is.na(df_test[c("b", "p", "beta_g2", "pval_g2")]))) { next }
  
  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta_g2, pvalues=df_test$pval_g2,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

#print(paste0("Processing on ", name,"'s chr",chr, ".", g, " (total is ",length(Genes), ")"))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  #snp.res <- as.data.frame(my.res$result)
  #snp.res$snp[snp.res$SNP.PP.H4==max(snp.res$SNP.PP.H4)]

  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$chrom = chr
  results$tissue = tis
  results$cell = cie
  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  fwrite(results, outfile, append=file.exists(outfile), col.names=!file.exists(outfile), sep="\t")
}else{
print(paste0("There is no overlap variants between ",tis," and ",name))
}
}
}
}
}
}
}
print(paste("All of cieQTL in",tis,"have been colocalized finished"))
