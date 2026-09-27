library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
tis = ARGS[1]

genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis)
sample=read.table(paste0(genodir,"/",tis,".autosomal.geno.fam"))
sample_size <- as.numeric(dim(sample)[1])

vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/")
Otdir_1=paste0("/data/coloc/Growth_iecoloc/",tis,"/")
dir.create(Otdir_1)
task_df <- read.table(paste0("/juice_data/zhh_independent/ZHH_update/GxE/Growth/",tis,".Growth_vGene.txt"),head=T)
trait_list <- unique(task_df$trait)
for(tr in trait_list){
ieDir=file.path("/ossfs/GxE/veQTL_Growth", tis, tr)
chr_list <- unique(task_df$Chr[task_df$trait==tr])
for(chr in chr_list){
file1 = paste0(vDir, tis, ".vqtl_chr",chr, ".all.query.gz") # nominal file of vQTL
file2 = paste0(ieDir, "/", tis, "_", tr, ".cis_qtl_pairs.",chr, ".txt.gz") # nominal file of ieQTL
nominal1 <- fread(file1, header = T)
nominal2 <- fread(file2, header = T)
Genes<-unique(task_df$Gene[task_df$trait==tr & task_df$Chr==chr])
for(g in Genes){
nominal1.1 <- nominal1[nominal1$Gene==g,]
nominal2.1 <- nominal2[nominal2$pheno_id==g,]

merged_data <- merge(nominal1.1, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)

if(!is.null(merged_data) && nrow(merged_data) > 0){  
  # ignore genes with missing data
  df_test <- as.data.frame(merged_data)
  if (any(is.na(df_test[c("b", "p", "beta_g2", "pval_g2")]))) { next }
  
  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta_g2, pvalues=df_test$pval_g2,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

print(paste("Processing on", tis,"'s", g, "in trait", tr, "(total is", length(Genes), ")"))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$chrom = chr
  results$tissue = tis
  results$trait = tr

  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  outfile = paste0(Otdir_1, tis,"_", tr, ".vQTL.pph4")
  write_header <- !file.exists(outfile)
  fwrite(results, outfile, append=TRUE, col.names=write_header, sep="\t")
}else{
print(paste0("There is no overlap variants between ",tis," and ",basename(f)))
}
}
}
}
print("All have been done")
