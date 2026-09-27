library(data.table)
library(dplyr)
library(purrr)
library(coloc)
#tis="Tibial.growth.plate"
ARGS <- commandArgs(trailingOnly = TRUE)
tis = ARGS[1]
time = ARGS[2]

genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis)
sample1=read.table(paste0(genodir,"/",tis,".autosomal.geno.fam"))
sample_size1 <- as.numeric(dim(sample1)[1])
sub_geno <- paste0(genodir,"/",tis,"_",time,".autosomal.geno.fam")
if(file.exists(sub_geno)){
sample2=read.table(sub_geno)
sample_size2 <- as.numeric(dim(sample2)[1])

vdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/")
eDir<-file.path("/ossfs/GxE/Etime",tis,time)
lefile<-file.path("/juice_data/zhh_independent/ZHH_update/GxE/Etime/Overlap_vGene",paste0(tis,"_",time,".overlap_vGene.list"))
if(file.exists(lefile)){
gene_df <- read.table(lefile,head=T)
inte_gene<-gene_df$Gene
Otdir="/data/coloc/Etime/"
outfile = paste0(Otdir, tis,"_", time, ".vQTL_eQTL.pph4")
if(file.exists(outfile)){
Otdata<-fread(outfile)
Otdata<-unique(Otdata)
write.table(Otdata, outfile, row.names=F, sep="\t",quote=F)
Otgene<-unique(Otdata[,7])
gene_df <- gene_df[!gene_df$Gene %in% Otgene,]
}
for(chr in unique(gene_df$Chr)){
print(paste0("Processing on ",chr))
file1 = paste0(vdir, tis, ".vqtl_chr",chr, ".all.query.gz")
file2 = paste0(eDir, "/", tis,"_", time, ".cis_qtl_pairs.",chr, ".txt.gz")
if(file.exists(file1) & file.exists(file2)){
nominal1 = fread(file1, header = T)
nominal2 = fread(file2, header = T)
Genes<-unique(gene_df$Gene[gene_df$Chr==chr])
for(g in Genes){
nominal1.1 <- nominal1[nominal1$Gene==g,]
nominal2.1 <- nominal2[nominal2$pheno_id==g,]
merged_data <- merge(nominal1.1, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)

if(nrow(merged_data)>0){
  
  # ignore genes with missing data
  df_test <- as.data.frame(merged_data)
  if (any(is.na(df_test[c("b", "p", "beta_g1", "pval_g1")]))) { next }
  
  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size1,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta_g1, pvalues=df_test$pval_g1,
                   snp=df_test$SNP, type="quant", N=sample_size2,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

  print(paste0("Processing on chr",chr, ".", g, " (total is ",length(Genes), ")"))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  #snp.res <- as.data.frame(my.res$result)
  #snp.res$snp[snp.res$SNP.PP.H4==max(snp.res$SNP.PP.H4)]

  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$chrom = chr
  results$tissue = tis
  results$time = time

  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  write_header <- !file.exists(outfile)

  fwrite(results, outfile, append=!write_header, col.names=write_header, sep="\t")
}else{
print(paste0("There is no overlap variants between ",tis))
}
}
}
}
}
}
print("Colocalization analysis has been done")
