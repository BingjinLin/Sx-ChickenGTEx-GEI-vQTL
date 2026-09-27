library(data.table)
library(dplyr)
library(purrr)
library(coloc)
#tiss="Tibial.growth.plate"
ARGS <- commandArgs(trailingOnly = TRUE)
tiss = ARGS[1]

vdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tiss,"/zscore/")
files<-list.files(vdir,pattern="\\.query$",full.name=T)
all_res<-data.frame()
for(f in files){
data<-fread(f)
all_res<-rbind(all_res,data)
}
vGene<-unique(all_res$Gene)

sDir<-paste0("/juice_data/zhh_independent/ZHH_update/Sex_Straitified/sex_biaed/tensorQTL/",tiss,"/")
sex_df<-fread("/juice_data/zhh_independent/ZHH_update/Sex_Straitified/sex_biaed_gene.list")
tis_df<-sex_df[sex_df$tissue==tiss, ]
inte_gene<-intersect(vGene, tis_df$pid)

genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tiss)
sample=read.table(paste0(genodir,"/",tiss,".autosomal.geno.fam"))
sample_size <- as.numeric(dim(sample)[1])
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tiss,"/zscore/coloc/")
dir.create(Otdir)
gene_df<-tis_df[tis_df$pid %in% inte_gene, ]
gene_df$chrom <- all_res$Chr[match(gene_df$pid, all_res$Gene)]

for(chr in unique(gene_df$chrom)){
print(paste0("Processing on ",chr))
file1 = paste0(vdir, tiss, ".vqtl_chr",chr, ".all.query.gz")
file2 = paste0(sDir, tiss, ".cis_qtl_pairs.",chr, ".txt.gz")
if(file.exists(file1) & file.exists(file2)){
nominal1 = fread(file1, header = T)
nominal2 = fread(file2, header = T)
Genes<-unique(gene_df$pid[gene_df$chrom==chr])
for(g in Genes){
nominal1.1 <- nominal1[nominal1$Gene==g,]
nominal2.1 <- nominal2[nominal2$phenotype_id==g,]
merged_data <- merge(nominal1.1, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)

if(nrow(merged_data)>0){
  
  # ignore genes with missing data
  df_test <- as.data.frame(merged_data)
  if (any(is.na(df_test[c("b", "p", "b_gi", "pval_gi")]))) { next }
  
  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$b_gi, pvalues=df_test$pval_gi,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

  print(paste0("Processing on chr",chr, ".", g, " (total is ",length(Genes), ")"))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  snp.res <- as.data.frame(my.res$result)
  snp.res$snp[snp.res$SNP.PP.H4==max(snp.res$SNP.PP.H4)]

  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$chrom = chr
  results$tissue = tiss
  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  outfile = paste0(Otdir, tiss, ".vQTL_sex-eQTL.pph4")
  fwrite(results, outfile, append=TRUE, col.names=FALSE, sep="\t")
}else{
print(paste0("There is no overlap variants on ",tiss))
}
}
}
}