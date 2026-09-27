library(data.table)
library(dplyr)
library(purrr)
library(coloc)
#tissue="Tibial.growth.plate"
ARGS <- commandArgs(trailingOnly = TRUE)
tissue = ARGS[1]

vdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore/")
files<-list.files(vdir,pattern="\\.query$",full.name=T)
all_res<-data.frame()
for(f in files){
data<-fread(f)
all_res<-rbind(all_res,data)
}
vGene<-unique(all_res$Gene)
LieDir="/juice_data/zhh_independent/ZHH_update/ieQTLmapping/ieGene_qval005"
iefile<-list.files(LieDir,full.name=T)
pattern <- paste(tissue, collapse = "|")
subset_file <- iefile[grep(pattern, iefile)]
for(f in subset_file){
print(paste0("Processing on ",f))
ieData1 <- fread(f)
inte_gene <- intersect(ieData1$pheno_id, vGene)
tis<-sub("_.*", "",basename(f))
cell<-sub("^[^_]+_(.*)\\.cis_qtl_qval005\\.txt$", "\\1",basename(f))
name=paste0(tis,"_",cell)

vdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/")
ieDir<-paste0("/ossfs/ieQTL_mapping/Update/",tis,"/",cell,"/")
genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis)
sample=read.table(paste0(genodir,"/",tis,".autosomal.geno.fam"))
sample_size <- as.numeric(dim(sample)[1])
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/coloc/")
#dir.create(Otdir)
outfile = paste0(Otdir, name, ".vQTL_ieQTL.pph4")
if(file.exists(outfile)){
Otdata<-fread(outfile)
Otdata<-unique(Otdata)
write.table(Otdata, outfile, row.names=F, sep="\t",quote=F)
Otgene<-unique(Otdata[,7])
gene_df<-data.frame(gene=inte_gene, chrom=all_res$Chr[match(inte_gene, all_res$Gene)])
gene_df <- gene_df[gene_df$gene %in% Otgene,]
}

for(chr in unique(gene_df$chrom)){
print(paste0("Processing on ",chr))
file1 = paste0(vdir, tis, ".vqtl_chr",chr, ".all.query.gz")
file2 = paste0(ieDir, name, ".cis_qtl_pairs.",chr, ".txt.gz")
if(file.exists(file1) & file.exists(file2)){
nominal1 = fread(file1, header = T)
nominal2 = fread(file2, header = T)
Genes<-unique(gene_df$gene[gene_df$chrom==chr])
for(g in Genes){
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
  snp.res <- as.data.frame(my.res$result)
  snp.res$snp[snp.res$SNP.PP.H4==max(snp.res$SNP.PP.H4)]

  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$chrom = chr
  results$tissue = tis
  results$cell = cell
  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  fwrite(results, outfile, append=TRUE, col.names=FALSE, sep="\t")
}else{
print(paste0("There is no overlap variants between ",tis," and ",name))
}
}
}
}
}