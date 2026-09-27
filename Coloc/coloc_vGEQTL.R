library(data.table)
library(dplyr)
library(purrr)
library(coloc)
Gvfile <- "/juice_data/zhh_independent/ZHH_update/GxE/veQTL/GxE.vGene_res.txt"
all_Gv <- fread(Gvfile)
stat_v_GE <- data.frame(table(all_Gv$tissue, all_Gv$trait))
stat_v_GE <- stat_v_GE[stat_v_GE$Freq!=0,]

###colocation###
tissue_task <- unique(stat_v_GE$Var1)
for(tis in tissue_task){
print(paste0("Processing on tissue: ",tis))
vdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/")
GEdir=paste0("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/",tis)
file1 = paste0(vdir, tis, ".vqtl_chr",chr, ".all.query.gz")
    if (!file.exists(file1)) {
      print(paste0("vQTL file not found: ", file1))
      next
    }
    nominal1 <- fread(file1, header = T)
	
##sample size1
genodir1=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis)
sample1=read.table(paste0(genodir1,"/",tis,".autosomal.geno.fam"))
sample_size1 <- as.numeric(dim(sample1)[1])
##Output
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore/coloc/")
##sub data
task_df1 <- all_Gv[all_Gv$tissue==tis,]
chr_list <- unique(task_df1$chrom)
for(chr in chr_list){
print(paste0("Processing on ",chr))

trait_task <- unique(task_df1$trait[task_df1$chrom==chr])
for(tr in trait_task){
print(paste0("Processing on trait: ",tr))
file2 = file.path(GEdir,tr,paste0(tis,"_",tr,".cis_qtl_pairs.",chr,".txt.gz"))
if(file.exists(file1) & file.exists(file2)){
nominal2 = fread(file2, header = T)

##sample size1
sample2=read.table(paste0("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/tmp/",tis,"_",tr,".sample.list"))
sample_size2 <- as.numeric(dim(sample2)[1])

Genes<-unique(task_df1$pheno_id[task_df1$chrom==chr & task_df1$trait==tr])
if(length(Genes) == 0) next

for(g in Genes){
nominal1.1 <- nominal1[nominal1$Gene==g,]
nominal2.1 <- nominal2[nominal2$pheno_id==g,]
merged_data <- merge(nominal1.1, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)
if(nrow(merged_data)>0){
  
  # ignore genes with missing data
  df_test <- as.data.frame(merged_data)
  if (any(is.na(df_test[c("b", "p", "beta_g2", "pval_g2")]))) { next }
  
  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size1,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta_g2, pvalues=df_test$pval_g2,
                   snp=df_test$SNP, type="quant", N=sample_size2,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

print(paste0("Processing on ", tr,"'s chr",chr, ".", g, " (total is ",length(Genes), ")"))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  snp.res <- as.data.frame(my.res$result)
  snp.res$snp[snp.res$SNP.PP.H4==max(snp.res$SNP.PP.H4)]

  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$chrom = chr
  results$tissue = tis
  results$trait = tr
  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  outfile = paste0(Otdir, tis, ".vQTL_GEQTL.coloc")
          write_header <- FALSE
        if (!file.exists(outfile)) {
          write_header <- TRUE
        }
  fwrite(results, outfile, append=TRUE, col.names=write_header, sep="\t")
  write_header <- FALSE 
}else{
print(paste0("There is no overlap variants on ",tis," ",tr, "between GEQTL and vQTL"))
}
}
}
}
}
}
