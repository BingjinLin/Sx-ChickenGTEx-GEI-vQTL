library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
start <- as.integer(ARGS[1])
end   <- as.integer(ARGS[2])
trait = ARGS[3]
locidr = ARGS[4]

Gdir=file.path("/ossfs/GWAS",trait,trait)
Gsample=read.table(paste0(Gdir,"/","GWAS1_", trait,".fam"))
sample_size2 <- as.numeric(dim(Gsample)[1])

tissue<-readLines("/soft/zhh/splicing/tiss.txt")
#tis="Liver"
tissue_list=tissue[start:end]
print(paste("Task for", tissue_list))

loci<-fread(paste0(locidr, "GWAS1_", trait, ".GWAS.1e-5.jma.cojo"))
loci[, Chr := as.character(Chr)]
locidr_1=paste0(locidr,"locus/1e-5.inde")
file2<-list.files(locidr_1, pattern="\\.txt$", full.name=T) # nominal file of GWAS
loci_stat <- data.frame()
for(f in file2){
loci_df = fread(f, header = T)
tmp_stat <- data.frame(chrom=unique(loci_df$chr), filename=f)
loci_stat <- rbind(loci_stat, tmp_stat)
}
loci_stat$chrom <- as.character(loci_stat$chrom)

chrom<-unique(loci$Chr)
for(chr in chrom){
for(tis in tissue_list){
genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis)
sample=read.table(paste0(genodir,"/",tis,".autosomal.geno.fam"))
sample_size1 <- as.numeric(dim(sample)[1])

vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tis,"/zscore_update/")
Otdir_1=paste0(vDir,"coloc/trait_GWAS")
dir.create(Otdir_1, recursive = TRUE, showWarnings = FALSE)

file1 = paste0(vDir, tis, ".vqtl_chr",chr, ".all.query.gz") # nominal file of vQTL
  if (!file.exists(file1)) next
  nominal1 <- fread(file1, header = T)
  sig_vdf <- fread(paste0(vDir, tis, "_Vgene.lead.txt"))
  sig_nominal1 <- nominal1[nominal1$Gene %in% sig_vdf$Gene,]

chr_files <- loci_stat$filename[loci_stat$chrom==chr]
for(f in chr_files){

GWAS = fread(f, header = T)
GWAS[, chr := as.character(chr)]
if (!all(GWAS$chr == chr)) {
next
}
nominal1.1 <- sig_nominal1[sig_nominal1$SNP %in% GWAS$rs, ]
pheno_id<-unique(nominal1.1$Probe)
for(g in pheno_id){
nominal1.2 <- nominal1.1[nominal1.1 $ Probe == g, ]
merged_data <- merge(nominal1.2, GWAS, by.x="SNP", by.y="rs", all=FALSE, suffixes=c("_vqtl","_gwas"))

if(!is.null(merged_data) && nrow(merged_data) > 0){  
  # ignore genes with missing data
  df_test <- as.data.frame(merged_data)
  df_test <- df_test[
    complete.cases(df_test[, c("b", "p", "beta", "p_wald")]),
]
  if (nrow(df_test) < 10) { next }

  #if (any(is.na(df_test[c("b", "p", "beta", "p_wald")]))) { next }
  
  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size1,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta, pvalues=df_test$p_wald,
                   snp=df_test$SNP, type="quant", N=sample_size2,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

print(paste("Processing on", tis,"'s", g, "in file", f, "(total is", length(pheno_id), ")"))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  results=list(my.res$summary)
  as.data.frame(t(unlist(results))) -> results
  results$gene = g
  results$bin = sub("(.*)\\..*\\..*", "\\1",basename(f))
  results$chrom = chr
  results$tissue = tis
  results$trait = trait

  # columns: nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf pheno1 pheno2 chrom tissue
  outfile = paste0(Otdir_1, "/", trait,"_", tis, ".vQTL_gwas_1e-5.update.pph4")
  write_header <- !file.exists(outfile)
  fwrite(results, outfile, append=TRUE, col.names=write_header, sep="\t")
}else{
print(paste0("There is no overlap variants between ",tis," and ",basename(f)))
}
}
}
}
}

print("All have been done")

