library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
tissue=ARGS[1]
pheno1=ARGS[2]
chr=ARGS[3]
trait=ARGS[4]
loci=ARGS[5] #trait,loci,level

Otdir="/data/vQTL_coloc_Plot/example/specific_eQTL/"

genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tissue)
sample=read.table(paste0(genodir,"/",tissue,".autosomal.geno.fam"))
sample_size <- as.numeric(nrow(sample))
GWAS_geno_dir<-paste0("/ossfs/GWAS/",trait,"/",trait,"/")
GWAS_data <- read.table(paste0(GWAS_geno_dir,"GWAS1_", trait, ".fam"))
GWAS_size <- as.numeric(dim(GWAS_data)[1])

eDir<-paste0("/ossfs/ZHH_update/ZHH_update/eQTLmapping/",tissue,"/cis_mlm/")
file1 = paste0(eDir, tissue, ".cis_qtl_pairs.",chr, ".txt.gz")
Gdir= paste0("/ossfs/GWAS/",trait, "/colocation/locus/1e-5.inde/")
file2 <- paste0(Gdir, loci,".loci.txt") # nominal file of GWAS

nominal1 <- fread(file1, header = T)
nominal1.1 <- nominal1[nominal1 $ pheno_id == pheno1, ]

GWAS = fread(file2 , header = T)
merged_data <- merge(nominal1.1, GWAS, by.x="variant_id", by.y="rs", all=FALSE, suffixes=c("_eqtl","_gwas"))
  
  # ignore pheno1s with missing data
  df_test <- as.data.frame(merged_data)

  df1_coloc = list(beta=df_test$beta_g1, pvalues=df_test$pval_g1,
                   snp=df_test$variant_id, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$af_eqtl<= 0.5, df_test$af_eqtl, 1-df_test$af_eqtl)))
  df2_coloc = list(beta=df_test$beta, pvalues=df_test$p_wald,
                   snp=df_test$variant_id, type="quant", N=GWAS_size,
                   MAF=as.numeric(ifelse(df_test$af_gwas <= 0.5, df_test$af_gwas, 1-df_test$af_gwas)))

  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  results=as.data.frame(my.res$results)
  hc_snp <- results$snp[results$SNP.PP.H4==max(results$SNP.PP.H4)][1]

###eqtl plot###
qtl_df <- data.frame(pheno_id = pheno1, variant_id=hc_snp)
Otfile1=paste0(Otdir, tissue, "_",pheno1,"_eqtl.extract.txt")
write.table(qtl_df, Otfile1, row.names=F, quote=F, sep="\t")
prefix <- paste0(tissue, "_chr", chr, "_", pheno1,".eqtl")
script="/soft/zhh/OmiGA/omiga_eQTL_plot.sh"
omiga_cmd <- paste("bash",script, tissue, Otfile1, Otdir, prefix)
system(omiga_cmd)

###gwas plot###
Otfile3=paste0(Otdir, loci, "_gwas.txt")
Gfile2 <- paste0("/ossfs/GWAS/", trait, "/", trait, "/output/GWAS1_", trait, ".GWAS.assoc.txt")
G_df <- fread(Gfile2) %>% as.data.frame()
GWAS_SNP<- intersect(G_df$rs, nominal1.1$variant_id)
sig_file=paste0(Otdir,loci,".GWAS_sig-snp.list")
write.table(GWAS_SNP,sig_file,quote=F,row.names=F,col.names=F)
outfile0=paste0(Otdir,trait,",",hc_snp,"_genotype")
outfile1=paste0(Otdir,loci,".sig")
outfile2=paste0(Otdir,loci)
GWAS_geno<-paste0(GWAS_geno_dir,"GWAS1_", trait)
plink0<-paste0("plink --bfile ", GWAS_geno, " --chr-set 39 --keep-allele-order --recode A --snp ", hc_snp, " --out ", outfile0)
system(plink0)
plink1<-paste0("plink --bfile ", GWAS_geno, " --chr-set 39 --keep-allele-order --make-bed --extract ", sig_file, " --out ", outfile1)
system(plink1)
plink2<-paste0("plink --bfile ", outfile1, " --chr-set 39 --r2 --ld-snp ",hc_snp, " --ld-window-r2 0 --ld-window 2000000 --out ", outfile2)
system(plink2)
ld_score<-fread(paste0(outfile2,".ld"))
G_df_1 <- G_df[G_df$rs %in% GWAS_SNP, ]
G_df_1$LD_r2<-ld_score$R2[match(G_df_1$rs,ld_score$SNP_B)]
G_df_1 <- na.omit(G_df_1)
write.table(G_df_1,Otfile3,quote=F,row.names=F,sep="\t")

library(ggplot2)
library(patchwork)
###eqtl plot
snpr_df = fread(paste0(Otdir,prefix ,".snpr.plot.txt.gz"))
df1 = data.frame(pos = snpr_df$pos, logp = -log10(nominal1.1$pval_g1), ld = snpr_df$snp_r^ 2)

###vqtl plot
vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore/")
vfile = paste0(vDir, tissue, ".vqtl_chr",chr, ".all.query.gz")
pair_df1 <- fread(vfile, header = T)
pair_df1.1 <- pair_df1[pair_df1 $ Probe == pheno1, ]
pair_df1.1$dist <- abs(pair_df1.1$BP-pair_df1.1$Probe_bp)
pair_df1.2 <- pair_df1.1[pair_df1.1$dist < 1e6, ]
df2 = data.frame(pos = pair_df1.2$BP, logp = -log10(pair_df1.2$p), ld = snpr_df$snp_r[match(pair_df1.2$BP, snpr_df$pos)]^ 2)

###gwas plot
G_df_1$logp = -log10(G_df_1$p_wald)

save(hc_snp, df1, df2, G_df_1, file=paste0(Otdir, tissue, "_",pheno1, "_", loci, "_plot.Rdata"))

###boxplot
plot_dt = fread(paste0(Otdir,prefix, ".pheno_geno.plot.txt.gz"))
plot1 = ggplot(plot_dt,aes(x=gt,y=pheno_adj))+
  geom_point()+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  xlab(plot_dt$variant_id[1])+
  ylab(paste0("Normalized expression (",pheno1,")"))+
  theme_classic()
  
snp_geno <- fread(paste0(outfile0,".raw"))
names(snp_geno)[7] <- "Genotype"
allele0 <- G_df_1$allele0[G_df_1$rs==hc_snp]
allele1 <- G_df_1$allele1[G_df_1$rs==hc_snp]
snp_geno$Genotype[snp_geno$Genotype==0] <- paste0(allele0,"|",allele0)
snp_geno$Genotype[snp_geno$Genotype==1] <- paste0(allele0,"|",allele1)
snp_geno$Genotype[snp_geno$Genotype==2] <- paste0(allele1,"|",allele1)
snp_geno$Genotype <- factor(snp_geno$Genotype)
plot2 <- ggplot(snp_geno,aes(x=Genotype,y=PHENOTYPE))+
  geom_point()+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  xlab(plot_dt$variant_id[1])+
  ylab(trait)+
  theme_classic()

Otbox <- paste0(Otdir, loci,"_",tissue, ",",pheno1, "_box.pdf")
ggsave(Otbox, plot1/plot2,width=2.5,height=5)