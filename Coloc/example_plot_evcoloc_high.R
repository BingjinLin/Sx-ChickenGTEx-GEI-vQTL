library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
tissue=ARGS[1]
pheno1=ARGS[2]
chr=ARGS[3]

Otdir="/juice_data/zhh_independent/ZHH_update/vQTL/res/Plot/colocation/example/eQTL"

genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tissue)
sample=read.table(paste0(genodir,"/",tissue,".autosomal.geno.fam"))
sample_size <- as.numeric(nrow(sample))

vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore/")
file1 = paste0(vDir, tissue, ".vqtl_chr",chr, ".all.query.gz")
nominal1 <- fread(file1, header = T)
nominal1.1 <- nominal1[nominal1 $ Probe == pheno1, ]
nominal1.1$dist <- abs(nominal1.1$BP-nominal1.1$Probe_bp)
nominal1.2 <- nominal1.1[nominal1.1$dist < 1e6, ]
Otfile2=paste0(Otdir, tissue, "_",pheno1,"_vqtl_plot.txt")
write.table(nominal1.2, Otfile2, row.names=F, quote=F, sep="\t")

eDir<-paste0("/ossfs/ZHH_update/ZHH_update/eQTLmapping/",tissue,"/cis_mlm/")
file2 = paste0(eDir, tissue, ".cis_qtl_pairs.",chr, ".txt.gz")
nominal2 = fread(file2, header = T)
nominal2.1 <- nominal2[nominal2 $ pheno_id == pheno1, ]

merged_data <- merge(nominal1.2, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)

if(nrow(merged_data) > 0){
  
  # ignore pheno1s with missing data
    df_test <- as.data.frame(merged_data)

  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta_g1, pvalues=df_test$pval_g1,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$af <= 0.5, df_test$af, 1-df_test$af)))

  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  results=as.data.frame(my.res$results)
  hc_snp <- results$snp[results$SNP.PP.H4==max(results$SNP.PP.H4)][1]
}
###eqtl pre###
qtl_df <- data.frame(pheno_id = pheno1, variant_id=hc_snp)
Otfile1=paste0(Otdir, tissue, "_",pheno1,"_eqtl.extract.txt")
write.table(qtl_df, Otfile1, row.names=F, quote=F, sep="\t")
prefix <- paste0(tissue, "_chr", chr, "_", pheno1,".eqtl")
script="/soft/zhh/OmiGA/omiga_eQTL_plot.sh"
omiga_cmd <- paste("bash",script, tissue, Otfile1, Otdir, prefix)
system(omiga_cmd)

library(ggplot2)
library(patchwork)
###eqtl plot
snpr_df = fread(paste0(Otdir,prefix ,".snpr.plot.txt.gz"))
df1 = data.frame(pos = snpr_df$pos, logp = -log10(nominal2.1$pval_g1), ld = snpr_df$snp_r^ 2)

p1 <- ggplot(df1)+
  geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T) +
  geom_point(data = subset(df1, ld == 1), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df1$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="eQTL",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df1, ld == 1)$pos / 1000000, y = subset(df1, ld == 1)$logp, label = hc_snp, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

###vqtl plot
df2 = data.frame(pos = nominal1.2$BP, logp = -log10(nominal1.2$p), ld = snpr_df$snp_r[match(nominal1.2$BP, snpr_df$pos)]^ 2)

p2 <- ggplot(df2)+
  geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T) +
  geom_point(data = subset(df2, ld == 1), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df2$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="vQTL",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df2, ld == 1)$pos / 1000000, y = subset(df2, ld == 1)$logp, label = hc_snp, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

save(hc_snp, df1,df2, file=paste0(Otdir, tissue, "_",pheno1, "_mahattan.plot.Rdata"))
Otplot <- paste0(Otdir, tissue, "_",pheno1, "_plot.pdf")
ggsave(Otplot, p2/p1,width=5,height=5)

box_df = fread(paste0(Otdir,prefix ,".pheno_geno.plot.txt.gz"))
p3 <- ggplot(box_df,aes(x=factor(gt),y=pheno_adj))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  geom_point(alpha = 0.7, size = 1, width = 0.1)+
  xlab(box_df$variant_id[1])+
  ylab(paste0("Normalized expression ","(", pheno1,")"))+
  theme_classic()
Otplot <- paste0(Otdir, tissue, "_",pheno1, "_boxplot.pdf")
ggsave(Otplot, p3,width=3,height=3)
