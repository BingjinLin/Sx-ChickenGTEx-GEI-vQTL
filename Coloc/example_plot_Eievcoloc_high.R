library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
tissue=ARGS[1]
pheno1=ARGS[2]
chr=ARGS[3]
time=ARGS[4]

Otdir="/juice_data/zhh_independent/ZHH_update/vQTL/res/Plot/colocation/example/Etime/"

genodir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tissue)
sample1=read.table(paste0(genodir,"/",tissue,".autosomal.geno.fam"))
sample_size1 <- as.numeric(nrow(sample1))
sample2=read.table(paste0(genodir,"/",tissue,"_",time,".autosomal.geno.fam"))
sample_size2 <- as.numeric(nrow(sample1))

vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore/")
file1 = paste0(vDir, tissue, ".vqtl_chr",chr, ".all.query.gz")
nominal1 <- fread(file1, header = T)
nominal1.1 <- nominal1[nominal1 $ Probe == pheno1, ]
nominal1.1$dist <- abs(nominal1.1$BP-nominal1.1$Probe_bp)
nominal1.2 <- nominal1.1[nominal1.1$dist < 1e6, ]
Otfile2=paste0(Otdir, tissue, "_",pheno1,"_vqtl_plot.txt")
write.table(nominal1.2, Otfile2, row.names=F, quote=F, sep="\t")

eDir<-paste0("/ossfs/GxE/Etime/",tissue,"/",time,"/")
file2 = paste0(eDir, tissue,"_",time, ".cis_qtl_pairs.",chr, ".txt.gz")
nominal2 = fread(file2, header = T)
nominal2.1 <- nominal2[nominal2 $ pheno_id == pheno1, ]

merged_data <- merge(nominal1.2, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE)

if(nrow(merged_data) > 0){
  
  # ignore pheno1s with missing data
    df_test <- as.data.frame(merged_data)

  df1_coloc = list(beta=df_test$b, pvalues=df_test$p,
                   snp=df_test$SNP, type="quant", N=sample_size1,
                   MAF=as.numeric(ifelse(df_test$Freq <= 0.5, df_test$Freq, 1-df_test$Freq)))
  df2_coloc = list(beta=df_test$beta_g1, pvalues=df_test$pval_g1,
                   snp=df_test$SNP, type="quant", N=sample_size2,
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
efile=file.path("/ossfs/GxE/Etime",tissue,time,paste0(tissue,"_",time,".cis_qtl_pairs.",chr,".txt.gz"))
pairs_df = fread(efile)
pairs_df1 <- pairs_df[pairs_df$pheno_id==pheno1,]
bim=fread(paste0(genodir,"/",tissue,"_",time,".autosomal.geno.bim"))
pairs_df1$pos <- bim$V4[match(pairs_df1$variant_id,bim$V2)]
snpr_df = fread(paste0(Otdir,prefix ,".snpr.plot.txt.gz"))
df1 = data.frame(pos = pairs_df1$pos, logp = -log10(pairs_df1$pval_g1), ld = (snpr_df$snp_r[match(pairs_df1$pos,snpr_df$pos)])^ 2)

p1 <- ggplot(df1)+
  geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T) +
  geom_point(data = subset(df1, ld == 1), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df1$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title=paste0(time," eQTL"),x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df1, ld == 1)$pos / 1000000, y = subset(df1, ld == 1)$logp, label = hc_snp, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

pattern=paste0("\\.cis_qtl_pairs\\.",chr,"\\.txt.gz$")
other_times <- list.files(paste0("/ossfs/GxE/Etime/",tissue),pattern=pattern, full.name=T,recursive=T)
other_times <- setdiff(other_times, efile)
for(i in seq(1,length(other_times))){
pairs_of = fread(other_times[i])
pairs_of1 <- pairs_of[pairs_of$pheno_id==pheno1,]
pairs_of1$pos <- bim$V4[match(pairs_of1$variant_id,bim$V2)]
pdf = data.frame(pos = pairs_of1$pos, logp = -log10(pairs_of1$pval_g1), ld = (snpr_df$snp_r[match(pairs_of1$pos,snpr_df$pos)])^ 2)

otime <- basename(dirname(other_times[i]))
pic <- ggplot(pdf)+
  geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T) +
  geom_point(data = subset(pdf, ld == 1), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df1$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title=paste0(otime," eQTL"),x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df1, ld == 1)$pos / 1000000, y = subset(df1, ld == 1)$logp, label = hc_snp, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))
p1 <- p1/pic
}
###vqtl plot
vfile=file.path("/ossfs/ZHH_update/ZHH_update/eQTLmapping",tissue,"cis_mlm",paste0(tissue,".cis_qtl_pairs.",chr,".txt.gz"))
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

Otplot <- paste0(Otdir, tissue, "_",pheno1, "_plot.pdf")
ggsave(Otplot, p2/p1+plot_layout(heights = c(1, 3)),width=5,height=7)

pheno_df = fread(paste0(Otdir,prefix ,".pheno_geno.plot.txt.gz"))
pheno_df$IID <- sample1$V2
otime_list <- basename(dirname(other_times))
Etime_sam <- data.frame(time=time, sample=sample2$V2)
for(otime in otime_list){
sample=read.table(paste0(genodir,"/",tissue,"_",otime,".autosomal.geno.fam"))
tmp_sam <- data.frame(time=otime, sample=sample$V2)
Etime_sam <- rbind(Etime_sam, tmp_sam)
}
pheno_df$time <- Etime_sam$time[match(pheno_df$IID,Etime_sam$sample)]
pic <- ggplot(pheno_df,aes(x=gt,y=pheno_adj))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  facet_wrap(facets = vars(time))+
  xlab(pheno_df$variant_id[1])+
  ylab(paste0("Normalized expression (",pheno1,")"))+
  theme_classic()

Otbox <- paste0(Otdir, tissue, "_",pheno1, "_boxplot.pdf")
ggsave(Otbox, pic,width=6,height=2.5)
