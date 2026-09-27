library(data.table)
library(dplyr)
library(purrr)
library(coloc)
ARGS <- commandArgs(trailingOnly = TRUE)
tissue=ARGS[1]
pheno1=ARGS[2]
chr=ARGS[3]

Otdir="/juice_data/zhh_independent/ZHH_update/vQTL/res/Plot/colocation/example/eQTL/"

vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore/")
file1 = paste0(vDir, tissue, ".vqtl_chr",chr, ".all.query.gz")
nominal1 <- fread(file1, header = T)
nominal1.1 <- nominal1[nominal1 $ Probe == pheno1, ]
nominal1.1$dist <- abs(nominal1.1$BP-nominal1.1$Probe_bp)
nominal1.2 <- nominal1.1[nominal1.1$dist < 1e6, ]
Otfile1.1=paste0(Otdir, tissue, "_",pheno1,"_vqtl_plot.txt")
write.table(nominal1.2, Otfile1.1, row.names=F, quote=F, sep="\t")
lead_vqtl=nominal1.2$SNP[nominal1.2$p==min(nominal1.2$p)][1]
vqtl_pos=nominal1.2$BP[nominal1.2$SNP==lead_vqtl]
vqtl_df <- data.frame(pheno_id = pheno1, variant_id=lead_vqtl)
Otfile1.2=paste0(Otdir, tissue, "_",pheno1,"_vqtl.extract.txt")
write.table(vqtl_df, Otfile1.2, row.names=F, quote=F, sep="\t")
prefix0 <- paste0(tissue, "_chr", chr, "_", pheno1,".vqtl")
script="/soft/zhh/OmiGA/omiga_eQTL_plot.sh"
omiga_cmd0 <- paste("bash",script, tissue, Otfile1.2, Otdir, prefix0)
system(omiga_cmd0)

eqtl_df <- fread(paste0("/ossfs/ZHH_update/ZHH_update/leading_eQTL/FDR/",tiss,".cis_qtl.0.05.txt"))
lead_eqtl=eqtl_df[eqtl_df$pheno_id==pheno1,c("pheno_id","variant_id")]
Otfile1.3=paste0(Otdir, tissue, "_",pheno1,"_eqtl.extract.txt")
write.table(lead_eqtl, Otfile1.3, row.names=F, quote=F, sep="\t")
prefix1 <- paste0(tissue, "_chr", chr, "_", pheno1,".eqtl")
script="/soft/zhh/OmiGA/omiga_eQTL_plot.sh"
omiga_cmd1 <- paste("bash",script, tissue, Otfile1.3, Otdir, prefix1)
system(omiga_cmd1)

###eqtl pre###
efile=file.path("/ossfs/ZHH_update/ZHH_update/eQTLmapping",tissue,"cis_mlm",paste0(tissue,".cis_qtl_pairs.",chr,".txt.gz"))
pairs_df = fread(efile)
pairs_df1 <- pairs_df[pairs_df$pheno_id==pheno1,]
lead_eqtl <- pairs_df1$variant_id[pairs_df1$pval_g1==min(pairs_df1$pval_g1)]
qtl_df <- data.frame(pheno_id = pheno1, variant_id=lead_eqtl)
Otfile2=paste0(Otdir, tissue, "_",pheno1,"_eqtl.extract.txt")
write.table(qtl_df, Otfile2, row.names=F, quote=F, sep="\t")
prefix <- paste0(tissue, "_chr", chr, "_", pheno1,".eqtl")
script="/soft/zhh/OmiGA/omiga_eQTL_plot.sh"
omiga_cmd <- paste("bash",script, tissue, Otfile2, Otdir, prefix)
system(omiga_cmd)


library(ggplot2)
library(patchwork)
###eqtl plot
lable_df <- pairs_df1[pairs_df1$variant_id %in% c(lead_vqtl,lead_eqtl),]
snpr_df = fread(paste0(Otdir,prefix1 ,".snpr.plot.txt.gz"))
df1 = data.frame(pos = snpr_df$pos, logp = -log10(pairs_df1$pval_g1), ld = snpr_df$snp_r^ 2)
df1$label <- "NA"
df1$label[which(pairs_df1$variant_id %in% c(lead_vqtl,lead_eqtl))] <- "highlight"
if(eqtl_pos < vqtl_pos){
labels <- c(lead_eqtl,lead_vqtl)
}else{
labels <- c(lead_vqtl,lead_eqtl)
}
p1 <- ggplot(df1)+
  geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T) +
  geom_point(data = subset(df1, label == "highlight"), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df1$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="eQTL",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df1, label == "highlight")$pos / 1000000, y = subset(df1, label == "highlight")$logp, label = lable_df$variant_id, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

###vqtl plot
snpr_df = fread(paste0(Otdir,prefix0 ,".snpr.plot.txt.gz"))
df2 = data.frame(pos = nominal1.2$BP, logp = -log10(nominal1.2$p), ld = snpr_df$snp_r[match(nominal1.2$BP, snpr_df$pos)]^ 2)
df2$label <- "NA"
df2$label[which(nominal1.2$SNP %in% c(lead_vqtl,lead_eqtl))] <- "highlight"
p2 <- ggplot(df2)+
  geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T) +
  geom_point(data = subset(df2, label == "highlight"), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df2$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="vQTL",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df2, label == "highlight")$pos / 1000000, y = subset(df2, label == "highlight")$logp, label = lable_df$variant_id, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

save(lable_df,df1,df2, file=paste0(Otdir, tissue, "_",pheno1, "_mahattan.plot.Rdata"))
Otplot <- paste0(Otdir, tissue, "_",pheno1, "_plot.pdf")
ggsave(Otplot, p2/p1,width=5,height=5)

box_df = fread(paste0(Otdir,prefix0 ,".pheno_geno.plot.txt.gz"))
p3 <- ggplot(box_df,aes(x=factor(gt),y=pheno_adj))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  geom_point(alpha = 0.7, size = 1, width = 0.1)+
  xlab(box_df$variant_id[1])+
  ylab(paste0("Normalized expression ","(", pheno1,")"))+
  theme_classic()

box_df1 = fread(paste0(Otdir,prefix1 ,".pheno_geno.plot.txt.gz"))
p4 <- ggplot(box_df1,aes(x=factor(gt),y=pheno_adj))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  geom_point(alpha = 0.7, size = 1, width = 0.1)+
  xlab(box_df1$variant_id[1])+
  ylab(paste0("Normalized expression ","(", pheno1,")"))+
  theme_classic()
Otplot1 <- paste0(Otdir, tissue, "_",pheno1, "_boxplot_eqtl.pdf")
ggsave(Otplot1, p4/p3,width=3,height=6)
