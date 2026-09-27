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
type=ARGS[6]
interaction=ARGS[7]


Otdir="/juice_data/zhh_independent/ZHH_update/vQTL/res/Plot/colocation/example/vQTL_mGWAS_ieQTL/"

genodir=paste0("/ossfs/ZHH_update/ZHH_update/Genotype/",tissue)
sample=read.table(paste0(genodir,"/",tissue,".autosomal.geno.fam"))
sample_size <- as.numeric(nrow(sample))

vDir<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore_update/")
file1 = paste0(vDir, tissue, ".vqtl_chr",chr, ".all.query.gz")
nominal1 <- fread(file1, header = T)
nominal1.1 <- nominal1[nominal1 $ Probe == pheno1, ]
#nominal1.1$dist <- abs(nominal1.1$BP-nominal1.1$Probe_bp)
#nominal1.1 <- nominal1.1[nominal1.1$dist < 1e6, ]
Otfile2=paste0(Otdir, tissue, "_",pheno1,"_vqtl_plot.txt")
write.table(nominal1.1, Otfile2, row.names=F, quote=F, sep="\t")

Gdir= "/ossfs/6.mGWAS/mGWAS_inde_loci/"
file2 <- paste0(Gdir, loci,".loci.txt")
GWAS = fread(file2 , header = T)
merged_data <- merge(nominal1.1, GWAS, by="SNP", all=FALSE, suffixes=c("_vqtl","_gwas"))

  df_test <- as.data.frame(merged_data)

  df1_coloc = list(beta=df_test$b_vqtl, pvalues=df_test$p_vqtl,
                   snp=df_test$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test$Freq_vqtl <= 0.5, df_test$Freq_vqtl, 1-df_test$Freq_vqtl)))
  df2_coloc = list(beta=df_test$b_gwas, pvalues=df_test$p_gwas,
                   snp=df_test$SNP, type="quant", N=999,
                   MAF=as.numeric(ifelse(df_test$Freq_gwas <= 0.5, df_test$Freq_gwas, 1-df_test$Freq_gwas)))
  my.res <- coloc.abf(dataset1=df1_coloc, dataset2=df2_coloc)
  results=as.data.frame(my.res$results)
  hc_snp <- results$snp[results$SNP.PP.H4==max(results$SNP.PP.H4)][1]

iefile <- switch(type,
  "Cell" = {
    cell_type <- interaction
    ieDir <- file.path("/ossfs/ieQTL_mapping/DWLS_Update", tissue, cell_type)
    paste0(ieDir, "/", cell_type,
           ".cis_qtl_pairs.", chr, ".txt.gz")
  },
  "Growth" = {
    ieDir <- file.path("/ossfs/GxE/veQTL_Growth", tissue, interaction)
    paste0(ieDir, "/", tissue, "_", interaction,
           ".cis_qtl_pairs.", chr, ".txt.gz")
  },
  "Slaughter" = {
    ieDir <- file.path("/ossfs/GxE/veQTL", tissue, interaction)
    paste0(ieDir, "/", tissue, "_", interaction,
           ".cis_qtl_pairs.", chr, ".txt.gz")
  },
  "Sex" = {
    ieDir <- file.path("/juice_data/zhh_independent/ZHH_update/GxE/tensorQTL", tissue)
    paste0(ieDir, "/", tissue, ".cis_qtl_pairs.", chr, ".txt.gz")
  },
  stop("Unknown type: ", type)
)

nominal2 <- fread(iefile, header = T)
names(nominal2)[1] <- "pheno_id"
nominal2.1 <- nominal2[nominal2$pheno_id == pheno1, ]
merged_data_ie <- merge(nominal1.1, nominal2.1, by.x="SNP", by.y="variant_id", all=FALSE, suffixes=c("_vqtl","_ieqtl"))

  df_test_ie <- as.data.frame(merged_data_ie)

  df1_coloc_ie = list(beta=df_test_ie$b, pvalues=df_test_ie$p,
                   snp=df_test_ie$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test_ie$Freq <= 0.5, df_test_ie$Freq, 1-df_test_ie$Freq)))
if(type=="Sex"){
  df2_coloc_ie = list(beta=df_test_ie$b_gi, pvalues=df_test_ie$pval_gi,
                   snp=df_test_ie$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test_ie$af <= 0.5, df_test_ie$af, 1-df_test_ie$af)))
}else{
  df2_coloc_ie = list(beta=df_test_ie$beta_g2, pvalues=df_test_ie$pval_g2,
                   snp=df_test_ie$SNP, type="quant", N=sample_size,
                   MAF=as.numeric(ifelse(df_test_ie$af <= 0.5, df_test_ie$af, 1-df_test_ie$af)))
}

  res_ie <- coloc.abf(dataset1=df1_coloc_ie, dataset2=df2_coloc_ie)
  results_ie=as.data.frame(res_ie$results)
  hc_snp_ie <- results_ie$snp[results_ie$SNP.PP.H4==max(results_ie$SNP.PP.H4)][1]

#merge_res <- inner_join(results,results_ie,by=c("snp"))
###eqtl plot###
qtl_df <- data.frame(pheno_id = pheno1, variant_id=c(hc_snp,hc_snp_ie))
Otfile1=paste0(Otdir, tissue, "_",pheno1,"_eqtl.extract.txt")
write.table(qtl_df, Otfile1, row.names=F, quote=F, sep="\t")
prefix <- paste0(tissue, "_chr", chr, "_", pheno1,".eqtl")
script="/soft/zhh/OmiGA/omiga_eQTL_plot.sh"
omiga_cmd <- paste("bash",script, tissue, Otfile1, Otdir, prefix)
system(omiga_cmd)

###ieqtl plot###
geno_file=paste0(genodir,"/",tissue,".autosomal.geno")
out_geno=paste0(Otdir,tissue,"_", hc_snp_ie,"_genotype")
plink_cmd <- paste0("plink --bfile ", geno_file, " --chr-set 39 --keep-allele-order --snp ", hc_snp_ie, " --recode tab --out ",out_geno)
system(plink_cmd)
geno_df<-fread(paste0(out_geno,".ped"))
names(geno_df)<-c("ID1","ID2","info1","info2","info3","trait","snp")
genotype<-geno_df$snp
genotype<-gsub(" ","|",genotype)
names(genotype)<-geno_df$ID2

###gwas plot###
Otfile3=paste0(Otdir, loci, "_gwas.txt")
#if(!file.exists(Otfile3)){
Gfile2 <- paste0("/ossfs/6.mGWAS/", trait, ".mlma")
G_df <- fread(Gfile2) %>% as.data.frame()
GWAS_SNP<- intersect(G_df$SNP, nominal1.1$SNP)
sig_file=paste0(Otdir,loci,".GWAS_sig-snp.list")
write.table(GWAS_SNP,sig_file,quote=F,row.names=F,col.names=F)
outfile0=paste0(Otdir,trait,",",hc_snp,"_genotype")
outfile1=paste0(Otdir,loci,".sig")
outfile2=paste0(Otdir,loci)
Otfile3=paste0(Otdir, loci, "_gwas.txt")
GWAS_geno<-"/juice_data/WGS/GWAS1/Metabolome/prepare_GRM/GWAS1_999_no_dup"
plink0<-paste0("plink --bfile ", GWAS_geno, " --chr-set 39 --keep-allele-order --recode A --snp ", hc_snp, " --out ", outfile0)
system(plink0)
plink1<-paste0("plink --bfile ", GWAS_geno, " --chr-set 39 --keep-allele-order --make-bed --extract ", sig_file, " --out ", outfile1)
system(plink1)
plink2<-paste0("plink --bfile ", outfile1, " --chr-set 39 --r2 --ld-snp ",hc_snp, " --ld-window-r2 0 --ld-window 2000000 --out ", outfile2)
system(plink2)
ld_score<-fread(paste0(outfile2,".ld"))
G_df_1 <- G_df[G_df$SNP %in% GWAS_SNP, ]
G_df_1$LD_r2<-ld_score$R2[match(G_df_1$SNP,ld_score$SNP_B)]
G_df_1 <- na.omit(G_df_1)
write.table(G_df_1,Otfile3,quote=F,row.names=F,sep="\t")
#}else{
#print(paste(Otfile3,"has been existed!"))
#}

library(ggplot2)
library(ggrastr)
library(patchwork)

###vqtl plot
snpr_df = fread(paste0(Otdir,prefix ,".snpr.plot.txt.gz"))
snpr_vdf <- snpr_df[snpr_df$variant_id==hc_snp,]
df2 = data.frame(SNP = nominal1.1$SNP, pos = nominal1.1$BP, logp = -log10(nominal1.1$p), ld = snpr_vdf$corr[match(nominal1.1$BP, snpr_vdf$pos)]^ 2)

p2 <- ggplot(df2)+
  rasterise(geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T),dpi=150) +
  geom_point(data = subset(df2, SNP%in%c(hc_snp,hc_snp_ie)), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df2$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="vQTL",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df2, SNP%in%c(hc_snp,hc_snp_ie))$pos / 1000000, y = subset(df2, SNP%in%c(hc_snp,hc_snp_ie))$logp, label = c(hc_snp,hc_snp_ie), hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

###ieqtl plot
snpr_idf <- snpr_df[snpr_df$variant_id==hc_snp_ie,]
df1 = data.frame(pos = df2$pos[match(nominal2.1$variant_id, df2$SNP)], logp = -log10(nominal2.1$pval_g2))
df1$ld = snpr_idf$corr[match(df1$pos, snpr_idf$pos)]^ 2
p1 <- ggplot(df1)+
  rasterise(geom_point(aes(x = pos / 1000000, y = logp, fill = ld), shape = 21, show.legend = T),dpi=150) +
  geom_point(data = subset(df1, ld == 1), aes(x = pos / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(df1$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="eQTL",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(df1, ld == 1)$pos / 1000000, y = subset(df1, ld == 1)$logp, label = hc_snp_ie, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

###gwas plot
G_df_1$logp = -log10(G_df_1$p)
p3 <- ggplot(G_df_1)+
  rasterise(geom_point(aes(x = bp / 1000000, y = logp, fill = LD_r2), shape = 21, show.legend = T),dpi=150) +
  geom_point(data = subset(G_df_1, LD_r2 == 1), aes(x = bp / 1000000, y = logp), shape = 23, fill = "#EB2326", size = 2.5) +
  theme_classic() +
  ylim(0,max(G_df_1$logp)*1.1)+
  ggplot2::binned_scale(aesthetics = "fill", scale_name = "custom", palette = ggplot2:::pal_binned(scales::manual_pal(values = c("#14128A", "#8CCCF0", "#6CBE44", "#F8A41C", "#EB2326"))), guide = "bins", breaks = c(0.2, 0.4, 0.6, 0.8), labels = function(x) { formatC(x, digits = 1, format = "f", flag = "#") }, limits = c(0, 1), show.limits = T, name = expression(r ^ 2)) +
  guides(fill = guide_colorsteps(show.limits = TRUE, ticks = T)) +
  labs(title="GWAS",x="",y=expression(-log[10](italic(P)))) +
  annotate("label", x = subset(G_df_1, LD_r2 == 1)$bp / 1000000, y = subset(G_df_1, LD_r2== 1)$logp, label = hc_snp, hjust = 1.1, size = 2, label.size=0, fill = alpha("white",0.5))

Otplot <- paste0(Otdir, interaction, "_",pheno1, "_", loci, "_plot.pdf")
ggsave(Otplot, p2/p1/p3,width=5.5,height=5)
save(df1,df2,G_df_1,hc_snp_ie,hc_snp,file=paste0(Otdir, interaction, "_",pheno1, "_", loci, "_plot.Rdata"))

plot_dt = fread(paste0(Otdir,prefix, ".pheno_geno.plot.txt.gz"))
plot_vdt <- plot_dt[plot_dt$variant_id==hc_snp,]
plot1 = ggplot(plot_vdt,aes(x=gt_base,y=pheno_adj))+
  geom_point()+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  xlab(plot_dt$variant_id[1])+
  ylab(paste0("Normalized expression (",pheno1,")"))+
  theme_classic()
  
mebo_level <- fread("/juice_data/WGS/GWAS1/Metabolome/vmQTL/trait_data/Metabolome.pheno.txt")
#mebo_level <- fread("/juice_data/WGS/GWAS1/Metabolome/heritability/Metabolome_data.txt")
mebo_label <- read.table("/juice_data/WGS/GWAS1/Metabolome/heritability/Metabolome_label.list")
sample_label <- read.table("/juice_data/WGS/GWAS1/Metabolome/heritability/sample_info.txt")
change_lable <- read.table("/soft/zhh/vQTL/Coloc/Metabolome_lable_change.list",head=T)
trait_lable <- change_lable$OSCA_lable[change_lable$GWAS_lable==trait]
snp_geno <- fread(paste0(outfile0,".raw"))
names(snp_geno)[7] <- "Genotype"
allele0 <- G_df_1$A1[G_df_1$SNP==hc_snp]
allele1 <- G_df_1$A2[G_df_1$SNP==hc_snp]
snp_geno$Genotype[snp_geno$Genotype==0] <- paste0(allele0,"|",allele0)
snp_geno$Genotype[snp_geno$Genotype==1] <- paste0(allele0,"|",allele1)
snp_geno$Genotype[snp_geno$Genotype==2] <- paste0(allele1,"|",allele1)
snp_geno$Genotype <- factor(snp_geno$Genotype)
mebo_select <- mebo_level[, c("IID", trait_lable), with = FALSE]
names(mebo_select)[2] <- "PHENOTYPE"
snp_geno$PHENOTYPE <- mebo_select$PHENOTYPE[match(snp_geno$IID, mebo_select$IID)]

model <- aov(PHENOTYPE ~ Genotype, data = snp_geno) 
summary(model)

library(dplyr)
comparisons <- list(
  c(levels(snp_geno$Genotype)[1],levels(snp_geno$Genotype)[2]),
  c(levels(snp_geno$Genotype)[2],levels(snp_geno$Genotype)[3])
)

results <- lapply(comparisons, function(comp){
  df <- subset(snp_geno, Genotype %in% comp)
  
  test <- wilcox.test(PHENOTYPE  ~ Genotype, data = df)
  
  data.frame(
    group1 = comp[1],
    group2 = comp[2],
    p.value = test$p.value
  )
})

results <- do.call(rbind, results)
results
mebo_level[[trait_lable]] -> mebo_select
names(mebo_select) <- mebo_level$IID
snp_geno$PHENOTYPE <- mebo_select[match(snp_geno$IID, names(mebo_select))]
plot2 <- ggplot(snp_geno,aes(x=Genotype,y=PHENOTYPE))+
  geom_point()+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  xlab(plot_dt$variant_id[1])+
  ylab(trait)+
  theme_classic()

Otbox <- paste0(Otdir, loci,"_",interaction, ",",pheno1, "_box.pdf")
ggsave(Otbox, plot1/plot2,width=2.5,height=5)

###ieQTL plot
if (type == "Cell") {
cell_dat<-fread(file.path("/juice_data/zhh_independent/ZHH_update/03_DWLS_Results_Use/DWLS_raw_74_TCPairs_MedianPro0.1",tissue,paste0(cell_type,"_result.txt"))) %>% as.data.frame()
#cell_dat<-cell_df[,names(cell_df)%in%c("sample_id",cell_type)]
names(cell_dat)<-c("sample_id","cell_prop")
plot_idt <- plot_dt[plot_dt$variant_id==hc_snp_ie,]
geno_df$ID1 <- substr(geno_df$ID2,3,nchar(geno_df$ID2))
plot_ie<-data.frame(TMM_inv=as.numeric(plot_idt$pheno_adj),genotype=plot_idt$gt_base, cell_prop=cell_dat$cell_prop[match(geno_df$ID1,cell_dat$sample_id)])
pp1 <- ggplot(plot_ie,aes(x=cell_prop,y=TMM_inv,color=factor(genotype))) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.2) +
  geom_smooth(aes(group = genotype), method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values=c("#f9b376","#58a0cd","#aed8a3")) +
  labs(x = paste(cell_type,"Enrichment"), y = paste0(pheno1," expresssion"),color="",title=paste(tissue,pheno1,hc_snp_ie)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "serif", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

pp2 <- ggplot(plot_ie,aes(x=factor(genotype),y=TMM_inv,color=factor(genotype))) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.2) +
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  scale_color_manual(values=c("#f9b376","#58a0cd","#aed8a3")) +
  labs(x = hc_snp_ie, y =paste0(pheno1," expresssion"),color="",title=paste(tissue,pheno1)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "serif", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )
pic_file<-paste0(Otdir,cell_type,".",pheno1,".vQTL.ieQTL.pheno_group.pdf")  
ggsave(pic_file,pp1+pp2,width = 10,height = 5)

} else if (type == "Growth") {
trait_df<-fread("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Growth_trait.txt") %>% as.data.frame()
trait_dat<-trait_df[,names(trait_df)%in%c("IID",interaction)]
trait_dat$sample_id <- paste0(substr(names(genotype),1,2),trait_dat$IID)
names(trait_dat)[2]<-"trait_level"
plot_ie<-data.frame(TMM_inv=as.numeric(plot_dt$pheno_adj),genotype=genotype, trait_level=trait_dat$trait_level[match(names(genotype),trait_dat$sample_id)])
pp1 <- ggplot(plot_ie,aes(x=trait_level,y=TMM_inv,color=factor(genotype))) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.2) +
  geom_smooth(aes(group = genotype), method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values=c("#f9b376","#58a0cd","#aed8a3")) +
  labs(x = interaction, y = paste0(pheno1," expresssion"),color="",title=paste(tissue,pheno1,hc_snp_ie)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "serif", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

pp2 <- ggplot(plot_ie,aes(x=factor(genotype),y=TMM_inv,color=factor(genotype))) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.2) +
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  scale_color_manual(values=c("#f9b376","#58a0cd","#aed8a3")) +
  labs(x = hc_snp_ie, y =paste0(pheno1," expresssion"),color="",title=paste(tissue,pheno1)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "serif", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )
pic_file<-paste0(Otdir,tissue,".",interaction,".",pheno1,".vQTL.ieQTL.pheno_group.pdf")  
ggsave(pic_file,pp1+pp2,width = 10,height = 5)

} else if (type == "Slaughter") {
trait_df<-fread("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Slaughter_trait.txt") %>% as.data.frame()
trait_dat<-trait_df[,names(trait_df)%in%c("IID",interaction)]
trait_dat$sample_id <- paste0(substr(names(genotype),1,2),trait_dat$IID)
names(trait_dat)[2]<-"trait_level"
plot_ie<-data.frame(TMM_inv=as.numeric(plot_dt$pheno_adj),genotype=genotype, trait_level=trait_dat$trait_level[match(names(genotype),trait_dat$sample_id)])
pp1 <- ggplot(plot_ie,aes(x=trait_level,y=TMM_inv,color=factor(genotype))) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.2) +
  geom_smooth(aes(group = genotype), method = "lm", se = FALSE, linewidth = 1) +
  scale_color_manual(values=c("#f9b376","#58a0cd","#aed8a3")) +
  labs(x = interaction, y = paste0(pheno1," expresssion"),color="",title=paste(tissue,pheno1,hc_snp_ie)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "serif", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )

pp2 <- ggplot(plot_ie,aes(x=factor(genotype),y=TMM_inv,color=factor(genotype))) +
  geom_jitter(alpha = 0.7, size = 2, width = 0.2) +
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  scale_color_manual(values=c("#f9b376","#58a0cd","#aed8a3")) +
  labs(x = hc_snp_ie, y =paste0(pheno1," expresssion"),color="",title=paste(tissue,pheno1)) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "serif", face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14),
    legend.position = "bottom"
  )
pic_file<-paste0(Otdir,tissue,".",interaction,".",pheno1,".vQTL.ieQTL.pheno_group.pdf")  
ggsave(pic_file,pp1+pp2,width = 10,height = 5)

} else if (type == "Sex") {
plot_dt$IID <- sample$V2
plot_dt$Sex_id <- sample$V5
plot_dt$Sex <- ifelse(plot_dt$Sex_id=="1","Male","Female")
ppic <- ggplot(plot_dt,aes(x=gt,y=pheno_adj))+
  geom_violin(linewidth=0.25)+
  geom_boxplot(linewidth=0.25, width=0.25, outlier.size = 0.5)+
  facet_wrap(facets = vars(Sex))+
  xlab(plot_dt$variant_id[1])+
  ylab(paste0("Normalized expression (",pheno1,")"))+
  theme_classic()
pic_file <- paste0(Otdir, tissue, "_",pheno1, ".vqtl_sxqtl_boxplot.pdf")
ggsave(pic_file, ppic,width=6,height=2.5)
}
}