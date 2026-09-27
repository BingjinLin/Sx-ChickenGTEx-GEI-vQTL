library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(grid)
library(gridExtra)
library(dendextend)
setwd("/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/result")
lfsr <- readRDS("lfsr_m.s.RDS")
pm<-readRDS("pm_m.s.RDS")
tissue_list=colnames(lfsr)
cor_res<-data.frame(tissue=tissue_list)
for (tissue in tissue_list) {
other_tis<-setdiff(tissue_list,tissue)
tis_cor<-c(1)
for (tis in other_tis) {
lfsr_tis<-lfsr[,c(tissue,tis)]
pm_tis<-pm[,c(tissue,tis)]
pm_sig <- pm_tis[apply(lfsr_tis < 0.05, 1, any), ]
pm_sig_corr <- cor(pm_sig, method = "spearman")
tis_cor<-c(tis_cor,pm_sig_corr[1,2])
}
tmp_res<-data.frame(tissue=c(tissue,other_tis),tis_cor)
names(tmp_res)[2]<-tissue
cor_res<-merge(cor_res,tmp_res)
}
gene_cor_df <- gather(cor_res, key = "tissue", value = "cor_eQTL")

library(tidyverse)
pm_sig_corr<-as.matrix(column_to_rownames(cor_res,"tissue"))

hc <- hclust(dist(pm_sig_corr, method = "euclidean"), method = "average")
od <- hc$order
pm_sig_corr <- pm_sig_corr[od, od]
row_hc <- hclust(dist(pm_sig_corr, method = "euclidean"), method = "average")
column_hc <- hclust(dist(t(pm_sig_corr), method = "euclidean"), method = "average")
row_dend <- as.dendrogram(row_hc)
column_dend <- as.dendrogram(column_hc)

tis_col <- read.table("/soft/zhh/vQTL/OSCA/tissue_colors.txt",header = T,sep = "\t",comment.char = "")
tis_col$tissue<-gsub(" ",".",tis_col$tissue)
all(colnames(pm_sig_corr) %in% tis_col$tissue)
tissue_col<-tis_col$colors
names(tissue_col)<-tis_col$tissue[match(tissue_col, tis_col$colors)]

###Adding left annotation
ha <- rowAnnotation(
  p = anno_density(
    pm_sig_corr, 
    type = "violin", 
    gp = gpar(fill = tissue_col[rownames(pm_sig_corr)]),
    border = FALSE,
    width = unit(3, "cm"),
  )
)
hf <- rowAnnotation(
  tissue1 = anno_points(
    x = rep(1, length(rownames(pm_sig_corr))),  # 使用常数来创建点
    axis = FALSE,
    gp = gpar(col = tissue_col[rownames(pm_sig_corr)], fill = tissue_col[rownames(pm_sig_corr)], fontsize = 10),  # 设置颜色和填充色
    pch = 16,  # 设置点的形状为圆形
    size = unit(3, "mm"),  # 设置点的大小
    border = FALSE,
    width = unit(0.25, "cm")
  ))

left_annotation <- c(ha, hf)

###Adding bottom annotation
bottom_annotation <- HeatmapAnnotation(
  tissue2 = anno_points(
    x = rep(1, length(colnames(pm_sig_corr))),  # 使用常数来创建点
    axis = FALSE,
    gp = gpar(col = tissue_col[colnames(pm_sig_corr)], fill = tissue_col[colnames(pm_sig_corr)], fontsize = 10),  # 设置颜色和填充色
    pch = 16,  # 设置点的形状为圆形
    size = unit(3, "mm"),  # 设置点的大小
    border = FALSE, 
    height = unit(0.25, "cm")
  ))

###plotting heatmap
ht <- Heatmap(
  pm_sig_corr,
  name = "cis-vQTL",
  cluster_rows = row_dend,   # 传递自定义行聚类树
  cluster_columns = column_dend,  # 传递自定义列聚类树
  show_row_dend = FALSE,
  show_column_dend = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  col = colorRamp2(c(0.4,0.5, 1), c("purple","white","firebrick3")),
  column_dend_side = "bottom",
  bottom_annotation = bottom_annotation,
  left_annotation = left_annotation,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(direction = "vertical"),
  rect_gp = gpar(type = "none"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (i > j) {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = fill))
    } else if (i == j) {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = NA))
      grid.text(
        colnames(pm_sig_corr)[j],
        x, 
        y, 
        gp = gpar(col = "black", fontsize = 9),
        just = "left"
      )
    } else {
      grid.rect(x, y, width, height, gp = gpar(col = NA, fill = NA))
    }
  },
  width = unit(ncol(pm_sig_corr)* 0.4, "cm"),  # 设置热图的宽度
  height = unit(nrow(pm_sig_corr)* 0.4, "cm")  # 设置热图的高度
)

pdf("cis-vQTL_heatmap_tissue-sharing.pdf", width = 12, height = 12)
draw(ht)
dev.off()

###############################################
###===============PICTURE2==================###
setwd("/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/result")
lfsr <- readRDS("lfsr_m.s.RDS")
lfsr_sig <- lfsr[apply(lfsr < 0.05, 1, any), ]
lfsr_count <- data.frame(count=apply(lfsr_sig, 1, function(x) sum(x < 0.05)))
breaks <- c(0, 5, 10, 15, 20, 25, 30, Inf)  # 使用 Inf 表示大于 30
labels <- c('1-5', '6-10', '11-15', '16-20', '21-25', '26-30', '>30')
lfsr_count$group <- cut(lfsr_count$count, breaks = breaks, labels = labels, right = TRUE)
write.table(lfsr_count,"lfsr_tissue_stat.txt",quote=F,sep="\t")
#####################
Otdir="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/result/Enrichment/"
dir.create(Otdir)
lfsr_count$SNP <- sub(".*,", "", row.names(lfsr_count))
spec_lfsr <- lfsr_count[lfsr_count$count==1,]
spec_file <- paste0(Otdir,"tis.specific_var.list")
write.table(spec_lfsr$SNP,spec_file, quote=F,row.names=F,col.names=F)

share_lfsr <- lfsr_count[lfsr_count$count>20,]
share_file <- paste0(Otdir,"tis.shared_var.list")
write.table(share_lfsr$SNP,share_file, quote=F,row.names=F,col.names=F)

spec_lfsr_sig <- lfsr_sig[rownames(spec_lfsr),]
spec_lfsr_sig$hit_col <- apply(spec_lfsr_sig, 1, function(x) {
  idx <- which(x < 0.05)
  if (length(idx) == 1) colnames(spec_lfsr_sig)[idx] else NA
})
spec_lfsr_sig$SNP <- sub(".*,", "", row.names(spec_lfsr_sig))
split_spec<-split(spec_lfsr_sig,spec_lfsr_sig$hit_col)

lapply(names(split_spec), function(group_name) {
  file_name <- paste0(Otdir, group_name, ".specify_var.list")
  write.table(split_spec[[group_name]]$SNP, file = file_name,quote=F,row.names=F,col.names=F)
}
)

Otdir_1="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/result/Enrichment/Group/"
dir.create(Otdir_1)

split_count<-split(lfsr_count,lfsr_count$count)

lapply(names(split_count), function(group_name) {
  file_name <- paste0(Otdir_1, "tissue.share", group_name, "._var.list")
  write.table(split_count[[group_name]]$SNP, file = file_name,quote=F,row.names=F,col.names=F)
}
)

#####################
count_table <- table(lfsr_count$group)
count_proportion <- prop.table(count_table)
count_df <- data.frame(
  count_group = names(count_proportion),
  proportion = as.numeric(count_proportion)
)
count_df$count_group<-factor(count_df$count_group,levels=labels)
count_df$type <- "vQTL"
library(ggplot2)
p1<-ggplot(count_df, aes(x = count_group, y = proportion)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black",width=.4) +
  labs(x = "Tissue with LFSR < 0.05", y = "Fraction of total cis-eQTL") +
   scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  theme(text=element_text(size=17))
ggsave("lfsr_tissue_stat.pdf", p1, width = 5, height = 5)

lfsr_eqtl <- readRDS("/ossfs/ZHH_update/ZHH_update/Tissue_sharing/result/output_strong_paris/lfsr_m.s.RDS")
lfsr_sig_e <- lfsr_eqtl[apply(lfsr_eqtl < 0.05, 1, any), ]
lfsr_count_e <- data.frame(count=apply(lfsr_sig_e, 1, function(x) sum(x < 0.05)))
breaks <- c(0, 5, 10, 15, 20, 25, 30, Inf)
labels <- c('1-5', '6-10', '11-15', '16-20', '21-25', '26-30', '>30')
lfsr_count_e$group <- cut(lfsr_count_e$count, breaks = breaks, labels = labels, right = TRUE)
count_table_e <- table(lfsr_count_e$group)
count_proportion_e <- prop.table(count_table_e)
count_df_e <- data.frame(
  count_group = names(count_proportion_e),
  proportion = as.numeric(count_proportion_e)
)
count_df_e$count_group<-factor(count_df_e$count_group,levels=labels)
count_df_e$type <- "eQTL"

merge_count <- rbind(count_df, count_df_e)
library(ggplot2)
pic<-ggplot(merge_count, aes(x = count_group, y = proportion, fill=type)) +
  geom_bar(stat = "identity", color = "black",width=.7,,position = position_dodge(width=0.8)) +
  scale_fill_manual(values = c("#F4F0B6","#D3BCED"))+
  labs(x = "Tissue with LFSR < 0.05", y = "Fraction of total cis-vQTL / cis-eQTL") +
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  theme(text=element_text(size=17))
ggsave("lfsr_tissue_stat_ev.pdf", pic, width = 8, height = 5)

count_table_1 <- table(lfsr_count$count)
count_proportion <- prop.table(count_table_1)
count_df <- data.frame(
  count_group = names(count_proportion),
  proportion = as.numeric(count_proportion)
)
count_df$count_group<-factor(count_df$count_group,levels=seq(1,32))
library(ggplot2)
p2<-ggplot(count_df, aes(x = count_group, y = proportion)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black",width=.45) +
  labs(x = "Tissue with LFSR < 0.05", y = "Fraction of total cis-eQTL") +
   scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  theme(text=element_text(size=16))
ggsave("lfsr_tissue_stat.32.pdf", p2, width = 8, height = 5)

###############################################
###===============PICTURE3==================###
lfsr_sig <- lfsr[apply(lfsr < 0.05, 1, any), ]
lfsr_count <- data.frame(count=apply(lfsr_sig, 1, function(x) sum(x < 0.05)))
split_rownames<-strsplit(rownames(lfsr_count),split = ",")
split_df <- as.data.frame(do.call(rbind, split_rownames))
colnames(split_df) <- c("SNP", "Gene")
lfsr_df <- cbind(split_df, lfsr_count)

library("data.table")
map<-fread("/home/hunau/Linbingjin/chicken_GTEx/Chicken.SNP_Autosomal.recode.map")
lfsr_df$POS<-map$V4[match(lfsr_df$SNP,map$V2)]
gene_info<-fread("/home/hunau/Linbingjin/chicken_GTEx/analyze_1/04_molQTLmapping/Input/TSS_annot.txt")
lfsr_df$TSS<-gene_info$TSS[match(lfsr_df$Gene, gene_info$gene_name)]
lfsr_df$dist<-abs(lfsr_df$POS-lfsr_df$TSS)/1000
freq<-fread("/home/hunau/Linbingjin/chicken_GTEx/Chicken.SNP.phased_Allele_Frequency.frq")
lfsr_df$MAF<-freq$MAF[match(lfsr_df$SNP,freq$SNP)]
write.table(lfsr_df,"cis-eQTL_lfsr_0.05_Info.txt",quote=F,sep="\t")

top_eQTL<-fread("all_top_cis-eQTL.txt")
lfsr_df$tissue<-top_eQTL$tissue[match(row.names(lfsr_df), top_eQTL$gene_SNP)]
tissue_list=unique(lfsr_df$tissue)
aFC_dir="/mnt/isidk43t/zhanghaihan/output_test/Liver/update/04_molQTLmapping/aFC/result/eQTL_aFC/"
new_res<-data.frame()
for (tissue in tissue_list){
tmp_df<-lfsr_df[lfsr_df$tissue==tissue,]
aFC<-fread(paste0(aFC_dir, tissue, ".OmiGA_leading.aFC.txt"))
tmp_df$aFC<-aFC$log2_aFC[match(tmp_df$Gene, aFC$pheno_id)]
new_res<-rbind(new_res,tmp_df)
}
write.table(new_res,"cis-eQTL_lfsr_0.05_Info.txt",quote=F,sep="\t")

###############################################
###===============PICTURE4==================###
library(dplyr)
tissue_list=colnames(lfsr)
cis_dir="/home/hunau/Linbingjin/chicken_GTEx/eQTL/cis_result/cis_res_fdr/"
geno_prefix="/home/hunau/Linbingjin/chicken_GTEx/Chicken.SNP_Autosomal.recode"
cor_res<-data.frame(tissue=tissue_list)
for (tissue1 in tissue_list) {
cis_file1=paste0(cis_dir, tissue1, ".cis_qtl_fdr0.05.txt.gz")
top_snp1<-fread(cis_file1)

other_tis<-setdiff(tissue_list,tissue1)
for (tissue2 in other_tis) {
cis_file2=paste0(cis_dir, tissue2, ".cis_qtl_fdr0.05.txt.gz")
top_snp2<-fread(cis_file2)

same_gene<-intersect(top_snp1$pheno_id,top_snp2$pheno_id)
same_snp<-intersect(top_snp1$variant_id,top_snp2$variant_id)
same_pair<-top_snp1 %>% inner_join(top_snp2,by=c("pheno_id","variant_id"))

test_snp<-union(top_snp1$variant_id,top_snp2$variant_id)
write.table(test_snp,"snps.list.tmp",quote=F,row.names=F,col.names=F)
system(paste0("plink --file ", geno_prefix, " --chr-set 39 --allow-extra-chr --keep-allele-order --extract snps.list.tmp --r square yes-really --out snps.test"))

}
tmp_res<-data.frame(tissue=c(tissue,other_tis),tis_cor)
names(tmp_res)[2]<-tissue
cor_res<-merge(cor_res,tmp_res)
}

