library("data.table")
setwd("/juice_data/zhh_independent/ZHH_update/vQTL")
info<-fread("/juice_data/zhh_independent/gene_annotion.txt")
tissue="Duodenum"
#tissue_list<-readLines("tiss.txt")
bed_file=paste0("/ossfs/bed/",tissue,".expr_tmm_inv.bed.gz")
bed_df<-fread(bed_file)
gene_info<-data.frame(bed_df[,c(1,4,3,4)],info$strand[match(bed_df$gene_id,info$gene_name)])
otdir=paste0("./pheno/",tissue,"/")
dir.create(otdir)
write.table(gene_info,paste0(otdir,tissue,".eed.opi"),quote=F,row.names=F,col.names=F,sep="\t")

expr<-data.frame(FID=0,IID=names(bed_df)[-c(1:4)],t(bed_df[,-c(1:4)]))
names(expr)[-c(1:2)]<-bed_df$gene_id
write.table(expr,paste0(otdir,tissue,".expr_tmm_inv.txt"),quote=F,row.names=F,sep="\t")

sex_info <- read.table(paste0("/juice_data/splicing_temp/covar/",tissue,"_known.covar.list"),head=T)
osca="/soft/osca-0.46.1-linux-x86_64/osca"
info_df <- data.frame(FID=0,IID=names(sex_info),ID1=0,ID2=0,sex=as.vector(t(sex_info)))
write.table(info_df,paste0(otdir,tissue,".eed.oii"),quote=F,row.names=F,col.names=F,sep="\t")
fam <- data.frame(info_df,trait=-9)
otfile=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tissue,"/",tissue,".autosomal.geno.fam")
write.table(fam,otfile,quote=F,row.names=F,col.names=F,sep="\t")
