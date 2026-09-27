# step1. output all top significant eQTL information from permutation results to a .txt file
library("data.table")
#SNP_info<-fread("/ossfs/tissue_sharing/MashR/GRCg7b_SNP.txt")
Indir="/juice_data/zhh_independent/ZHH_update/vQTL/res"
files<-list.files(Indir,pattern="_Vgene\\.lead\\.txt$",full.names=T, recursive = TRUE)
top_var<-data.frame()
for ( f in files ) {
data<-fread(f)
data_1<-data[, c("Chr", "SNP", "BP","Gene")]
top_var<-rbind(top_var, data_1)
}
top_var<-top_var[!duplicated(top_var),]
#names(top_var)<-c("chr","variant_id","pos","pheno_id")
otdir1<-"/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/top_pairs/"
dir.create(otdir1)
write.table(top_var,paste0(otdir1,"top_pairs.combined_signifpairs.txt"),sep="\t",quote=F,row.names=F)
