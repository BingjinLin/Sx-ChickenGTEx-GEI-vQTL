library(data.table)

ARGS <- commandArgs(trailingOnly = TRUE)
tis = ARGS[1]


extract_df<-fread("/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs/nominal_pairs.combined_signifpairs.txt.gz")
extract_df$pair_id<-paste0(extract_df$Gene,",",extract_df$SNP)
chr_list<-unique(extract_df$Chr)
for( chr in chr_list){
dir_output="/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs"
dir_nominal="/juice_data/zhh_independent/ZHH_update/vQTL/res"
Infile<-paste0(dir_nominal,"/",tis,"/zscore/",tis,".vqtl_chr",chr,".all.query.gz")
if ( file.exists(Infile)){
Indat<-fread(Infile)
Indat$pair_id<-paste0(Indat$Gene,",",Indat$SNP)
extract_df$beta[extract_df$Chr==chr] <- Indat$b[match(extract_df$pair_id[extract_df$Chr==chr], Indat$pair_id)]
extract_df$beta_se[extract_df$Chr==chr] <- Indat$SE[match(extract_df$pair_id[extract_df$Chr==chr], Indat$pair_id)]
extract_df$beta[extract_df$Chr==chr] <- Indat$b[match(extract_df$pair_id[extract_df$Chr==chr], Indat$pair_id)]
}
}
Otfile<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/tissue_sharing/Random_pairs/",tis,"_nominal_pairs.extracted_pairs.txt")
write.table(extract_df, Otfile,sep="\t",quote=F,row.names=F)
system(paste0("bgzip ",Otfile))
