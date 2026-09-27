library(data.table)
library(dplyr)

exposure_df2 <- read.table("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Slaughter_trait.txt",head=T,sep="\t")
datadir="/ossfs/bed"
Otdir="/data/bed"
exposure_df2$SEX[exposure_df2$SEX=="2"] <- "0"

tissue_list <- readLines("/soft/zhh/vQTL/OSCA/tissue.list")
for(tissue in tissue_list){
print(paste("Processing on",tissue))
phenofile=file.path(datadir,paste0(tissue,".expr_tmm_inv.bed.gz"))
pheno_df <- fread(phenofile)
info_df <- pheno_df[,c(1:4)]
pheno_data <- data.frame(t(pheno_df[,-c(1:4)]))
pheno_data$ID <- substr(rownames(pheno_data),3,nchar(rownames(pheno_data)))
pheno_data$Etime <- exposure_df2$Etime[match(pheno_data$ID, exposure_df2$IID)]
pheno_data$Etime <- factor(pheno_data$Etime,levels=c("1","2","3"))
pheno_list <- split(pheno_data, pheno_data$Etime)
for(name in names(pheno_list)){
Otfile=file.path(Otdir,paste0(tissue,"_Etime",name,".expr_tmm_inv.bed"))
Otdata <- pheno_list[[name]]
if(nrow(Otdata)>0){
Otdata_df <- as.data.frame(Otdata)
Otdata_1 <- Otdata_df[, !(colnames(Otdata_df) %in% c("ID","Etime"))]
Otdata_1 <- data.frame(cbind(info_df,t(Otdata_1)))
names(Otdata_1)[1] <- "#Chr"
write.table(Otdata_1, Otfile, quote=F, row.names=F,sep="\t")
bgzip_cmd <- paste("bgzip", Otfile)
system(bgzip_cmd)
tabix_cmd <- paste0("tabix -p bed ", Otfile,".gz")
system(tabix_cmd)

covar_file=paste0("/juice_data/splicing_temp/covar/",tissue,"_known.covar.list")
if (file.exists(covar_file)){
covar_df <- fread(covar_file,head=T) %>% as.data.frame()
target_cols <- c("V1", rownames(Otdata))
matched_cols <- colnames(covar_df)[match(target_cols,colnames(covar_df))]
covar_new <- covar_df[, matched_cols, drop = FALSE]
names(covar_new)[1] <- ""
if (length(unique(unlist(covar_new))) >1){
Otfile2 <- file.path(Otdir,paste0(tissue,"_Etime",name,"_known.covar.list"))
write.table(covar_new, Otfile2, row.names=F,quote=F,sep="\t")
}
}
geno_dir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tissue)
sample_list=data.frame(Otdata[,1:2])
Otfile3 <- file.path(geno_dir,paste0(tissue,"_Etime",name,"_sample.list"))
write.table(sample_list, Otfile3, quote=F,row.names=F,col.names=F,sep="\t")
}
}
}

Slaughter_df <- fread("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Slaughter_trait.txt") %>% as.data.frame()
tissue_list <- readLines("/soft/zhh/vQTL/OSCA/tissue.list")
for(tissue in tissue_list){
fam_file=file.path("/juice_data/zhh_independent/ZHH_update/Genotype",tissue,paste0(tissue,"_vqtl.autosomal.geno.fam"))
fam_df <- read.table(fam_file)
prefix <- substr(fam_df$V2,1,2)
new_tf <- Slaughter_df
new_tf$IID <- paste0(prefix, new_tf$IID)
new_tf <- new_tf[match(fam_df$V2, new_tf$IID),]
Otfile <- paste0("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Tissue/", tissue, "_Slaughter_trait.txt")
write.table(new_tf, Otfile, quote=F, row.names=F,sep="\t")
}