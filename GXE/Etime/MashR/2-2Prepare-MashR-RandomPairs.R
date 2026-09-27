ARGS <- commandArgs(trailingOnly = TRUE)
tissue = ARGS[1]
library("data.table")
info<-fread("/ossfs/tissue_sharing/MashR/GRCg7b_SNP.txt")
dir_output=paste0("/data/Etime/",tissue,"/Random_pairs")
data<-fread(file.path(dir_output,"nominal_paires.shuf.txt"),head=F)
df<-data.frame(chr=info$V1[match(data$V2,info$V2)],variant_id=data$V2, pos=info$V3[match(data$V2,info$V2)], pheno_id=data$V1)
write.table(df, file.path(dir_output,"nominal_pairs.combined_signifpairs.txt"),sep="\t",quote=F,row.names=F)

