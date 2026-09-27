library("data.table")
tissue_list<-readLines("/soft/zhh/vQTL/OSCA/tissue.list")

###cieQTL
for(tiss in tissue_list){
Indir1<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tiss,"/zscore")
Infile1<-paste0(Indir1,"/",tiss,"_Vgene.lead.txt")
Vgene_df<-fread(Infile1)

ieGene_df<-data.frame()
ieDir="/juice_data/zhh_independent/ZHH_update/ieQTLmapping/ieGene_qval005"
iefile<-list.files(ieDir,full.name=T)
pattern <- paste(tiss, collapse = "|")
subset_file <- iefile[grep(pattern, iefile)]
for(f in subset_file){
df<-fread(f)
string1<-gsub(tiss,"",basename(f))
string2<-gsub(".cis_qtl_qval005.txt","",string1)
cell<-gsub("_","",string2)
tmp1<-data.frame(tissue=tiss,cell=cell,ieGene=df$pheno_id)
ieGene_df<-rbind(ieGene_df,tmp1)
}
if(nrow(ieGene_df) >0 ){
ieGene_df2 <- inner_join(ieGene_df, Vgene_df[,c("Chr","Probe")], by=c("ieGene"="Probe"))
write.table(ieGene_df2, file.path(Indir1, paste0(tiss,"_cieQTL_intersection.txt")), quote=F, sep="\t", row.names=F)
}
}

Indir="/juice_data/zhh_independent/ZHH_update/vQTL/res"
Infiles0 <- list.files(Indir,pattern="\\.vQTL_ieQTL\\.pph4$",full.names=T,recursive = TRUE)
res_cie<-data.frame()
for(f0 in Infiles0){
cie_df<-fread(f0)
names(cie_df) <- c("nsnp", "pph0","pph1","pph2","pph3", "pph4", "gene","chr","tissue","cell")
tmp_cie<-data.frame(cie_df[cie_df$pph4 > 0.8, c("gene","tissue","cell")])
res_cie<-rbind(res_cie,tmp_cie)
}
res_cie$pair <- paste(res_cie$tissue, res_cie$gene, sep=",")
res_cie$type <- "cie_coloc"

###sex-biaed eQTL
Indir="/juice_data/zhh_independent/ZHH_update/vQTL/res"
Infiles1 <- list.files(Indir,pattern="\\.vQTL_sex-eQTL\\.pph4$",full.names=T,recursive = TRUE)
res_sex<-data.frame()
for(f1 in Infiles1){
sex_df<-fread(f1)
names(sex_df) <- c("nsnp", "pph0","pph1","pph2","pph3", "pph4", "gene","chr","tissue")
tmp_sex<-data.frame(sex_df[sex_df$pph4 > 0.8, c("gene","tissue")])
res_sex<-rbind(res_sex,tmp_sex)
}
res_sex$pair <- paste(res_sex$tissue, res_sex$gene, sep=",")
res_sex$type <- "sex_coloc"

###exposure-biaed eQTL
Indir="/juice_data/zhh_independent/ZHH_update/vQTL/res"
Infiles2 <- list.files(Indir,pattern="\\.vQTL_GEQTL\\.coloc$",full.names=T, recursive = TRUE)
res_GE<-data.frame()
for(f2 in Infiles2){
GE_df<-fread(f2)
tmp_GE<-data.frame(GE_df[GE_df$PP.H4.abf > 0.8, c("gene","tissue","trait")])
res_GE<-rbind(res_GE,tmp_GE)
}
res_GE$pair <- paste(res_GE$tissue, res_GE$gene, sep=",")
res_GE$type <- "GE_coloc"
