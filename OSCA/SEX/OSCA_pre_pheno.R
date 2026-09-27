library(data.table)

exposure_df2 <- read.table("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Slaughter_trait.txt",head=T,sep="\t")
datadir="/juice_data/zhh_independent/ZHH_update/vQTL/pheno"
exposure_df2$SEX[exposure_df2$SEX=="2"] <- "0"

tissue_list <- readLines("/soft/zhh/vQTL/OSCA/tissue.list")
for(tissue in tissue_list[-1]){
phenofile=file.path(datadir,tissue,paste0(tissue,".expr_tmm_inv.zscore"))
pheno_df <- fread(phenofile)
pheno_df$ID <- substr(pheno_df$IID,3,nchar(pheno_df$IID))
pheno_df$SEX<- exposure_df2$SEX[match(pheno_df$ID, exposure_df2$IID)]
pheno_df$SEX <- factor(pheno_df$SEX,levels=c("0","1"))
pheno_list <- split(pheno_df, pheno_df$SEX)
for(name in names(pheno_list)){
SEX <- ifelse(name=="0","female","male")
Otfile=file.path(datadir,tissue,paste0(tissue,"_",SEX,".expr_tmm_inv.zscore"))
Otdata <- pheno_list[[name]]
if(nrow(Otdata)>0){
Otdata_df <- as.data.frame(Otdata)
Otdata_1 <- Otdata_df[, !(colnames(Otdata_df) %in% c("ID","SEX"))]
write.table(Otdata_1, Otfile, quote=F, row.names=F,sep="\t")

sample_df<-data.frame(Otdata[,1:2],ID1=0,ID2=0,sexID=Otdata$SEX )
Otfile2=file.path(datadir,tissue,paste0(tissue,"_",SEX,".oii"))
write.table(sample_df, Otfile2, quote=F,row.names=F,col.names=F,sep="\t")

geno_dir=paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tissue)
sample_list=data.frame(Otdata[,1:2])
Otfile3 <- file.path(geno_dir,paste0(tissue,"_",SEX,"_sample.list"))
write.table(sample_list, Otfile3, quote=F,row.names=F,col.names=F,sep="\t")
}
}
}

library(dplyr)

female_ind <- exposure_df2$IID[exposure_df2$SEX==0]
male_ind <- exposure_df2$IID[exposure_df2$SEX==1]
# expression thresholds
count_threshold = 6
tpm_threshold = 0.1
sample_frac_threshold = 0.2
sample_count_threshold = 10

# Transform rows to a standard normal distribution
inverse_normal_transform = function(x) {
    qnorm(rank(x) / (length(x)+1))}

load("/ossfs/ZHH_update/ZHH_update/Expr_filtered.Rdata")
sample_info <- read.table("/ossfs/ZHH_update/ZHH_update/7969sample_Info.txt",head=T,sep="\t")
sample_info$tissue <- gsub(" ",".",sample_info$tissue)

for(tis in tissue_list[2:30]){
Tfile=paste0("/juice_data/zhh_independent/ZHH_update/TMM_matrix/", tis, "/", tis, ".expr_TMM.bed.gz")
TMM_df <- fread(Tfile)
TMM_df <- as.data.frame(TMM_df)
female_sam <- sample_info$sample[sample_info$tissue==tis & sample_info$sex=="female"]
male_sam <- sample_info$sample[sample_info$tissue==tis & sample_info$sex=="male"]
TMM_fdf <- TMM_df[, names(TMM_df) %in% c("#Chr", "start", "end", "gene_id", female_sam)]
TMM_mdf <- TMM_df[, names(TMM_df) %in% c("#Chr", "start", "end", "gene_id", male_sam)]

#keep the genes with >=0.1 tpm and >=6 read counts in >=20% samples.
expr_fdf = TPM[,names(TPM) %in% female_sam]
tpm_fth = rowSums (expr_fdf >= tpm_threshold)
expr_counts_f = Counts[rownames (expr_fdf),names(Counts) %in% female_sam]
count_fth = rowSums(expr_counts_f >= count_threshold)
ctrl1 = tpm_fth >= (sample_frac_threshold * ncol(expr_fdf))
ctrl2 = count_fth >= (sample_frac_threshold * ncol(expr_counts_f))
mask = ctrl1 & ctrl2
keep_gene_f <- rownames (expr_fdf)[mask]
TMM_fdf[TMM_fdf$gene_id %in% keep_gene_f,] -> TMM_fdf_1
TMM_fdf_inv = t(apply(TMM_fdf_1[,-c(1:4)], MARGIN = 1, FUN = inverse_normal_transform))
bed_f = data.frame(TMM_fdf_1[,c(1:4)], TMM_fdf_inv)
bed_f[,1] = as.numeric(bed_f[,1])
bed_f = bed_f[order(bed_f[,1],bed_f[,2]),]
colnames(bed_f)[1] = "#Chr"
#fam_f=read.table(paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis,"/",tis,"_female.autosomal.geno.fam"))
Otfile_f <- paste0("/juice_data/zhh_independent/ZHH_update/TMM_matrix/",tis,"/",tis,"_female.expr_tmm_inv.bed")
write.table(bed_f, Otfile_f, quote=F,row.names=F, sep="\t")
system(paste("bgzip", Otfile_f))
system(paste0("tabix -p bed ", Otfile_f,".gz"))

expr_mdf = TPM[,names(TPM) %in% male_sam]
tpm_mth = rowSums (expr_mdf >= tpm_threshold)
expr_counts_m = Counts[rownames (expr_mdf),names(Counts) %in% male_sam]
count_mth = rowSums(expr_counts_m >= count_threshold)
ctrl1 = tpm_mth >= (sample_frac_threshold * ncol(expr_mdf))
ctrl2 = count_mth >= (sample_frac_threshold * ncol(expr_counts_m))
mask = ctrl1 & ctrl2
keep_gene_m <- rownames (expr_mdf)[mask]
TMM_mdf[TMM_mdf$gene_id %in% keep_gene_m,] -> TMM_mdf_1
TMM_mdf_inv = t(apply(TMM_mdf_1[,-c(1:4)], MARGIN = 1, FUN = inverse_normal_transform))
bed_m = data.frame(TMM_mdf_1[,c(1:4)], TMM_mdf_inv)
bed_m[,1] = as.numeric(bed_m[,1])
bed_m = bed_m[order(bed_m[,1],bed_m[,2]),]
colnames(bed_m)[1] = "#Chr"
#fam_m=read.table(paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis,"/",tis,"_male.autosomal.geno.fam"))
Otfile_m <- paste0("/juice_data/zhh_independent/ZHH_update/TMM_matrix/",tis,"/",tis,"_male.expr_tmm_inv.bed")
write.table(bed_m, Otfile_m, quote=F,row.names=F, sep="\t")
system(paste("bgzip", Otfile_m))
system(paste0("tabix -p bed ", Otfile_m,".gz"))
}

#function
process_one_gene <- function(gene_values, covariates) {
  df <- data.frame(TMM_inv = gene_values, covariates)
  predictor_names <- names(df)[-1]
  
  formula_dynamic <- reformulate(termlabels = predictor_names, response = "TMM_inv")
  model <- lm(formula_dynamic, data = df)
  df$residuals <- residuals(model)

  mean_res <- mean(df$residuals)
  sd_res   <- sd(df$residuals)
  abs_z_score_residuals <- abs(df$residuals - mean_res) / sd_res
  cleaned_data <- df[abs_z_score_residuals <= 5, , drop = FALSE]
  
  final_data <- cleaned_data %>%
    mutate(
      z_score = (residuals - mean(residuals)) / sd(residuals)
    )
  rownames(final_data)<-row.names(cleaned_data)

  z_full <- rep(NA, nrow(covariates))
  names(z_full) <- rownames(covariates)
  z_full[rownames(final_data)] <- final_data$z_score
  return(z_full)
}

Otdir="/juice_data/zhh_independent/ZHH_update/vQTL/pheno"

for(tis in tissue_list[2:30]){
cat("Process on task: ", tis,"\n")
Infile_f <- paste0("/juice_data/zhh_independent/ZHH_update/TMM_matrix/",tis,"/",tis,"_female.expr_tmm_inv.bed.gz")
df_f <- fread(Infile_f)
expr_f <- data.frame(t(df_f[,-(1:4)]))
rownames(expr_f) = df_f$gene_id
covar_f <- paste0("/data/sex_stratified_eQTL/",tis,"/",tis,"_female.auto_covar")
covar_fdf = read.table(covar_f, header=T, row.names=1)
covariatesToUse_f <- data.frame(t(covar_fdf))
z_matrix_f <- sapply(1:ncol(expr_f), function(i) {
  process_one_gene(expr_f[, i], covariatesToUse_f)
})
colnames(z_matrix_f) <- df_f$pheno_id
na_per_colf <- colSums(is.na(z_matrix_f))
num_cols_with_na <- sum(na_per_colf > 0)
cat(tis," female num_cols_with_n: ", num_cols_with_na, "\n")

z_matrix_f_no_na <- z_matrix_f[, colSums(is.na(z_matrix_f)) == 0]

f_fam <- read.table(paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis,"/",tis,"_female.autosomal.geno.fam"))
otdata_f<-data.frame(FID=0,IID=f_fam$V2,z_matrix_f_no_na[f_fam$V2, ])
Otfile_f=file.path(Otdir,tis,paste0(tis,"_female.expr_tmm_inv.zscore"))
write.table(otdata_f, Otfile_f, quote=F,row.names=F,sep="\t")

Infile_m <- paste0("/juice_data/zhh_independent/ZHH_update/TMM_matrix/",tis,"/",tis,"_male.expr_tmm_inv.bed.gz")
df_m <- fread(Infile_m)
expr_m <- data.frame(t(df_m[,-(1:4)]))
rownames(expr_m) = df_m$gene_id
covar_m <- paste0("/data/sex_stratified_eQTL/",tis,"/",tis,"_male.auto_covar")
covar_mdf = read.table(covar_m, header=T, row.names=1)
covariatesToUse_m <- data.frame(t(covar_mdf))

z_matrix_m <- sapply(1:ncol(expr_m), function(i) {
  process_one_gene(expr_m[, i], covariatesToUse_m)
})
colnames(z_matrix_m) <- df_m$pheno_id
na_per_colm <- colSums(is.na(z_matrix_m))
num_cols_with_na <- sum(na_per_colm > 0)
cat(tis," male num_cols_with_n: ", num_cols_with_na, "\n")
z_matrix_m_no_na <- z_matrix_m[, colSums(is.na(z_matrix_m)) == 0]

m_fam <- read.table(paste0("/juice_data/zhh_independent/ZHH_update/Genotype/",tis,"/",tis,"_male.autosomal.geno.fam"))
otdata_m<-data.frame(FID=0,IID=m_fam$V2,z_matrix_m_no_na[m_fam$V2, ])
Otfile_m=file.path(Otdir,tis,paste0(tis,"_male.expr_tmm_inv.zscore"))
write.table(otdata_m, Otfile_m, quote=F,row.names=F,sep="\t")
}
