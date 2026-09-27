library(car)
library(data.table)
encode_genotype_auto <- function(genotype_vector) {
  all_alleles <- unique(unlist(strsplit(as.character(genotype_vector), " ", fixed = TRUE)))
  
  if (length(all_alleles) < 2) {
    warning("Only one unique allele found. All genotypes will be coded as 0.")
    return(rep(0, length(genotype_vector)))
  }
  
  sorted_alleles <- sort(all_alleles)
  allele_to_count <- sorted_alleles[2]
  
  numeric_genotypes <- sapply(strsplit(as.character(genotype_vector), " ", fixed = TRUE), 
                               function(x) sum(x == allele_to_count))
  
  return(numeric_genotypes)
}

run_interaction_model <- function(lead_df, pheno, Meta, geno_df, trait_name) {
  pval_res <- data.frame()
  for (i in seq(1, nrow(lead_df))) {
    vGene <- lead_df$Gene[i]
    SNP <- lead_df$SNP[i]

    tmp_df <- data.frame(IID = substr(pheno$IID, 3, nchar(pheno$IID)), expr = pheno[[vGene]])
    tmp_df[[trait_name]] <- Meta[[trait_name]][match(tmp_df$IID, Meta$IID)]
    #genotype_list <- geno_df[[SNP]]
    #tmp_df$genotype <- encode_genotype_auto(genotype_list)

    matching_cols <- colnames(geno_df)[grepl(SNP, colnames(geno_df))]
    tmp_df$genotype <- geno_df[[matching_cols]]

    form <- as.formula(paste("expr ~", trait_name, "* genotype"))
    model <- lm(form, data = tmp_df)

    coef_matrix <- summary(model)$coefficients
    interaction_name <- paste0(trait_name, ":genotype")

    if (interaction_name %in% rownames(coef_matrix)) {
      interaction_p_value <- coef_matrix[interaction_name, "Pr(>|t|)"]
    } else {
      interaction_p_value <- NA
    }
    #interaction_p_value <- coef_matrix[4, "Pr(>|t|)"]
    form2 <- as.formula(paste("expr ~", trait_name, "+ genotype"))
    model_no_interaction <- lm(form2, data = tmp_df)
    anova_res <- anova(model_no_interaction, model)
    anova_pval <- anova_res$"Pr(>F)"[2]
    tmp_res <- data.frame(gene = vGene, pval_inte = interaction_p_value, pval_anova = anova_pval)
    pval_res <- rbind(pval_res, tmp_res)
  }
  pval_res$interaction <- trait_name
  return(pval_res)
}

#Meta<-read.table("/juice_data/zhh_independent/ZHH_update/vQTL/Meta_info.txt",head=T)
#traits <- c("Thymus_weight", "Spleen_weight", "Bursa_weight", "Adipose_weight")

tissue_list<-readLines("/soft/zhh/vQTL/OSCA/tissue.list")
tissue<-tissue_list[1]
for (tissue in tissue_list){
print(paste0("Processing on ",tissue))
#phenofile=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/pheno/",tissue,"/",tissue,".expr_tmm_inv.zscore")
#pheno<-fread(phenofile)
Indir1<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore_update")
#files1<-list.files(Indir1,pattern="\\.query$",full.names=T)
lead_df <- read.table(file.path(Indir1,paste0(tissue,"_Vgene.lead.txt")),head=T)
lead_vqtl <- unique(lead_df$SNP)
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/GXE/")
#dir.create(Otdir)
vcf=paste0("/ossfs/ZHH_update/ZHH_update/Genotype/",tissue,"/",tissue,".autosomal.geno.vcf.gz")
geno_prefix=paste0("/ossfs/ZHH_update/ZHH_update/Genotype/",tissue,"/",tissue,".autosomal.geno")
snp_list=paste0(Otdir,tissue,"_lead.vQTL.list")
write.table(lead_vqtl,snp_list,quote=F,col.names=F,row.names=F)
output_prefix<-gsub(".list","",snp_list)
geno_file<-gsub(".list",".genotype",snp_list)
plink_cmd <- paste0("plink --vcf ", vcf, " --chr-set 39 --keep-allele-order --extract ", snp_list, " --recode A --out ",geno_file)
system(plink_cmd)
}

#plink_command <- paste("plink --bfile", geno_prefix,
#                       "--extract", snp_list,
#                       "--r2 --chr-set 39 --keep-allele-order",
#                       "--ld-snp-list", snp_list,
#                       "--ld-window-kb", "5000",
#                       "--ld-window-r2", "0",
#                       "--out", output_prefix)

#system(plink_command)

#ld_matrix_df <- read.table(paste0(output_prefix, ".ld"), header = TRUE)
#ld_matrix_df1 <- ld_matrix_df[ld_matrix_df$R2!=1, ]
#ld_matrix_df2 <- ld_matrix_df1[ld_matrix_df1$R2>0.01, ]
#SNP_remove<-unique(c(ld_matrix_df2$SNP_A, ld_matrix_df2$SNP_B))
#SNP_keep<-setdiff(unique(lead_snps$SNP), SNP_remove)
#ld_matrix <- ld_matrix_df1[ld_matrix_df1$SNP_A %in% SNP_keep | ld_matrix_df1$SNP_B %in% SNP_keep, ]

Slaughter_df <- read.table("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Slaughter_trait.txt",sep="\t",head=T)
interactions <- readLines("/soft/zhh/vQTL/GXE/slaugghter_trait.1.list")

all_inte<-data.frame()
for (tissue in tissue_list){
print(paste0("Processing on ",tissue))
phenofile=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/pheno/",tissue,"/",tissue,".expr_tmm_inv.zscore")
pheno<-fread(phenofile)
Indir1<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore_update")
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/GXE/")
snp_list=paste0(Otdir,tissue,"_lead.vQTL.list")
geno_file<-gsub(".list",".genotype",snp_list)
geno_df <- fread(paste0(geno_file,".raw"))

lead_df <- read.table(file.path(Indir1,paste0(tissue,"_Vgene.lead.txt")),head=T)
#lead_snps1 <- lead_snps[lead_snps$SNP %in% SNP_keep, ]

  all_pval_list <- lapply(interactions, function(interactions) {
    run_interaction_model(lead_df, pheno, Slaughter_df, geno_df, interactions)
  })

all_pval <- do.call(rbind, all_pval_list)
all_pval$tissue <- tissue
all_pval$FDR<-p.adjust(all_pval$pval_anova,method="BH")
#all_pval1<-all_pval[all_pval$FDR<0.05,]
threshold <- 0.05/(nrow(lead_df)*4)
all_pval1<-all_pval[all_pval$pval_inte<threshold | all_pval$pval_anova<threshold,]
Output2=paste0(Otdir,tissue,"_slaughter.GXE.all.res.txt")
write.table(all_pval, Output2, quote=F,row.names=F,sep="\t")

if (nrow(all_pval1) > 0) { 
Output1=paste0(Otdir,tissue,"_slaughtern.GXE.res.txt")
write.table(all_pval1, Output1, quote=F,row.names=F,sep="\t")
all_inte<-rbind(all_pval1, all_inte)
} else {
}
}
write.table(all_inte, "/juice_data/zhh_independent/ZHH_update/vQTL/res/vQTL_slaughter.GXE.res.txt",quote=F,row.names=F,sep="\t")

###Growth traits intersection###
Growth_df <- read.table("/juice_data/zhh_independent/ZHH_update/GxE/veQTL/interaction_file/Growth_trait.txt",sep="\t",head=T)
interactions <- names(Growth_df)[-c(1:2)]
all_Growth<-data.frame()
for (tissue in tissue_list){
print(paste0("Processing on ",tissue))
phenofile=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/pheno/",tissue,"/",tissue,".expr_tmm_inv.zscore")
pheno<-fread(phenofile)
Indir1<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore_update")
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/GXE/")
snp_list=paste0(Otdir,tissue,"_lead.vQTL.list")
geno_file<-gsub(".list",".genotype",snp_list)
geno_df <- fread(paste0(geno_file,".raw"))

lead_df <- read.table(file.path(Indir1,paste0(tissue,"_Vgene.lead.txt")),head=T)
#lead_snps1 <- lead_snps[lead_snps$SNP %in% SNP_keep, ]

  all_pval_list <- lapply(interactions, function(interactions) {
    run_interaction_model(lead_df, pheno, Growth_df, geno_df, interactions)
  })

all_pval <- do.call(rbind, all_pval_list)
all_pval$tissue <- tissue
all_pval$FDR<-p.adjust(all_pval$pval_anova,method="BH")
#all_pval1<-all_pval[all_pval$FDR<0.05,]
threshold <- 0.05/(nrow(lead_df)*4)
all_pval1<-all_pval[all_pval$pval_inte<threshold | all_pval$pval_anova<threshold,]
Output2=paste0(Otdir,tissue,"_Growth.GXE.all.res.txt")
write.table(all_pval, Output2, quote=F,row.names=F,sep="\t")

if (nrow(all_pval1) > 0) {
Output1=paste0(Otdir,tissue,"_Growth.GXE.res.txt")
write.table(all_pval1, Output1, quote=F,row.names=F,sep="\t")
all_Growth<-rbind(all_pval1, all_Growth)
} else {
}
}
write.table(all_Growth, "/juice_data/zhh_independent/ZHH_update/vQTL/res/vQTL_Growth.GXE.res.txt",quote=F,row.names=F,sep="\t")

###Cell proportion intersection###

all_inte_cell<-data.frame()
for (tissue in tissue_list[-1]){
print(paste0("Processing on ",tissue))
phenofile=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/pheno/",tissue,"/",tissue,".expr_tmm_inv.zscore")
pheno<-fread(phenofile)
Indir1<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore_update")
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/GXE/")
snp_list=paste0(Otdir,tissue,"_lead.vQTL.list")
geno_file<-gsub(".list",".genotype",snp_list)
geno_df <- fread(paste0(geno_file,".raw"))

lead_df <- read.table(file.path(Indir1,paste0(tissue,"_Vgene.lead.txt")),head=T)
#lead_snps1 <- lead_snps[lead_snps$SNP %in% SNP_keep, ]

cell_dir <- paste0("/juice_data/zhh_independent/ZHH_update/03_DWLS_Results_Use/DWLS_INT_74_TCPairs_MedianPro0.1/",tissue)
cell_file <- list.files(cell_dir,pattern="_result\\.txt$",full.names=T)
if(length(cell_file)>0){
cell_types <- gsub("_result.txt","",basename(cell_file))
cell_df <- read.table(cell_file[1])
if(length(cell_file)>1){
for(f in cell_file[-1]){
tmp_cell <- read.table(f)
cell_df <- merge(cell_df,tmp_cell,by="V1")
}
}
names(cell_df) <- c("SID",cell_types)
cell_df$IID <- substr(cell_df$SID,3,nchar(cell_df$SID))

  all_pval_list <- lapply(cell_types, function(cell_types) {
    run_interaction_model(lead_df, pheno, cell_df, geno_df, cell_types)
  })

all_pval <- do.call(rbind, all_pval_list)
all_pval$tissue <- tissue
all_pval$FDR<-p.adjust(all_pval$pval_anova,method="BH")
#all_pval1<-all_pval[all_pval$FDR<0.05,]
threshold <- 0.05/(nrow(lead_df)*4)
all_pval1<-all_pval[all_pval$pval_inte<threshold | all_pval$pval_anova<threshold,]
Output2=paste0(Otdir,tissue,"_cell.GXE.all.res.txt")
write.table(all_pval, Output2, quote=F,row.names=F,sep="\t")

if (nrow(all_pval1) > 0) {
Output1=paste0(Otdir,tissue,"_cell.GXE.res.txt")
write.table(all_pval1, Output1, quote=F,row.names=F,sep="\t")
all_inte_cell<-rbind(all_pval1, all_inte_cell)
} else {
}
}else{
}
}
write.table(all_inte_cell, "/juice_data/zhh_independent/ZHH_update/vQTL/res/vQTL_cell.GXE.res.txt",quote=F,row.names=F,sep="\t")

###Etime, Sex intersection analysis###
run_interaction_factor <- function(lead_df, pheno, Meta, geno_df, trait_name) {
  pval_res <- data.frame(gene = character(), 
                          pval_inte = numeric(), 
                          pval_anova = numeric(), 
                          stringsAsFactors = FALSE)
  for (i in seq(1, nrow(lead_df))) {
    vGene <- lead_df$Gene[i]
    SNP <- lead_df$SNP[i]

    tryCatch({

    tmp_df <- data.frame(IID = substr(pheno$IID, 3, nchar(pheno$IID)), expr = pheno[[vGene]])
    tmp_df[[trait_name]] <- Meta[[trait_name]][match(tmp_df$IID, Meta$IID)]
    tmp_df[[trait_name]] <- as.factor(tmp_df[[trait_name]])

    if(nlevels(tmp_df[[trait_name]])<2){next}
    matching_cols <- colnames(geno_df)[grepl(SNP, colnames(geno_df))]
    tmp_df$genotype <- geno_df[[matching_cols]]
    model <- lm(expr ~ factor(get(trait_name)) * genotype, data = tmp_df)
    coef_matrix <- summary(model)$coefficients

    interaction_terms <- grep(":genotype", rownames(coef_matrix))
    interaction_p_values <- coef_matrix[interaction_terms, "Pr(>|t|)"]
    interaction_p_value <- min(interaction_p_values, na.rm = TRUE)

    model_no_interaction <- lm(expr ~ factor(get(trait_name)) + genotype, data = tmp_df)
    anova_res <- anova(model_no_interaction, model)
    anova_pval <- anova_res$"Pr(>F)"[2]
    tmp_res <- data.frame(gene = vGene, pval_inte = interaction_p_value, pval_anova = anova_pval)
    pval_res <- rbind(pval_res, tmp_res)
    }, error = function(e) {
      warning(paste("Error in iteration", i, "Gene:", vGene, ":", e$message))
    })
}
  pval_res$interaction <- rep(trait_name, nrow(pval_res))
  return(pval_res)
}

Factors <- c("SEX", "Etime")
all_inte_factor<-data.frame()
for (tissue in tissue_list){
print(paste0("Processing on ",tissue))
phenofile=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/pheno/",tissue,"/",tissue,".expr_tmm_inv.zscore")
pheno<-fread(phenofile)
Indir1<-paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/zscore_update")
Otdir=paste0("/juice_data/zhh_independent/ZHH_update/vQTL/res/",tissue,"/GXE/")
snp_list=paste0(Otdir,tissue,"_lead.vQTL.list")
geno_file<-gsub(".list",".genotype",snp_list)
geno_df <- fread(paste0(geno_file,".raw"))

lead_df <- read.table(file.path(Indir1,paste0(tissue,"_Vgene.lead.txt")),head=T)
#lead_snps1 <- lead_snps[lead_snps$SNP %in% SNP_keep, ]

  all_pval_list <- lapply(Factors, function(Factors) {
    run_interaction_factor(lead_df, pheno, Slaughter_df, geno_df, Factors)
  })

all_pval <- do.call(rbind, all_pval_list)
all_pval$tissue <- tissue
all_pval$FDR<-p.adjust(all_pval$pval_anova,method="BH")
#all_pval1<-all_pval[all_pval$FDR<0.05,]
threshold <- 0.05/(nrow(lead_df)*4)
all_pval1<-all_pval[all_pval$pval_inte<threshold | all_pval$pval_anova<threshold,]
Output2=paste0(Otdir,tissue,"_factor.GXE.all.res.txt")
write.table(all_pval, Output2, quote=F,row.names=F,sep="\t")

if (nrow(all_pval1) > 0) {
Output1=paste0(Otdir,tissue,"_factor.GXE.res.txt")
write.table(all_pval1, Output1, quote=F,row.names=F,sep="\t")
all_inte_factor<-rbind(all_pval1, all_inte_factor)
} else {
}
}
write.table(all_inte_factor, "/juice_data/zhh_independent/ZHH_update/vQTL/res/vQTL_factor.GXE.res.txt",quote=F,row.names=F,sep="\t")
