args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript genesis.R <variable_type> <population>\nExample: Rscript genesis.R gaussian Score")
}

variable_type <- args[1]
population <- args[2]
  
if (variable_type == "continuous") {
  family_variable <- "gaussian"
  test_variable <- "Score"
} else if (variable_type == "categorical") {
  family_variable <- "binomial"
  test_variable <- "Score.SPA"
} 

if (population == "AFR") { dir <- "/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/AFR" 
  genotype_file <- paste0(dir, "/AFR.genotype.gds")
  kinship_file <- paste0(dir, "/pcrelate_kinship_AFR.RData")
} else if (population == "EAS") {
  dir <- "/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/EAS" 
  genotype_file <- paste0(dir, "/EAS.genotype.gds")
  kinship_file <- paste0(dir, "/pcrelate_kinship_EAS.RData")
} else if (population == "SAS") {
  dir <- "/mnt/storage/private/mrcieu/research/UKBIOBANK_GWAS_Pipeline/data/SAS" 
  genotype_file <- paste0(dir, "/SAS.genotype.gds")
  kinship_file <- paste0(dir, "/pcrelate_kinship_SAS.RData")
}

num_pcs <- as.integer(args[1]) 
pc_columns <- paste0("PC", 1:num_pcs)

#Load packages

library(GENESIS)
library(GWASTools)
library(gdsfmt)
library(SNPRelate)
library(data.table)
library(dplyr)
library(SeqVarTools)
library(qqman)

#Read in input data
phen <- file.path(input_path, "data", "phenotypes", user, "input", pheno_file)
cov <- file.path(input_path, "data", "phenotypes", user, "input", cov_file) 
pc <- file.path(input_path, "data", "phenotypes", user, "input", cov_file)

#Join phenotype, covariate and PC data files and prepare input data
cov_select <- cov %>% select(IID, age, sex) 
cov_replace <- cov_select %>% mutate(sex = ifelse(sex == 1,"M","F"))
phen_select_con <- phen_con %>% select("IID", "dbp") #change phenotype
dat <- merge(phen_select_con, cov_replace, by="IID") %>% merge(., pc, by="IID") 
dat <- dat %>% na.omit()

input_dat <- data.frame(scanID = dat$ieu, pheno = dat$dbp, age = 
                        dat$age, sex = dat$sex,
                        dat[, pc_columns, drop = FALSE])

scanAnnot <- ScanAnnotationDataFrame(input_dat)
scanAnnot

### Make a genetic relationship matrix (GRM) from kinship estimates calculated from PC-Relate analysis
load(file = kinship_file)
kinship_scaled <- pcrelateToMatrix(mypcrelate, thresh = 2^(-11/2), scaleKin = 2)

#Fit null model 
nullmod <- fitNullModel(scanAnnot, outcome = "pheno", covars = c("age", "sex", pc_columns), family = family_variable)

# create a genotype object
geno <- GdsGenotypeReader(filename = genotype_file)
genoData <- GenotypeData(geno)

#Run SNP-phenotype association test
genoIterator <- GenotypeBlockIterator(genoData)
assoc <- assocTestSingle(genoIterator, null.model = nullmod, test= test_variable, BPPARAM = BiocParallel::SerialParam())
ea <- effectAllele(genoData, variant.id=assoc$variant.id)
assoc$chr <- as.numeric(assoc$chr)
assoc_out <- merge(assoc, ea, by="variant.id")
assoc_out_sorted <- assoc_out[order(assoc_out$chr),]

#Save output files 
output_file <- file.path(input_path, "data", "phenotypes", user, "output", pheno_name, uid, paste0(pheno_name, "_out.txt.gz")) 
write.table(assoc_out_sorted, gzfile(output_file), sep = "\t", row.names = FALSE, quote = FALSE)

#Create Manhattan and Q-Q plot 
png(file= "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/AFR/GENESIS_dbp_assoc_manhattan.png", width=8, height=7, units="in", res=300)
manhattan(assoc_out_sorted, chr="chr", bp="pos", snp="variant.id", p="Score.pval" )
dev.off()

png(file= "/user/work/ac14629/MRC_network_project/results/UKB/GWAS/GENESIS/AFR/GENESIS_dbp_assoc_qq.png", width=8, height=7, units="in", res=300)
qq(assoc_out_sorted$Score.pval)
dev.off()
