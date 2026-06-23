#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(rMVP)
  library(CMplot)
  library(bigmemory)
})

work_dir <- Sys.getenv("HPGWAS_RMVP_DIR", unset = ".")
vcf_file <- Sys.getenv("HPGWAS_VCF", unset = "output.vcf")
phenotype_file <- Sys.getenv("HPGWAS_PHENOTYPE", unset = "phen.txt")
output_prefix <- Sys.getenv("HPGWAS_RMVP_PREFIX", unset = "mvp.vcf")
n_cpus <- as.integer(Sys.getenv("HPGWAS_NCPUS", unset = "10"))

setwd(work_dir)

MVP.Data(
  fileVCF = vcf_file,
  filePhe = phenotype_file,
  fileKin = TRUE,
  filePC = TRUE,
  out = output_prefix
)

genotype <- attach.big.matrix(paste0(output_prefix, ".geno.desc"))
phenotype <- read.table(paste0(output_prefix, ".phe"), header = TRUE)
map <- read.table(paste0(output_prefix, ".geno.map"), header = TRUE)
kinship <- attach.big.matrix(paste0(output_prefix, ".kin.desc"))

gwas_result <- MVP(
  phe = phenotype,
  geno = genotype,
  map = map,
  K = kinship,
  nPC.GLM = 5,
  nPC.MLM = 3,
  nPC.FarmCPU = 3,
  ncpus = n_cpus,
  vc.method = "BRENT",
  maxLoop = 10,
  method.bin = "static",
  threshold = 0.1,
  method = c("GLM", "MLM", "FarmCPU"),
  file.output = TRUE
)

MVP.Report(
  gwas_result,
  plot.type = "q",
  col = c("dodgerblue1", "olivedrab3", "darkgoldenrod1"),
  threshold = 1e10,
  signal.pch = 19,
  signal.cex = 1.5,
  signal.col = "red",
  box = FALSE,
  multracks = TRUE,
  file.type = "pdf",
  memo = "",
  dpi = 300
)

pca <- attach.big.matrix(paste0(output_prefix, ".pc.desc"))[, 1:3]
MVP.PCAplot(
  PCA = pca,
  Ncluster = 2,
  legend.pos = "topright",
  class = NULL,
  col = c("red", "green"),
  file.type = "jpg"
)

