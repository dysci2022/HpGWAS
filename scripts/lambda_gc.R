#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else "SNPlmm_PCA_GRM_allSNPs.txt"

gwas <- fread(input_file, data.table = FALSE)

possible_p_cols <- c("lrt-pvalue", "pvalue", "p-value", "P", "p")
p_col <- possible_p_cols[possible_p_cols %in% colnames(gwas)][1]

if (is.na(p_col)) {
  stop("No p-value column detected. Please check the GWAS output header.")
}

p <- as.numeric(gwas[[p_col]])
p <- p[!is.na(p) & p > 0 & p <= 1]

chisq <- qchisq(1 - p, df = 1)
lambda_gc <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

cat("Input file:", input_file, "\n")
cat("P-value column:", p_col, "\n")
cat("Number of variants used:", length(p), "\n")
cat("Lambda GC:", lambda_gc, "\n")

