#!/usr/bin/env bash
set -euo pipefail

input_dir="${1:-${HPGWAS_FASTA_DIR:-$PWD}}"
reference="${2:-${HPGWAS_REFERENCE:-sequence.gb}}"
cpus="${HPGWAS_CPUS:-128}"
plink_bin="${PLINK_BIN:-plink}"

cd "$input_dir"

find "$PWD" -maxdepth 1 -name "*.fna" -type f | sort > file_paths.txt
sed 's#.*/##; s/\.fna$//' file_paths.txt > sample_ids.txt
paste sample_ids.txt file_paths.txt > snippy_samples.tsv

snippy-multi snippy_samples.tsv --ref "$reference" --cpus "$cpus" > runme.sh
bash runme.sh

# For rRNA-focused mapping, a permissive mapping-quality setting can be used:
# snippy-multi snippy_samples.tsv --ref "$reference" --cpus "$cpus" --mapqual 0 > runme.sh

bcftools view -m2 -M2 -v snps core.vcf -o output.vcf --threads "$cpus"

"$plink_bin" \
  --vcf core.vcf \
  --out output \
  --maf 0.05 \
  --allow-extra-chr \
  --chr-set 1 \
  --recode vcf \
  --indep-pairwise 50 10 0.2 \
  --keep-allele-order

snippy-clean_full_aln core.aln > clean.core.aln
iqtree2 -s clean.core.aln -m GTR+F -T AUTO
iqtree2 -s clean.core.aln -m GTR+F+R5 -T AUTO -fast

