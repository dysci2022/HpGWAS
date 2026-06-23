#!/usr/bin/env bash
set -euo pipefail

snpeff_dir="${1:-${SNPEFF_DIR:-$PWD/snpEff}}"
genome_id="${HPGWAS_SNPEFF_GENOME:-hp_26695}"
genome_name="${HPGWAS_SNPEFF_GENOME_NAME:-Helicobacter pylori 26695}"
genbank_file="${2:-${HPGWAS_GENBANK:-sequence.gb}}"
snplist="${3:-${HPGWAS_SNPLIST:-snplist.txt}}"

cd "$snpeff_dir"

grep -qxF "$genome_id.genome: $genome_name" snpEff.config ||
  echo "$genome_id.genome: $genome_name" >> snpEff.config

mkdir -p "data/$genome_id"
cp "$genbank_file" "data/$genome_id/genes.gbk"

java -Xmx4g -jar snpEff.jar build -genbank -v "$genome_id"

{
  printf '##fileformat=VCFv4.2\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
  awk '{print $1"\t"$2"\t.\t"$3"\t"$4"\t.\t.\t."}' "$snplist"
} > snp.vcf

java -Xmx4g -jar snpEff.jar -v "$genome_id" snp.vcf > annotated.vcf

