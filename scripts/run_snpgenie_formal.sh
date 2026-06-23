#!/usr/bin/env bash
set -euo pipefail

base="${1:-${HPGWAS_DNDS_DIR:-$PWD}}"
script_dir="${2:-${SNPGENIE_DIR:-$base/SNPGenie-master}}"
analysis="$base/formal_analysis"
plus="$analysis/plus"
minus="$analysis/minus"
seq_len="$(awk 'BEGIN{n=0} /^>/{next} {gsub(/[[:space:]]/,""); n+=length($0)} END{print n}' "$base/sequence.fasta")"

mkdir -p "$plus" "$minus" "$analysis/logs"

echo "[START] $(date '+%F %T')"
echo "[INFO] base=$base"
echo "[INFO] SNPGenie directory=$script_dir"
echo "[INFO] sequence_length=$seq_len"

python3 "$(dirname "$0")/prepare_snpgenie_inputs.py" --base "$base" --analysis "$analysis"

echo "[INFO] Running SNPGenie plus-strand analysis"
perl "$script_dir/snpgenie.pl" \
  --vcfformat=1 \
  --snpreport="$plus/core_snpgenie_format1.vcf" \
  --fastafile="$plus/sequence.fasta" \
  --gtffile="$plus/genomic_snpgenie_clean.gtf" \
  --workdir="$plus" \
  --outdir="$analysis/SNPGenie_Results_plus_clean" \
  > "$analysis/logs/snpgenie_plus.stdout.log" 2> "$analysis/logs/snpgenie_plus.stderr.log"

echo "[INFO] Preparing reverse-complement inputs for minus-strand analysis"
cp "$plus/core_snpgenie_format1.vcf" "$minus/core_snpgenie_format1.vcf"
cp "$plus/genomic_snpgenie_clean.gtf" "$minus/genomic_snpgenie_clean.gtf"
cp "$plus/sequence.fasta" "$minus/sequence.fasta"

(
  cd "$minus"
  rm -f core_snpgenie_format1_revcom.vcf genomic_snpgenie_clean_revcom.gtf sequence_revcom.fasta
  perl "$script_dir/vcf2revcom.pl" core_snpgenie_format1.vcf "$seq_len" > "$analysis/logs/vcf2revcom.stdout.log" 2> "$analysis/logs/vcf2revcom.stderr.log"
  perl "$script_dir/gtf2revcom.pl" genomic_snpgenie_clean.gtf "$seq_len" > "$analysis/logs/gtf2revcom.stdout.log" 2> "$analysis/logs/gtf2revcom.stderr.log"
  perl "$script_dir/fasta2revcom.pl" sequence.fasta > "$analysis/logs/fasta2revcom.stdout.log" 2> "$analysis/logs/fasta2revcom.stderr.log"
)

echo "[INFO] Running SNPGenie minus-strand analysis"
perl "$script_dir/snpgenie.pl" \
  --vcfformat=1 \
  --snpreport="$minus/core_snpgenie_format1_revcom.vcf" \
  --fastafile="$minus/sequence_revcom.fasta" \
  --gtffile="$minus/genomic_snpgenie_clean_revcom.gtf" \
  --workdir="$minus" \
  --outdir="$analysis/SNPGenie_Results_minus_clean" \
  > "$analysis/logs/snpgenie_minus.stdout.log" 2> "$analysis/logs/snpgenie_minus.stderr.log"

echo "[SUMMARY] plus_dir=$analysis/SNPGenie_Results_plus_clean"
echo "[SUMMARY] minus_dir=$analysis/SNPGenie_Results_minus_clean"
echo "[END] $(date '+%F %T')"

