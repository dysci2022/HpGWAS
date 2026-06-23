#!/usr/bin/env bash
set -euo pipefail

input_dir="${1:-${HPGWAS_FS_DIR:-$PWD}}"
reference="${2:-${HPGWAS_REFERENCE:-sequence.gb}}"
cpus="${HPGWAS_CPUS:-128}"
parallel_stage2="${HPGWAS_FS_STAGE2_JOBS:-32}"
parallel_stage3="${HPGWAS_FS_STAGE3_JOBS:-64}"
fs_bin="${FINESTRUCTURE_BIN:-fs}"
vcf2cp="${VCF2CP_PL:-vcf2cp.pl}"
makeuniformrecfile="${MAKEUNIFORMRECFILE_PL:-makeuniformrecfile.pl}"

cd "$input_dir"

find "$PWD" -maxdepth 1 -name "*.fna" -type f | sort > file_paths.txt
sed 's#.*/##; s/\.fna$//; s/_//g' file_paths.txt > sample_ids.txt
paste sample_ids.txt file_paths.txt > snippy_samples.tsv

snippy-multi snippy_samples.tsv --ref "$reference" --cpus "$cpus" > runme.sh
bash runme.sh

bcftools query -l core.vcf > old_samples.txt
sed 's/_//g; s/.fna$//' old_samples.txt > new_samples.txt
paste old_samples.txt new_samples.txt > rename_map.txt
bcftools reheader -s rename_map.txt -o core.renamed.vcf core.vcf

bcftools view -m2 -M2 -v snps core.renamed.vcf |
  bcftools view -i 'COUNT(GT="mis")==0' -Oz -o core.filtered.vcf

perl "$vcf2cp" -p 1 -g 10000000 core.filtered.vcf output_project
perl "$makeuniformrecfile" output_project.phase output_project.recombfile

"$fs_bin" project.cp -hpc 1 -idfile output_project.ids -phasefiles output_project.phase -recombfiles output_project.recombfile -ploidy 1 -s3iters 200000 -go

cat project/commandfiles/commandfile1.txt | parallel
"$fs_bin" project.cp -go

cat project/commandfiles/commandfile2.txt | parallel -j "$parallel_stage2"
"$fs_bin" project.cp -go

cat project/commandfiles/commandfile3.txt | parallel -j "$parallel_stage3"
"$fs_bin" project.cp -go

cat project/commandfiles/commandfile4.txt | parallel
"$fs_bin" project.cp -go

"$fs_bin" fs -X -Y -e X2 project_linked.chunkcounts.out project_linked_tree.xml project_linked.mapstate.csv
"$fs_bin" fs -X -Y -e meancoincidence project_linked.chunkcounts.out project_linked_mcmc.xml project_linked.meancoincidence.csv

