#!/usr/bin/env python3
from pathlib import Path
import argparse


def normalize_gt(gt):
    gt = gt.split(":", 1)[0]
    if gt in {".", "./.", ".|."}:
        return []
    alleles = []
    for part in gt.replace("|", "/").split("/"):
        if part in {".", ""}:
            continue
        try:
            alleles.append(int(part))
        except ValueError:
            continue
    return alleles


def get_gene_id(attrs):
    for part in attrs.split(";"):
        part = part.strip()
        if part.startswith("gene_id"):
            pieces = part.split(" ", 1)
            if len(pieces) == 2:
                return pieces[1].strip().strip('"')
    return None


def make_clean_gtf(in_gtf, out_gtf, excluded_path):
    cds_lengths = {}
    records = []
    with in_gtf.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                records.append((line, None, None))
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                records.append((line, None, None))
                continue
            gene_id = get_gene_id(fields[8])
            feature = fields[2]
            if feature == "CDS" and gene_id:
                cds_lengths[gene_id] = cds_lengths.get(gene_id, 0) + int(fields[4]) - int(fields[3]) + 1
            records.append((line, feature, gene_id))

    bad = {gene for gene, length in cds_lengths.items() if length % 3 != 0}
    with out_gtf.open("w") as out:
        for line, feature, gene_id in records:
            if feature == "CDS" and gene_id in bad:
                continue
            if feature == "CDS":
                out.write(line)
    with excluded_path.open("w") as out:
        out.write("gene_id\tcds_length\tlength_mod_3\n")
        for gene in sorted(bad):
            out.write(f"{gene}\t{cds_lengths[gene]}\t{cds_lengths[gene] % 3}\n")


def convert_vcf(in_vcf, out_vcf):
    converted = 0
    skipped_non_snp = 0
    with in_vcf.open() as inp, out_vcf.open("w") as out:
        for line in inp:
            if line.startswith("##"):
                out.write(line)
                if line.startswith("##INFO=<ID=TYPE"):
                    out.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Number of called alleles calculated from GT fields for SNPGenie format 1">\n')
                    out.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency calculated from GT fields for SNPGenie format 1">\n')
                continue
            if line.startswith("#CHROM"):
                out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = fields[:8]
            alts = alt.split(",")
            if len(ref) != 1 or any(len(a) != 1 for a in alts):
                skipped_non_snp += 1
                continue
            if ref.upper() not in "ACGT" or any(a.upper() not in "ACGT" for a in alts):
                skipped_non_snp += 1
                continue

            alt_counts = [0] * len(alts)
            an = 0
            for sample_gt in fields[9:]:
                alleles = normalize_gt(sample_gt)
                for allele in alleles:
                    if allele == 0:
                        an += 1
                    elif 1 <= allele <= len(alts):
                        an += 1
                        alt_counts[allele - 1] += 1
            if an == 0:
                continue
            af = ",".join(f"{count / an:.10g}" for count in alt_counts)
            info_parts = [p for p in info.split(";") if p and not p.startswith("AN=") and not p.startswith("AF=")]
            info_parts.extend([f"AN={an}", f"AF={af}"])
            out.write("\t".join([chrom, pos, vid, ref, alt, qual, filt, ";".join(info_parts)]) + "\n")
            converted += 1
    return converted, skipped_non_snp


def main():
    parser = argparse.ArgumentParser(description="Prepare core VCF/FASTA/GTF inputs for SNPGenie format 1.")
    parser.add_argument("--base", default=".", help="Base directory containing core.vcf, sequence.fasta and genomic.gtf.")
    parser.add_argument("--analysis", default=None, help="Output analysis directory. Default: <base>/formal_analysis.")
    args = parser.parse_args()

    base = Path(args.base).resolve()
    analysis = Path(args.analysis).resolve() if args.analysis else base / "formal_analysis"
    plus = analysis / "plus"
    plus.mkdir(parents=True, exist_ok=True)

    in_vcf = base / "core.vcf"
    in_fasta = base / "sequence.fasta"
    in_gtf = base / "genomic.gtf"
    out_vcf = plus / "core_snpgenie_format1.vcf"
    out_fasta = plus / "sequence.fasta"
    out_gtf = plus / "genomic.gtf"
    out_clean_gtf = plus / "genomic_snpgenie_clean.gtf"
    excluded_cds = analysis / "excluded_incomplete_cds_gene_ids.txt"

    out_fasta.write_text(in_fasta.read_text())
    out_gtf.write_text(in_gtf.read_text())
    make_clean_gtf(in_gtf, out_clean_gtf, excluded_cds)
    converted, skipped_non_snp = convert_vcf(in_vcf, out_vcf)

    summary = analysis / "input_conversion_summary.txt"
    summary.write_text(
        "\n".join([
            "SNPGenie input conversion",
            f"source_vcf={in_vcf}",
            f"output_vcf={out_vcf}",
            "method=ALT allele frequencies calculated from per-sample GT fields",
            "vcfformat=1",
            f"converted_snp_records={converted}",
            f"skipped_non_snp_records={skipped_non_snp}",
            f"clean_gtf={out_clean_gtf}",
            f"excluded_incomplete_cds_gene_ids={excluded_cds}",
        ]) + "\n"
    )
    print(summary)


if __name__ == "__main__":
    main()

