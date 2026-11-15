#!/usr/bin/env python3

import csv
import gzip
import re
import sys
from pathlib import Path

if len(sys.argv) != 5:
    raise SystemExit(
        "Uso: 16_vcf_to_iupac_fasta.py "
        "entrada.vcf.gz saida.fasta resumo.tsv missing_por_amostra.tsv"
    )

vcf_path = Path(sys.argv[1])
fasta_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
sample_qc_path = Path(sys.argv[4])

iupac = {
    frozenset(["A"]): "A",
    frozenset(["C"]): "C",
    frozenset(["G"]): "G",
    frozenset(["T"]): "T",
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
}

states = {
    "A": set("A"),
    "C": set("C"),
    "G": set("G"),
    "T": set("T"),
    "R": set("AG"),
    "Y": set("CT"),
    "S": set("GC"),
    "W": set("AT"),
    "K": set("GT"),
    "M": set("AC"),
}

opener = gzip.open if str(vcf_path).endswith(".gz") else open

samples = []
sequences = []

input_records = 0
biallelic_snp_records = 0
retained_sites = 0
constant_compatible_removed = 0
non_biallelic_or_non_snp_removed = 0
missing_gt = []

with opener(vcf_path, "rt") as handle:
    for line in handle:
        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            samples = header[9:]
            sequences = [bytearray() for _ in samples]
            missing_gt = [0 for _ in samples]
            continue

        if line.startswith("#"):
            continue

        input_records += 1
        fields = line.rstrip("\n").split("\t")

        ref = fields[3].upper()
        alt_values = fields[4].upper().split(",")

        if (
            len(ref) != 1
            or ref not in "ACGT"
            or len(alt_values) != 1
            or len(alt_values[0]) != 1
            or alt_values[0] not in "ACGT"
        ):
            non_biallelic_or_non_snp_removed += 1
            continue

        alt = alt_values[0]
        biallelic_snp_records += 1

        format_keys = fields[8].split(":")

        try:
            gt_index = format_keys.index("GT")
        except ValueError:
            raise SystemExit(
                f"ERRO: campo GT ausente em {fields[0]}:{fields[1]}"
            )

        alleles = [ref, alt]
        site_codes = []
        compatible_intersection = None

        for sample_field in fields[9:]:
            values = sample_field.split(":")

            if gt_index >= len(values):
                code = "N"
            else:
                genotype = values[gt_index]

                if "." in genotype:
                    code = "N"
                else:
                    genotype_alleles = re.split(r"[/|]", genotype)

                    try:
                        bases = {
                            alleles[int(index)]
                            for index in genotype_alleles
                        }
                        code = iupac.get(frozenset(bases), "N")
                    except (ValueError, IndexError):
                        code = "N"

            site_codes.append(code)

            if code != "N":
                current_states = states[code]

                if compatible_intersection is None:
                    compatible_intersection = set(current_states)
                else:
                    compatible_intersection &= current_states

        variable_for_iqtree = (
            compatible_intersection is not None
            and len(compatible_intersection) == 0
        )

        if not variable_for_iqtree:
            constant_compatible_removed += 1
            continue

        for index, code in enumerate(site_codes):
            sequences[index].append(ord(code))

            if code == "N":
                missing_gt[index] += 1

        retained_sites += 1

if not samples:
    raise SystemExit("ERRO: nenhuma amostra encontrada no VCF")

if retained_sites == 0:
    raise SystemExit("ERRO: nenhum sítio variável foi retido")

fasta_path.parent.mkdir(parents=True, exist_ok=True)

with fasta_path.open("w") as output:
    for sample, sequence in zip(samples, sequences):
        dna = bytes(sequence).decode("ascii")
        output.write(f">{sample}\n")

        for start in range(0, len(dna), 80):
            output.write(dna[start:start + 80] + "\n")

with summary_path.open("w", newline="") as output:
    writer = csv.writer(
        output,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow(["Metric", "Value"])
    writer.writerow(["Samples", len(samples)])
    writer.writerow(["Input_VCF_records", input_records])
    writer.writerow(["Biallelic_SNP_records", biallelic_snp_records])
    writer.writerow(["Retained_alignment_sites", retained_sites])
    writer.writerow([
        "Constant_compatible_sites_removed",
        constant_compatible_removed
    ])
    writer.writerow([
        "Non_biallelic_or_non_SNP_removed",
        non_biallelic_or_non_snp_removed
    ])

with sample_qc_path.open("w", newline="") as output:
    writer = csv.writer(
        output,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writerow([
        "Analysis_ID",
        "Alignment_sites",
        "Missing_N",
        "Missing_N_pct"
    ])

    for sample, number_missing in zip(samples, missing_gt):
        missing_pct = 100 * number_missing / retained_sites

        writer.writerow([
            sample,
            retained_sites,
            number_missing,
            f"{missing_pct:.4f}"
        ])

print(f"Amostras: {len(samples)}")
print(f"SNPs no VCF: {input_records}")
print(f"Sítios no alinhamento: {retained_sites}")
print(
    "Sítios constantes/compatíveis removidos: "
    f"{constant_compatible_removed}"
)
print(f"FASTA: {fasta_path}")
print(f"Resumo: {summary_path}")
print(f"Missing por amostra: {sample_qc_path}")
