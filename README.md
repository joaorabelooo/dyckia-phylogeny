# Dyckia phylogeny: custom analysis scripts

> **AI disclosure:** This repository was generated and organized with assistance from OpenAI Codex. The custom scripts and documentation were developed with AI assistance during the research workflow.

Custom scripts supporting GBS-based phylogenetic analyses of *Dyckia* and sampled *Encholirium* lineages, with *Deuterocohnia meziana* as the outgroup.

## Scope and publication status

This repository documents analysis from an existing genotype-masked VCF through SNP filtering, alignment preparation, ML and SVDquartets inference, and comparison of the final trees. It includes the complete converter and executors retrieved from the cluster, the SNP-spacing script, and the final-tree audit. The original project paths are preserved, so paths must be adapted to another installation. Input VCFs, sample lists, third-party software, and inferred trees must be supplied separately.

Read mapping and variant calling are outside the scope of this repository. Raw reads and genotypes are not included in this deposit. Sampling metadata are available in Table S1.

## Supplementary materials

- [Figure S1: SVDquartets 10-kb tree, SVG panels A–E](supplementary/README.md#figure-s1).
- [Table S1: sampling information for the 202 accessions (Excel)](supplementary/Table_S1_Dyckia_sampling.xlsx).

See the [supplementary index and figure legend](supplementary/README.md) for details.

## Available scripts

### `18_random_min_distance_sites.py`

Selects SNP positions subject to a minimum physical distance within each reference pseudochromosome. The deposited code preserves the selection logic of the script inspected on the cluster.

Requires Python 3; uses only the standard library (`bisect`, `random`, `sys`, and `collections`).

Input: a headerless, tab-separated position list with chromosome in column 1 and one-based position in column 2. Additional columns are ignored. Output: selected chromosome-position pairs, ordered by coordinate within each chromosome; chromosomes retain their first-appearance order in the input.

```bash
bcftools query -f '%CHROM\t%POS\n' filtered.vcf.gz > all.positions.tsv

python 18_random_min_distance_sites.py all.positions.tsv 100 20261021 panel100.positions.tsv
python 18_random_min_distance_sites.py all.positions.tsv 10000 20261020 panel10000.positions.tsv

bcftools view -T panel100.positions.tsv -Oz -o panel100.vcf.gz filtered.vcf.gz
bcftools index panel100.vcf.gz
```

For each chromosome, candidate positions are shuffled using a seeded Python random-number generator. A candidate is accepted only if its distance from its nearest retained neighbors on both sides is at least the specified threshold. Each panel is selected separately from the complete input position list. This greedy procedure does not guarantee the largest possible number of retained sites or uniform sampling across all possible valid subsets.

The same ordered input and seed should be retained for reproducibility, together with the Python version. The input SNP list and filtered VCF are external inputs, not included here. In the final 202-accession analysis, the 100-bp and 10-kb panels contained 8,291 and 4,934 SNPs, respectively. Physical spacing does not explicitly delimit GBS loci or directly estimate linkage disequilibrium.

### `audit_final202.py`

Audits ML and SVDquartets trees for the final 202-accession dataset. Uses only the Python standard library. Checks terminal sets and exclusions, extracts nontrivial unrooted bipartitions, summarizes support, and compares shared and incompatible bipartitions.

The original script contains project-specific paths in `PROJECT` and `RUN`. Edit these constants for another installation. The `main()` function defines the expected tree and metadata paths. The expected outgroup label and six excluded accession IDs are also specified in the script; historical analysis IDs are retained for traceability.

```bash
python audit_final202.py
```

Outputs are written to a new timestamped audit directory. The script can generate tables for verification even when those tables are not included as manuscript supplements.

Support criteria used for comparisons:

- ML: SH-aLRT >= 80 AND UFBoot >= 95.
- SVDquartets: standard bootstrap >= 70.
- Shared bipartitions: identical accession composition after orienting splits away from the outgroup.
- Supported conflict: incompatible bipartitions that each meet the support criteria of their respective analysis.

Absence of an exact split is not treated as bootstrap zero. Different support statistics remain distinct measurements.

## Converter and inference executors

| File | Purpose and provenance |
| --- | --- |
| `16_vcf_to_iupac_fasta.py` | Original cluster copy. Converts biallelic SNP genotypes into IUPAC-coded FASTA, removes constant-compatible sites, and writes alignment summaries and missingness per accession. |
| `17_run_main208_ml.sh` | Original cluster copy. Template used by the final executor to generate its 202-accession ML launcher; not a separate final analysis. |
| `25_run_svdq_bootstrap100.sh` | Previously inspected original SVD template. Applying the substitutions in `run_final202.sh` reproduces the downloaded final SVD executor, apart from trailing whitespace. |
| `run_final202.sh` | Original cluster copy. Validates the final input, prepares ML and SVD matrices, adapts the two templates, and runs all three analyses. |
| `run_svd.sh` | Original generated executor downloaded from the final cluster run; included to document the actual PAUP commands. |
| `PROVENANCE.json` | Source records and SHA256 hashes of the deposited analysis scripts. |

Python scripts use the standard library. Execution of the inference workflow additionally requires Bash, BCFtools (including the `fill-tags` plugin for filtering), IQ-TREE, PAUP*, vcf2phylip, GNU time at `/usr/bin/time`, and standard Unix utilities. The repository does not redistribute third-party tools.

Example converter invocation:

```bash
python 16_vcf_to_iupac_fasta.py filtered.vcf.gz alignment.fasta alignment_summary.tsv alignment_missingness.tsv
```

The converter retains a site only if the intersection of nucleotide states compatible with all nonmissing IUPAC calls is empty. Missing genotypes are encoded as `N`. This implements removal of sites compatible with a constant nucleotide state before using ascertainment-bias correction.

The ML template uses `-m MFP+ASC -B 1000 --alrt 1000 --bnni -T 16 --seed 20261007`. The final data selected GTR+F+ASC+R5. The SVD launcher evaluates all quartets, uses QFM, distributes ambiguous states across compatible nucleotides, and runs 100 standard bootstrap replicates with 16 threads. Both SNP selection and PAUP use seed 20261021 for 100 bp and 20261020 for 10 kb.

## Final VCF filtering protocol

The following commands document the validated filtering procedure. They are not variant-calling commands and start with the existing biallelic, genotype-masked VCF. Genotype masking has already been performed on that input.

```bash
INPUT_VCF=dyckia_210_samples.biallelic_snps.gtmasked.vcf.gz
KEEP=samples_final_202.txt
OUTPUT_VCF=dyckia_final202.qual30.miss10.dpmax17888.mq40.mac2.vcf.gz

set -euo pipefail
bcftools view -S "$KEEP" -Ou "$INPUT_VCF" |
    bcftools +fill-tags -Ou -- -t AC,AN,AF,MAF,F_MISSING |
    bcftools view -i 'QUAL>=30 && INFO/F_MISSING<=0.10 && INFO/DP<=17888 && INFO/MQ>=40' -Ou |
    bcftools view -c 2:minor -Oz -o "$OUTPUT_VCF"
bcftools index "$OUTPUT_VCF"
```

Allele-count and missingness tags are recalculated after selecting accessions. Site-level `INFO/DP` and `INFO/MQ` are inherited from the existing input; this procedure does not recalculate them from alignments. The documented final dataset contained 202 accessions and 95,059 variants.

## Running the original workflow

1. Supply the filtered VCF and the final accession list. The list must match the sample names in the VCF.
2. Place `16_vcf_to_iupac_fasta.py`, `17_run_main208_ml.sh`, `18_random_min_distance_sites.py`, and `25_run_svdq_bootstrap100.sh` in the project's `scripts/` directory.
3. Edit `PROJECT` in the executors and the input `VCF`, `KEEP`, and third-party converter path in `run_final202.sh`. The original 202-sample and 95,059-variant assertions are intentionally retained.
4. Create a run directory and execute the launcher, for example:

```bash
mkdir -p /path/to/project/reproduced_final202
bash run_final202.sh /path/to/project/reproduced_final202
```

5. Edit `PROJECT` and `RUN` in `audit_final202.py` to match the project and resulting run, supply its expected metadata files, and run the audit.

`run_final202.sh` generates its own `run_ml.sh` and `run_svd.sh` from the templates. The separately deposited `run_svd.sh` records the historical final run and contains that run's absolute directory. The scripts are archival research executors, not a path-independent package.

