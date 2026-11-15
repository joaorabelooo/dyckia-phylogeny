# Dyckia phylogeny: custom analysis scripts

> **AI disclosure:** This repository was generated and organized with assistance from OpenAI Codex. The custom scripts were developed with AI assistance during the research workflow. Inclusion here does not imply independent software validation or independent confirmation of biological conclusions.

Custom scripts supporting GBS-based phylogenetic analyses of *Dyckia* and sampled *Encholirium* lineages, with *Deuterocohnia meziana* as the outgroup.

## Scope and publication status

This initial deposit contains the complete SNP-spacing script and final-tree audit script available locally. Additional analysis scripts will be deposited after the complete versions actually used on the cluster have been retrieved. This is not yet a complete executable pipeline.

Read mapping and variant calling are outside the scope of this repository. Raw reads, genotypes, and sampling metadata are not included in this deposit.

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

## Further scripts awaiting retrieval

- Complete VCF-to-IUPAC FASTA converter used in the analysis.
- Final VCF-subsetting and filtering procedure, beginning with the existing genotype-masked VCF.
- Final ML/SVD executors and the scripts they invoke.

These are not represented by reconstructed or truncated substitutes in this initial deposit.
