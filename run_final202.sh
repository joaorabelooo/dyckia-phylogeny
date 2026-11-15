#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/jrabelo/tese/filogenia_dyckia
RUN=$1
SCRIPTS="$PROJECT/scripts"
VCF="$PROJECT/03_variants/final/final202.GLHLnV/dyckia_final202.qual30.miss10.dpmax17888.mq40.mac2.vcf.gz"
KEEP="$PROJECT/00_metadata/samples_final_202.txt"
CONVERTER="$PROJECT/software/vcf2phylip/vcf2phylip.py"
ML="$RUN/ml"
SVD="$RUN/svdquartets"

trap 'printf "ERRO na linha %s. Confira o log em %s\n" "$LINENO" "$RUN" >&2' ERR

mkdir -p "$ML" "$SVD/input/vcf2phylip"
[[ $(bcftools index -n "$VCF") -eq 95059 ]]
LC_ALL=C sort "$KEEP" > "$RUN/expected.samples.txt"
[[ $(wc -l < "$RUN/expected.samples.txt") -eq 202 ]]
bcftools query -l "$VCF" | LC_ALL=C sort > "$RUN/vcf.samples.txt"
cmp "$RUN/expected.samples.txt" "$RUN/vcf.samples.txt"

# Adapta cópias dos executores originais para este diretório.
sed \
    -e "s|^MLDIR=.*|MLDIR=\"$ML\"|" \
    -e 's/dyckia_main208/dyckia_final202/g' \
    "$SCRIPTS/17_run_main208_ml.sh" > "$RUN/run_ml.sh"

sed \
    -e "s|^SVD=.*|SVD=\"$SVD\"|" \
    -e 's/dyckia_main208/dyckia_final202/g' \
    -e 's/^for DIST in 10000 100$/for DIST in 100 10000/' \
    -e 's/75\.760\.620/67.331.650/g' \
    "$SCRIPTS/25_run_svdq_bootstrap100.sh" > "$RUN/run_svd.sh"

printf '\n### Matriz de ML\n'
FASTA="$ML/dyckia_final202.iupac.variable.fasta"
python "$SCRIPTS/16_vcf_to_iupac_fasta.py" \
    "$VCF" "$FASTA" \
    "$ML/alignment_summary.tsv" "$ML/alignment_missingness.tsv"

awk '/^>/ {sub(/^>/, ""); print $1}' "$FASTA" |
    LC_ALL=C sort > "$ML/fasta.samples.txt"
cmp "$RUN/expected.samples.txt" "$ML/fasta.samples.txt"
cat "$ML/alignment_summary.tsv"

printf '\n### Preparação dos painéis SVDquartets\n'
bcftools query -f '%CHROM\t%POS\n' "$VCF" > "$SVD/input/all.positions.tsv"

for DIST in 100 10000
do
    SEED=20261021
    if [[ "$DIST" -eq 10000 ]]; then
        SEED=20261020
    fi

    STEM="dyckia_final202.mindist${DIST}bp"
    SITES="$SVD/input/$STEM.positions.tsv"
    PANEL="$SVD/input/$STEM.vcf.gz"
    NEXUS="$SVD/input/vcf2phylip/$STEM.min1.nexus"

    printf '\nPainel %s bp | semente da seleção: %s\n' "$DIST" "$SEED"
    python "$SCRIPTS/18_random_min_distance_sites.py" \
        "$SVD/input/all.positions.tsv" "$DIST" "$SEED" "$SITES"

    bcftools view -T "$SITES" -Oz -o "$PANEL" "$VCF"
    bcftools index "$PANEL"
    bcftools query -f '%CHROM\t%POS\n' "$PANEL" > "$SVD/input/$STEM.observed.tsv"
    cmp "$SITES" "$SVD/input/$STEM.observed.tsv"

    NSITES=$(bcftools index -n "$PANEL")
    [[ "$NSITES" -gt 0 ]]
    python "$CONVERTER" -i "$PANEL" -m 1 -n -p \
        --output-folder "$SVD/input/vcf2phylip" \
        --output-prefix "$STEM"

    grep -Eiq 'ntax[[:space:]]*=[[:space:]]*202([^0-9]|$)' "$NEXUS"
    grep -Eiq "nchar[[:space:]]*=[[:space:]]*${NSITES}([^0-9]|$)" "$NEXUS"
    printf 'Painel validado: 202 terminais e %s SNPs\n' "$NSITES"
done

printf '\n### ML final\n'
bash "$RUN/run_ml.sh"
[[ -s "$ML/dyckia_final202.mfp_asc.treefile" ]]
[[ -s "$ML/dyckia_final202.mfp_asc.contree" ]]

printf '\n### SVDquartets final\n'
bash "$RUN/run_svd.sh"

printf '\nTRÊS ANÁLISES CONCLUÍDAS: %s\n' "$RUN"
