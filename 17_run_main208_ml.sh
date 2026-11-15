#!/usr/bin/env bash
set -uo pipefail

PROJECT=/home/jrabelo/tese/filogenia_dyckia
MLDIR="$PROJECT/05_gene_trees/concatenated_ml"

FASTA="$MLDIR/dyckia_main208.iupac.variable.fasta"
PREFIX="$MLDIR/dyckia_main208.mfp_asc"
THREADS=16
SEED=20261007

IQTREE=$(command -v iqtree2 || command -v iqtree)

if [[ -z "${IQTREE:-}" ]]
then
    echo "ERRO: IQ-TREE não encontrado." >&2
    exit 1
fi

if [[ ! -s "$FASTA" ]]
then
    echo "ERRO: alinhamento ausente ou vazio: $FASTA" >&2
    exit 1
fi

if [[ -s "${PREFIX}.treefile" && -s "${PREFIX}.contree" ]]
then
    echo "Árvore principal já concluída."
    echo "Treefile: ${PREFIX}.treefile"
    echo "Contree: ${PREFIX}.contree"
    exit 0
fi

echo "=================================================="
echo "FILOGENIA ML PRINCIPAL — DYCKIA"
echo "Início: $(date)"
echo "IQ-TREE: $IQTREE"
echo "Alinhamento: $FASTA"
echo "Threads: $THREADS"
echo "Modelo: MFP+ASC"
echo "UFBoot: 1000"
echo "SH-aLRT: 1000"
echo "=================================================="

echo
echo "Comando executado:"
echo "$IQTREE -s $FASTA -st DNA -m MFP+ASC -B 1000 --alrt 1000 --bnni -T $THREADS --seed $SEED --prefix $PREFIX"
echo

/usr/bin/time -v "$IQTREE" \
    -s "$FASTA" \
    -st DNA \
    -m MFP+ASC \
    -B 1000 \
    --alrt 1000 \
    --bnni \
    -T "$THREADS" \
    --seed "$SEED" \
    --prefix "$PREFIX"

STATUS=$?

echo
echo "=================================================="
echo "Fim: $(date)"
echo "Status: $STATUS"
echo "=================================================="

if [[ "$STATUS" -eq 0 ]]
then
    echo "ÁRVORE PRINCIPAL CONCLUÍDA"
    echo "Treefile: ${PREFIX}.treefile"
    echo "Contree: ${PREFIX}.contree"
    echo "Relatório: ${PREFIX}.iqtree"
else
    echo "ERRO NA EXECUÇÃO — status $STATUS" >&2
fi

exit "$STATUS"
