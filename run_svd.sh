#!/usr/bin/env bash
set -euo pipefail

PROJECT=/home/jrabelo/tese/filogenia_dyckia
SVD="/home/jrabelo/tese/filogenia_dyckia/final202_trees.2jgge2/svdquartets"
PAUP="$PROJECT/software/paup/paup"

RESULTDIR="$SVD/results/bootstrap100"
LOGDIR="$SVD/logs/bootstrap100"
CMDDIR="$SVD/input/bootstrap100"

mkdir -p "$RESULTDIR" "$LOGDIR" "$CMDDIR"

for DIST in 100 10000
do
    NEXUS="$SVD/input/vcf2phylip/dyckia_final202.mindist${DIST}bp.min1.nexus"
    TREE="$RESULTDIR/dyckia_final202.mindist${DIST}bp.bootstrap100.tre"
    LOG="$LOGDIR/mindist${DIST}bp.bootstrap100.log"
    CMD="$CMDDIR/mindist${DIST}bp.bootstrap100.paup"

    if [[ -s "$TREE" ]]
    then
        echo "JÁ CONCLUÍDO: mindist${DIST}bp"
        continue
    fi

    if [[ "$DIST" == "10000" ]]
    then
        SEED=20261020
    else
        SEED=20261021
    fi

    cat > "$CMD" <<EOF
execute $NEXUS;

svdquartets
    evalquartets=all
    taxpartition=none
    treeinf=QFM
    bootstrap=standard
    nreps=100
    seed=$SEED
    nthreads=16
    ambigs=distribute
;

savetrees
    file=$TREE
    replace=yes
;

quit;
EOF

    echo "=================================================="
    echo "INÍCIO: mindist${DIST}bp — $(date)"
    echo "100 bootstraps; todos os 67.331.650 quartetos"
    echo "=================================================="

    /usr/bin/time -v "$PAUP" -n \
        < "$CMD" \
        > "$LOG" 2>&1

    if [[ ! -s "$TREE" ]]
    then
        echo "ERRO: árvore não produzida: $TREE" >&2
        exit 1
    fi

    if ! grep -q \
        'Bootstrap consensus tree stored' \
        "$LOG"
    then
        echo "ERRO: conclusão do bootstrap não confirmada." >&2
        exit 1
    fi

    echo "CONCLUÍDO: mindist${DIST}bp — $(date)"
done

echo
echo "TODOS OS BOOTSTRAPS FORAM CONCLUÍDOS."
