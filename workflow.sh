#!/bin/bash

PROJECT_ROOT=$(git rev-parse --show-toplevel)

CLUSTER_NAME="$1"
CLUSTER="${CLUSTER_NAME// /_}"



DATA_DIR=$PROJECT_ROOT/data/

mkdir -p $DATA_DIR
BASE_OUT=$PROJECT_ROOT/output/$CLUSTER/
PROKKA=$BASE_OUT/prokka
PANA=$BASE_OUT/panaroo
TMP=$BASE_OUT/tmp

mkdir -p $BASE_OUT
mkdir -p $PROKKA
mkdir -p $PANA
mkdir -p $TMP

# gtdb-dl --gtdb r226 --taxon "$CLUSTER_NAME" --flat species --output $DATA_DIR -v --flag-rep 


# for g in $DATA_DIR/$CLUSTER/*.fna.gz; do
#     base=$(basename "$g" .fna.gz)

#     TMP_FA=$TMP/$base.fna
#     zcat $g > $TMP_FA

#     echo "Output dir: $PROKKA"
#     echo "Prefix:     $base"

#     prokka \
#         --outdir $PROKKA \
#         --prefix $base \
#         --kingdom Bacteria \
#         --metagenome \
#         --cpus 8 \
#         --force \
#         "$TMP_FA"
#     rm $TMP_FA
# done

# rm -rf $TMP

# panaroo -i $PROKKA/*.gff -o $PANA --clean-mode moderate -t 8

# Filter panaroo FASTA files to keep only genes in representative genome
echo ""
echo "Building panaroo class FASTA files from representative genome sequences..."
python3 "$PROJECT_ROOT/build_speciesrep_class_fastas.py" "$PANA" "$PROKKA"

# mv $PANA/pan_genome_reference.* $CLUSTER

echo ""
echo "Pipeline complete!"
echo "Results:"
echo "  - Panaroo output: $PANA"
echo "  - Core genes table: $CORE_GENES"
echo "  - Filtered FASTA files:"
echo "    - $PANA/pan_genome_reference.core.speciesrep.fasta"
echo "    - $PANA/pan_genome_reference.shell.speciesrep.fasta"
echo "    - $PANA/pan_genome_reference.cloud.speciesrep.fasta"
echo "    - $PANA/pan_genome_reference.soft_core.speciesrep.fasta"


