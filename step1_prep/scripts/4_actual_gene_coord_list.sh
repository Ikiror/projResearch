#!/bin/bash
#SBATCH --array=1-20
#SBATCH --time=00:20:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sampleListGenerator
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/4_actual_gene_coord_list/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/4_actual_gene_coord_list/error_%A_%a.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#generates sample list with info about genes -> contig geneStart geneStop gene geneID
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of items being analyzed
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

WORKDIR="/data/users/aikiror/researchProject/pipeline"

GFF="${WORKDIR}/step1_prep/outputFiles/2_sortFilteredGFF/sorted_filteredOutput.gff"
OUT="${WORKDIR}/biscutellaVariaData/geneCoord.tsv"
GENEANNOTFTMAPDIR="${WORKDIR}/step1_prep/outputFiles/3_filterForGenes"

awk -F'\t' '
  BEGIN { OFS="\t" }

  # ---------- PASS 1: read GFF, store coordinates keyed by VARIA ID ----------
  FNR==NR {
    if ($0 ~ /^#/ ) next
    if ($3 != "gene") next

    attr=$9
    id=""
    if (match(attr, /ID=[^;]+/)) id=substr(attr, RSTART+3, RLENGTH-3)
    else next

    contig[id]=$1
    start[id]=$4
    stop[id]=$5
    next
  }

  # ---------- PASS 2: read featuremap lines ----------
  {
    # column 1: VARIA_... or VARIA_...-RA
    raw=$1
    id=raw
    sub(/-R[A-Z0-9]+$/, "", id)   # strips -RA, -RB, etc if present

    # Extract gene name from filename: filteredForAP1_P35631.tsv
    gene = FILENAME
    sub(/^.*filteredFor/, "", gene)
    sub(/_.*/, "", gene)

    if (!(id in contig)) next  # skip IDs not found in GFF

    print contig[id], start[id], stop[id], gene, id
  }
' "$GFF" $GENEANNOTFTMAPDIR/filteredFor*.tsv \
| awk 'BEGIN{OFS="\t"} {key=$5; if(!seen[key]++){print}}' \
> "$OUT"

#to run: sbatch path/to/file