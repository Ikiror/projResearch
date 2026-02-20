#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --array=1-20%2
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=8
#SBATCH --job-name=samtoolsDepth
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step2_calcCoverageDepth/logFiles/1_samtools_depth_per_region/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step2_calcCoverageDepth/logFiles/1_samtools_depth_per_region/error_%A_%a.e
#SBATCH --partition=pibu_el8

#calculates the coverage depth per 5-gene window
#alter SBATCH --output,error, workdir, tmpdir to match file paths
#alter SBATCH --array to match number of items to be analyzed
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

set -euo pipefail

#tmpdir set up
export TMPDIR="/data/users/aikiror/tmp"
mkdir -p "$TMPDIR"

#container
CONTAINER="/containers/apptainer/samtools-1.19.sif"

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUT_DIR="${WORKDIR}/step2_calcCoverageDepth/outputFiles/1_samtools_depth_per_region"
#bam file
BAMFILE="${WORKDIR}/biscutellaVariaData/VARI_short_reads.bam"
REGION_LIST="${WORKDIR}/biscutellaVariaData/flanking_2genes_gene_coords.tsv"

mkdir -p "${OUTPUT_DIR}"

cd $OUTPUT_DIR

CHR=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i {print $1}' $REGION_LIST)
FROM=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i {print $2}' $REGION_LIST)
TO=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i {print $3}' $REGION_LIST)
GENE=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i {print $4}' $REGION_LIST)
SEQNAME=$(awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i {print $5}' $REGION_LIST)

OUTFILE="${OUTPUT_DIR}/${SEQNAME}_${CHR}_${GENE}_coverage.tsv"


echo "[$(date)] Running samtools depth for $CHR:$FROM-$TO ..."

apptainer exec \
    --bind /data \
    --bind $TMPDIR:/tmp \
    "$CONTAINER" \
    samtools depth -r ${CHR}:${FROM}-${TO} "$BAMFILE" \
    > "$OUTFILE"

echo "[$(date)] Finished region: $CHR:$FROM-$TO"
echo "Output: $OUTFILE"