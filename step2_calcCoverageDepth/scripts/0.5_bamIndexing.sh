#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --job-name=bamIndex
#SBATCH --output=/data/users/aikiror/researchProject/step3_coveragePlots/logFiles_Reports/02.70_bamIndexing/output_02.70_bamIndexing_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/step3_coveragePlots/logFiles_Reports/02.70_bamIndexing/error_02.70_bamIndexing_%J.e
#SBATCH --partition=pibu_el8

#creates index for bam file
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of proteins
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#output goes to location of .bam file
set -euo pipefail

#tmpdir set up
export TMPDIR="/data/users/aikiror/tmp"
mkdir -p "$TMPDIR"

CONTAINER="/containers/apptainer/samtools-1.19.sif"
WORKDIR="/data/users/aikiror/researchProject/pipeline"
BAMFILE="${WORKDIR}/biscutellaVariaData/VARI_short_reads.bam"

apptainer exec --bind /data --bind $TMPDIR:/tmp $CONTAINER samtools index $BAMFILE