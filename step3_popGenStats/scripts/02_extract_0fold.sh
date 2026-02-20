#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=0foldExtract
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/02_extract_0fold/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/02_extract_0fold/error_%J.e
#SBATCH --partition=pibu_el8


#to run: sbatch path/to/file
#this job filters for the 0-fold degenerate sites in the degenerate bed file
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

# Define variables
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/02_extract_0fold"
mkdir -p $OUTPUTDIR

#degeneracy bed file
INPUT="${WORKDIR}/degeneracyInfo/degeneracy-all-sites.bed"

#filter for 0-fold
awk '$5 == "0"' $INPUT > ${OUTPUTDIR}/0fold_degenerate_sites.bed

