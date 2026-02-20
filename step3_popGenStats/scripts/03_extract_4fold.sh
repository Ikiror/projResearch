#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=4foldExtract
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/03_extract_4fold/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/03_extract_4fold/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this job filters for the 4-fold degenerate sites in the degenerate bed file
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

# Define variables
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/03_extract_4fold"
mkdir -p $OUTPUTDIR

#degeneracy bed file
INPUT="${WORKDIR}/degeneracyInfo/degeneracy-all-sites.bed"

#filter for 4-fold
awk '$5 == "4"' $INPUT > ${OUTPUTDIR}/4fold_degenerate_sites.bed

