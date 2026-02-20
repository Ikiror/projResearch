#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sortGFF
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/2_sortFilteredGFF/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/2_sortFilteredGFF/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#Sort filteredOutput.gff by scaffold/contig and then by genomic coordinates
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTDIR="${WORKDIR}/step1_prep/outputFiles/2_sortFilteredGFF"
INPUT="${WORKDIR}/step1_prep/outputFiles/1_filterGFF/filteredOutput.gff"
#create outputdir if it doenst exist
mkdir -p $OUTDIR

cd $OUTDIR

#sort by contig; and by start coord
sort -k1,1V -k4,4n $INPUT > sorted_filteredOutput.gff
