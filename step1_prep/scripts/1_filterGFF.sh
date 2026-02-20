#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=filterGFF
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/1_filterGFF/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/1_filterGFF/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this script will filter a gff file to only the gene features
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"

#output dir
OUTDIR="${WORKDIR}/step1_prep/outputFiles/1_filterGFF"
mkdir -p $OUTDIR

#gff 
INPUT="${WORKDIR}/biscutellaVariaData/Varia.all.maker.renamed.gff"

#filter gff file
awk '$3 == "gene"' $INPUT > $OUTDIR/filteredOutput.gff
