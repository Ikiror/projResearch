#!/bin/bash
#SBATCH --array=1-12
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=GeneFilt
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/3_filterForGenes/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/3_filterForGenes/error_%A_%a.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#array job to extract rows from an annotation feature map that match a protein ID listed in a 2-column sample list.
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of proteins
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTDIR="${WORKDIR}/step1_prep/outputFiles/3_filterForGenes"
mkdir -p $OUTDIR

#annot feature map to search in
INPUT="${WORKDIR}/biscutellaVariaData/annot_feature_map_VARIA.txt"

#list with proteins and corresponding uniProtIDs -> need to create manually -> gene uniProtID
SAMPLELIST="${WORKDIR}/biscutellaVariaData/ProteinID_sampleList.tsv"

#extract the gene + protein ID for this array task
GENE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
PROTEIN_ID=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`

#filter annot map for matching protein id 
grep ${PROTEIN_ID} ${INPUT} > ${OUTDIR}/filteredFor${GENE}_${PROTEIN_ID}.tsv
