#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --job-name=fai
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/5_fasta_indexing/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/5_fasta_indexing/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#array job to generate the index of the fasta file
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail
#to run: sbatch path/to/file

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"

#container
CONTAINER="/containers/apptainer/gatk4_4.6.2.0_samtools_1.21.sif"

#generate index
apptainer exec --bind /data/ ${CONTAINER} samtools faidx ${FASTA}

#echo path to .fai
echo ".fai file in same directory as ${FASTA}"

