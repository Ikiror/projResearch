#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=dictFile
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/04_dictFileCreation/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/04_dictFileCreation/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this job generates the fasta file's dict file 
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#define variables
WORKDIR="/data/users/aikiror/researchProject/pipeline"
#container
CONTAINER="/containers/apptainer/gatk4_4.6.2.0_samtools_1.21.sif"
#fasta file
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"
#OUTPUT will go to same location as FASTA file

#generate dict file
apptainer exec --bind /data/ ${CONTAINER} gatk CreateSequenceDictionary -R ${FASTA}