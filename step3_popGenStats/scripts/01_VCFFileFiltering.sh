#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=8
#SBATCH --job-name=VCFfilter
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/01_VCFFileFiltering/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/01_VCFFileFiltering/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this job filters the vcf file for missingness
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/01_VCFFileFiltering"
mkdir -p $OUTPUTDIR

#container
CONTAINER="/containers/apptainer/bcftools:1.21--h8b25389_0"

#path to vcf file
VCFFILE="${WORKDIR}/vcf_file/merged.vcf.gz"

#filter vcf for missingness
apptainer exec --bind /data/ ${CONTAINER} bash -c "bcftools view --threads 8 -O v -i 'F_MISSING <= 0.5' ${VCFFILE} | bcftools view --threads 8 -O z > ${OUTPUTDIR}/missingness_filtered.vcf.gz"