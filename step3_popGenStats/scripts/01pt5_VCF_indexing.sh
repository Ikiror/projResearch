#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=VCFindex
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/01pt5_VCF_indexing/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/01pt5_VCF_indexing/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this job indexes the missingness filtered vcf file
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
#container
CONTAINER="/containers/apptainer/bcftools:1.21--h8b25389_0"
#OUTPUT will go to same location as VCF file
VCF="${WORKDIR}/step3_popGenStats/outputFiles/01_VCFFileFiltering/missingness_filtered.vcf.gz"

#to generate .csi
apptainer exec --bind /data/ ${CONTAINER} bcftools index ${VCF}

#to generate .tbi
apptainer exec --bind /data/ ${CONTAINER} bcftools index -t ${VCF}