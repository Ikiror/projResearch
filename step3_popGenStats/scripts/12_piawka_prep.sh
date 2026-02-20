#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --job-name=piawkaPrep
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/12_piawka_prep/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/12_piawka_prep/error_%J.e
#SBATCH --partition=pibu_el8

#extracts population sample names - 
#will need to manually add other relevant information after generation
#e.g ploidy: diploid, and elevation: high elevation
#into format: sampleName ploidy_elevation (i.e RAM-3	high_tetraploid) - file also in biscutellaVariaData : sample_name_combined_elevation_ploidy.tsv
#to run: sbatch path/to/file
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#directories
WORKDIR="/data/users/aikiror/researchProject/pipeline"
CONTAINER="/containers/apptainer/bcftools:1.21--h8b25389_0"
OUTPUT_DIR="${WORKDIR}/step3_popGenStats/outputFiles/12_piawka_prep"
mkdir -p ${OUTPUT_DIR}

#files
VCF_FILE=${WORKDIR}/step3_popGenStatss/outputFiles/05_VCF_0fold_extract/0_fold_degenerate_filtered.vcf.gz
TBI_FILE=${WORKDIR}/step3_popGenStats/outputFiles/05_VCF_0fold_extract/*.tbi

SAMPLE_NAME_TXT="${OUTPUT_DIR}/sample_names.txt"

#to generate sample name txt
apptainer exec --bind /data/ ${CONTAINER} bcftools query -l ${VCF_FILE} > ${SAMPLE_NAME_TXT}

#then manually added elevation information to sample_names.txt
