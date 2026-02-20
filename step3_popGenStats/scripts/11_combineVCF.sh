#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=combineVCF
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/11_combineVCF/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/11_combineVCF/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this script takes multiple VCFs and combines them; 
#and takes multiple bed files and combines them
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#load bcftools module
module load BCFtools/1.12-GCC-10.3.0

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/11_combineVCF"
mkdir -p "${OUTPUTDIR}" #make outputdir if it doesnt already exist

#location of VCF files to be combined
VCF_FILE_0FOLD="${WORKDIR}/step3_popGenStats/outputFiles/08_VCF_0fold_extract_by_contig"
VCF_FILE_4FOLD="${WORKDIR}/step3_popGenStats/outputFiles/10_VCF_4fold_extract_by_contig"

#output prefix
MERGED_0FOLD="merged_0fold_bv1_9_and_scaffold31.vcf.gz"
MERGED_4FOLD="merged_4fold_bv1_9_and_scaffold31.vcf.gz"

#location of bed files to be combined
BED_FILES="${WORKDIR}/step1_prep/outputFiles/11_bedtools_blast_alignment_plot/bedtools_output/extracted_sequences_geneic_window/*.bed"

cd "${OUTPUTDIR}"

#concat 0fold vcf files and index them
bcftools concat --allow-overlaps -Oz --threads 4 -o "${MERGED_0FOLD}" "${VCF_FILE_0FOLD}"/*vcf.gz
bcftools index -t "${MERGED_0FOLD}"

#concat 4fold vcf files and index them
bcftools concat --allow-overlaps -Oz --threads 4 -o "${MERGED_4FOLD}" "${VCF_FILE_4FOLD}"/*vcf.gz
bcftools index -t "${MERGED_4FOLD}"

#combine bed files
cat ${BED_FILES} | sort -k1,1 -k2,2n > combined.bed