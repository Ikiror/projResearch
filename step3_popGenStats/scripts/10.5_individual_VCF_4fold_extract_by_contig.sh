#!/bin/bash
#SBATCH --time=48:00:00                
#SBATCH --cpus-per-task=4
#SBATCH --job-name=VCF4foldContig 
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/10.5_individual_VCF_4fold_extract_by_contig/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/10.5_individual_VCF_4fold_extract_by_contig/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch /path/to/script contig_name
#this script will run on an individual contig/scaffold as opposed to multiple as in 10_VCF_4fold_extract_by_contig.sh
#job to extract vcf at 4-fold degenerate sites
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#safety
set -euo pipefail

#directories
WORKDIR="/data/users/aikiror/researchProject/pipeline"
CONTAINER="/containers/apptainer/gatk4_4.6.2.0_samtools_1.21.sif"

# BED, FASTA, and VCF inputs
BED_DIR="${WORKDIR}/step3_popGenStats/outputFiles/09_generate4FoldBedByContig/bed_per_contig"
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"
VCF="${WORKDIR}/step3_popGenStats/outputFiles/01_VCFFileFiltering/missingness_filtered.vcf.gz"

# Output directories
#script is 10.5 but output same as 10 script
OUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/10_VCF_4fold_extract_by_contig"
mkdir -p "${OUTDIR}"

#argument passed in command line script will be the contig used.
CONTIG=$1
BED="${BED_DIR}/${CONTIG}.bed"
OUT="${OUTDIR}/${CONTIG}_4fold.vcf.gz"

echo "[$(date)] Extracting 4-fold variants for ${CONTIG}"

apptainer exec --bind /data/ ${CONTAINER} gatk SelectVariants \
  -R "${FASTA}" \
  -V "${VCF}" \
  -O "${OUT}" \
  -L "${BED}" \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  --verbosity INFO

echo "[$(date)] Done ${CONTIG}"

