#!/bin/bash
#SBATCH --array=1-9                    #run for biggest contigs: Bv1-Bv9
#SBATCH --time=48:00:00                
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=VCF0foldContig 
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/08_VCF_0fold_extract_by_contig/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/08_VCF_0fold_extract_by_contig/error_%A_%a.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#array job to extract vcf at 0-fold degenerate sites
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of contigs
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#safety
set -euo pipefail

#directories
WORKDIR="/data/users/aikiror/researchProject/pipeline"
CONTAINER="/containers/apptainer/gatk4_4.6.2.0_samtools_1.21.sif"

# BED, FASTA, and VCF inputs
BED_DIR="${WORKDIR}/step3_popGenStats/outputFiles/07_generate0FoldBedByContig/bed_per_contig"
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"
VCF="${WORKDIR}/step3_popGenStats/outputFiles/01_VCFFileFiltering/missingness_filtered.vcf.gz"

# Output directories
OUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/08_VCF_0fold_extract_by_contig"
mkdir -p "${OUTDIR}" 


# Contig name list (generated from reference)
CONTIG_NAMES_LIST="${WORKDIR}/step3_popGenStats/outputFiles/07_generate0FoldBedByContig/contigNamesList.txt"

#path to contig fasta = query
#assign the line that corresponds to the array_task_id: taskID=1 => 1st line 
CONTIG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${CONTIG_NAMES_LIST}")
BED="${BED_DIR}/${CONTIG}.bed"
OUT="${OUTDIR}/${CONTIG}_0fold.vcf.gz"

echo "[$(date)] Extracting 0-fold variants for ${CONTIG}"

apptainer exec --bind /data/ ${CONTAINER} gatk SelectVariants \
  -R "${FASTA}" \
  -V "${VCF}" \
  -O "${OUT}" \
  -L "${BED}" \
  --select-type-to-include SNP \
  --restrict-alleles-to BIALLELIC \
  --verbosity INFO

echo "[$(date)] Done ${CONTIG}"
