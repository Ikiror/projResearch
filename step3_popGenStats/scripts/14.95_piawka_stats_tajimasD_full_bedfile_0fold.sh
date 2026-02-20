#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=256gb
#SBATCH --cpus-per-task=8
#SBATCH --job-name=TajimasD0fold
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/14.95_piawka_stats_tajimasD_full_bedfile_0fold/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/14.95_piawka_stats_tajimasD_full_bedfile_0fold/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#job to run piawka on 0-fold degenerate sites
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

set -euo pipefail

#make piawka scripts available in jobs; make gawk5 available
module load BCFtools/1.12-GCC-10.3.0
module load HTSlib/1.12-GCC-10.3.0 #error indicated need for bgzip and tabix which htslib should have
export PATH="/data/users/aikiror/researchProject/pipeline/tools/gawk5/bin:${PATH}"
export PATH="/data/users/aikiror/researchProject/pipeline/tools/piawka/scripts:${PATH}"

#directories
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUT_DIR="${WORKDIR}/step3_popGenStats/outputFiles/14.95_piawka_stats_tajimasD_full_bedfile_0fold"
mkdir -p ${OUTPUT_DIR}

SAMPLE_ID="${WORKDIR}/step3_popGenStats/outputFiles/12_piawka_prep/sample_name_combined_elevation_ploidy.tsv"

BED_FILE_0FOLD="$WORKDIR/step3_popGenStats/outputFiles/11.5_alt_combine_bed_file/combined_0fold_Bv1-9_scaffold31.bed"

#match vcf file to same contig as bed_file
VCF_0FOLD_FILE="${WORKDIR}/step3_popGenStats/outputFiles/11_combineVCF/merged_0fold_bv1_9_and_scaffold31.vcf.gz"

cd "${OUTPUT_DIR}"

# Tajimas D, Weir and Cockerham Fst, rho, watterson theta, max missing GT per group @ site, output val for each site,
echo "[$(date)] Starting 0-fold analysis..."

piawka --groups ${SAMPLE_ID} \
    -v ${VCF_0FOLD_FILE} \
    -b ${BED_FILE_0FOLD} \
    -F \
    -r \
    -w \
    -T \
    -M 0.2 \
    -j 8 > ${OUTPUT_DIR}/stats_0fold.tsv

echo "[$(date)] Finishing 0-fold analysis..."

status0=$?
echo "[$(date)] 0-fold exit code: ${status0}"

#if [[ $status0 -eq 0 && $status4 -eq 0 ]]; then
if [[ $status0 -eq 0 ]]; then
    echo "[$(date)] Piawka finished successfully."
else
    echo "[$(date)] WARNING: Piawka exited with non-zero code(s): 0-fold=$status0 "
fi

#calc piawka for 0fold - calculate average stats = this will correspond to a line.
#extract piawka for geneic regions - plot that (boxplot)
#if geneic regions above line(avg for genome) - significant