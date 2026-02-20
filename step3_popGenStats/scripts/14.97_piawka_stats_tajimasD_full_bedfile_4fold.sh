#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --mem=512gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=TajimasD4fold
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/14.97_piawka_stats_tajimasD_full_bedfile_4fold/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/14.97_piawka_stats_tajimasD_full_bedfile_4fold/error_%J.e
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#to run: sbatch path/to/file
#job to run piawka on 4-fold degenerate sites
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#safety
set -euo pipefail

#make piawka scripts available in jobs; make gawk5 available
module load BCFtools/1.12-GCC-10.3.0
module load HTSlib/1.12-GCC-10.3.0 #error indicated need for bgzip and tabix which htslib should have
export PATH="/data/users/aikiror/researchProject/pipeline/tools/gawk5/bin:${PATH}"
export PATH="/data/users/aikiror/researchProject/pipeline/tools/piawka/scripts:${PATH}"

#directories
WORKDIR="/data/users/aikiror/researchProject/pipeline"

OUTPUT_DIR="${WORKDIR}/step3_popGenStats/outputFiles/14.97_piawka_stats_tajimasD_full_bedfile_4fold"

SAMPLE_ID="${WORKDIR}/step3_popGenStats/outputFiles/12_piawka_prep/sample_name_combined_elevation_ploidy.tsv"

BED_FILE_4FOLD="$WORKDIR/step3_popGenStats/outputFiles/11.5_alt_combine_bed_file/combined_4fold_Bv1-9_scaffold31.bed"

#match vcf file to same contig as bed_file
VCF_4FOLD_FILE="${WORKDIR}/step3_popGenStats/outputFiles/11_combineVCF/merged_4fold_bv1_9_and_scaffold31.vcf.gz"

cd "${OUTPUT_DIR}"

# Weir and Cockerham Fst, rho, watterson theta, max missing GT per group @ site, output val for each site,
echo "[$(date)] Starting 4-fold analysis..."

piawka --groups ${SAMPLE_ID} \
    -v ${VCF_4FOLD_FILE} \
    -b ${BED_FILE_4FOLD} \
    -F \
    -r \
    -w \
    -T \
    -M 0.2 \
    -j 16 > ${OUTPUT_DIR}/stats_4fold.tsv

echo "[$(date)] Finishing 4-fold analysis..."

status0=$?
echo "[$(date)] 4-fold exit code: ${status0}"

#if [[ $status0 -eq 0 && $status4 -eq 0 ]]; then
if [[ $status0 -eq 0 ]]; then
    echo "[$(date)] Piawka finished successfully."
else
    echo "[$(date)] WARNING: Piawka exited with non-zero code(s): 0-fold=$status0 "
fi

#calc piawka for 0fold - calculate average stats = this will correspond to a line.
#extract piawka for geneic regions - plot that (boxplot)
#if geneic regions above line(avg for genome) - significant