#!/bin/bash
#SBATCH --job-name=piawkaPlotting
#SBATCH --partition=pibu_el8 
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G 
#SBATCH --time=3:00:00
#SBATCH --array=0-19%4
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step4_plotting/logFiles/1_multirun_piawka_R_visualization_w_boxplots/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step4_plotting/logFiles/1_multirun_piawka_R_visualization_w_boxplots/error_%A_%a.e

#runs analyses on mulitple genes' piawka output in R and visualizes it
#to run: sbatch path/to/file
#alter SBATCH --array as needed
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

WORKDIR="/data/users/aikiror/researchProject/pipeline"
RSCRIPT="${WORKDIR}/step4_plotting/scripts/1_piawka_R_visualization_w_boxplots.R"

#alter if need be
#output folder (BASE - will create subfolders per VARIA ID)
OUTPUT_DIR="${WORKDIR}/step4_plotting/outputFiles/1_multirun_piawka_R_visualization_w_boxplots"

#output from piawka
PIAWKA_0FOLD="${WORKDIR}/step3_popGenStats/outputFiles/14.95_piawka_stats_tajimasD_full_bedfile_0fold/stats_0fold.tsv"
PIAWKA_4FOLD="${WORKDIR}/step3_popGenStats/outputFiles/14.97_piawka_stats_tajimasD_full_bedfile_4fold/stats_4fold.tsv"

#bedfiles 0fold and 4fold
BED_0FOLD="${WORKDIR}/step3_popGenStats/outputFiles/11.5_alt_combine_bed_file/combined_0fold_Bv1-9_scaffold31.bed"
BED_4FOLD="${WORKDIR}/step3_popGenStats/outputFiles/11.5_alt_combine_bed_file/combined_4fold_Bv1-9_scaffold31.bed"

#actual gene info (coords: start, stop, etc)
ACTUAL_GENE_BOUNDARY_LIST="${WORKDIR}/biscutellaVariaData/geneCoord.tsv"

#table with geneic region w flanks start and end info
GENEIC_REGION_BOUNDARY_LIST="${WORKDIR}/biscutellaVariaData/flanking_2genes_gene_coords.tsv"

#flanking list for varia of interest -> list of genes in window
FLANKINGGENES_LIST="${WORKDIR}/biscutellaVariaData/flanking_2genes_gene_coords.tsv"


#coverage dir
COVERAGE_DIR="${WORKDIR}/step2_calcCoverageDepth/outputFiles/1_samtools_depth_per_region"

#build list of flanking files (adjust pattern if needed)
# e.g. flanking_list_for_VARIA_0007106.tsv, flanking_list_for_VARIA_XXXXXXX.tsv
mapfile -t FLANK_FILES < <(ls -1 "${FLANKINGGENES_LIST}"/flanking_list_for_VARIA_*.tsv | sort)

#pick file for this array task
FLANKING_FILE="${FLANK_FILES[$SLURM_ARRAY_TASK_ID]}"

# Safety checks
if [[ -z "${FLANKING_FILE}" ]]; then
  echo "ERROR: No file for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
  echo "Counted ${#FLANK_FILES[@]} files in ${FLANKINGGENES_LIST}"
  exit 1
fi

# *** EXTRACT VARIA ID FROM FILENAME ***
# From: flanking_list_for_VARIA_0007106.tsv
# Extract: VARIA_0007106
VARIA_ID=$(basename "${FLANKING_FILE}" | sed 's/flanking_list_for_//; s/.tsv$//')

# *** CREATE VARIA-SPECIFIC OUTPUT SUBFOLDER ***
VARIA_OUTPUT_DIR="${OUTPUT_DIR}/${VARIA_ID}"
mkdir -p "${VARIA_OUTPUT_DIR}"

echo "========================================"
echo "Task ${SLURM_ARRAY_TASK_ID}"
echo "VARIA ID: ${VARIA_ID}"
echo "Flanking file: ${FLANKING_FILE}"
echo "Output directory: ${VARIA_OUTPUT_DIR}"
echo "========================================"

#master_flanking_windows -> all genes in all the windows
master_flanking_list="${WORKDIR}/step1_prep/outputFiles/10_flanking_list/master_flanking_list.tsv"

#load R module - packages used in Rscript should be available on CRAN
module load R-bundle-CRAN/2023.11-foss-2021a

#run - pass VARIA-SPECIFIC output directory
Rscript "${RSCRIPT}" \
    "${PIAWKA_0FOLD}" \
    "${BED_0FOLD}" \
    "${PIAWKA_4FOLD}" \
    "${BED_4FOLD}" \
    "${ACTUAL_GENE_BOUNDARY_LIST}" \
    "${GENEIC_REGION_BOUNDARY_LIST}" \
    "${FLANKING_FILE}" \
    "${master_flanking_list}" \
    "${COVERAGE_DIR}" \
    "${VARIA_OUTPUT_DIR}"