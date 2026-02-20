#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --job-name=contig4FoldBed
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/09_generate4FoldBedByContig/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/09_generate4FoldBedByContig/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this job filters for the 4-fold degenerate sites in the degenerate bed file
#it makes bed files per contig
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail


##getting 4fold bed file by contig
set -euo pipefail

#define file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"

OUTPUT_DIR="${WORKDIR}/step3_popGenStats/outputFiles/09_generate4FoldBedByContig"
mkdir -p ${OUTPUT_DIR}

#output for the contigs' bed files
BED_BY_CONTIG_OUTPUT="${OUTPUT_DIR}/bed_per_contig"
mkdir -p ${BED_BY_CONTIG_OUTPUT}

#output for the contig names list - already generated while generating 0fold bed by contig
CONTIG_NAMES="${WORKDIR}/step3_popGenStats/outputFiles/07_generate0FoldBedByContig/contigNamesList.txt"

# #if not already generated:
# CONTIG_NAMES="${OUTPUT_DIR}/contigNamesList.txt"
# #get names of contigs
# echo "[$(date)] Extracting contig names from FASTA..."
# grep ">" "${FASTA}" | cut -d' ' -f1 | sed 's/>//' > "${CONTIG_NAMES}"
# echo "[$(date)] Found $(wc -l < "${CONTIG_NAMES}") contigs in reference."

#fasta seq
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"

#4fold bed sites
BED4FOLD="${WORKDIR}/step3_popGenStats/outputFiles/03_extract_4fold/4fold_degenerate_sites.bed"



echo "[$(date)] Splitting BED into one file per contig..."
awk -v outdir="${BED_BY_CONTIG_OUTPUT}" 'BEGIN{OFS="\t"} {print > (outdir "/" $1 ".bed")}' "${BED4FOLD}"
#BEGIN{OFS="\t"} - keeps \t (tab spacing) intact
echo "[$(date)] BED splitting complete!"
echo "[$(date)] Example output files:"
ls -lh "${BED_BY_CONTIG_OUTPUT}" | head
