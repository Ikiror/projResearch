#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --job-name=combineBed
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/11.5_alt_combine_bed_file/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/11.5_alt_combine_bed_file/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#this script takes multiple bed files and combines them
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail


#load bcftools module
module load BCFtools/1.12-GCC-10.3.0

WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/11.5_alt_combine_bed_file"
mkdir -p "${OUTPUTDIR}" #make outputdir if it doesnt already exist

#list of contig names
CONTIG_NAMES_LIST="${WORKDIR}/step3_popGenStats/outputFiles/07_generate0FoldBedByContig/contigNamesList.txt"
TARGET_CONTIGS=("Bv1" "Bv2" "Bv3" "Bv4" "Bv5" "Bv6" "Bv7" "Bv8" "Bv9" "scaffold_31")

#location of bed files to be combined
BEDFILE_0FOLD="$WORKDIR/step3_popGenStats/outputFiles/07_generate0FoldBedByContig/bed_per_contig"
BEDFILE_4FOLD="$WORKDIR/step3_popGenStats/outputFiles/09_generate4FoldBedByContig/bed_per_contig"

cd "${OUTPUTDIR}"

#combine bed files

#0fold
OUT_0FOLD="${OUTPUTDIR}/combined_0fold_Bv1-9_scaffold31.bed" > "${OUT_0FOLD}"  # create empty file

for contig in "${TARGET_CONTIGS[@]}"; do
    bedfile="${BEDFILE_0FOLD}/${contig}.bed"
    if [[ -f "$bedfile" ]]; then
        echo "Adding: $bedfile"
        cat "$bedfile" >> "${OUT_0FOLD}"
    else
        echo "missing: $bedfile" >&2
    fi
done

sort -k1,1 -k2,2n "${OUT_0FOLD}" -o "${OUT_0FOLD}"
echo "0fold combined bed file saved to ${OUT_0FOLD}"

#4fold
OUT_4FOLD="${OUTPUTDIR}/combined_4fold_Bv1-9_scaffold31.bed"
> "${OUT_4FOLD}"

for contig in "${TARGET_CONTIGS[@]}"; do
    bedfile="${BEDFILE_4FOLD}/${contig}.bed"
    if [[ -f "$bedfile" ]]; then
        echo "Adding: $bedfile"
        cat "$bedfile" >> "${OUT_4FOLD}"
    else
        echo "missing: $bedfile" >&2
    fi
done

sort -k1,1 -k2,2n "${OUT_4FOLD}" -o "${OUT_4FOLD}"
echo "4fold combined file saved to ${OUT_4FOLD}"

echo "Finished combining BEDs for Bv1–Bv9 and scaffold_31"

