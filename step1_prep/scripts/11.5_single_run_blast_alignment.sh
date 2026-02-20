#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=blastPlot
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/11.5_single_run_blast_alignment/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/11.5_single_run_blast_alignment/error_%J.e
#SBATCH --partition=pibu_el8
#SBATCH --mail-user=amo.ikiror@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

#to run: sbatch path/to/file
#job for a single gene: to create a bed file, extract sequence, compare sequence window to blast database, plot alignment
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

set -euo pipefail

echo "loading files..."

#target
TARGET_GENE_ID="VARIA_0021468"


WORKDIR="/data/users/aikiror/researchProject/pipeline"
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"
OUTPUTDIR="${WORKDIR}/step1_prep/outputFiles/11.5_single_run_blast_alignment"

BLAST_OUTPUTDIR="${WORKDIR}/step1_prep/outputFiles/9_blast_db/blast_output" #location of /genomedb;
BLAST_DB="${BLAST_OUTPUTDIR}/genome_db"

BLAST_ALIGNMENT_OUTPUTDIR="${OUTPUTDIR}/blast_alignments" #self_alignment.tsv's
BEDTOOLS_OUTPUTDIR="${OUTPUTDIR}/bedtools_output" #location of bedfiles
SEQUENCE_OUTPUTDIR="${BEDTOOLS_OUTPUTDIR}/extracted_sequences_geneic_window" #location of extracted seqs
R_OUTDIR="${OUTPUTDIR}/visualization_of_alignment"

mkdir -p $BEDTOOLS_OUTPUTDIR
mkdir -p $BLAST_ALIGNMENT_OUTPUTDIR
mkdir -p $SEQUENCE_OUTPUTDIR
mkdir -p $R_OUTDIR

#check if it blast db exists
if [ ! -f "${BLAST_DB}.nhr" ]; then
    echo "ERROR: BLAST database not found at ${BLAST_DB}"
    echo "Please run step2a_blast_db_setup.sh first!"
    exit 1
fi

#info from master list geneic window 
SAMPLELIST="${WORKDIR}/biscutellaVariaData/flanking_2genes_gene_coords.tsv"



#extract info from id
read -r seqName start stop geneName geneID < <(
  awk -v id="$TARGET_GENE_ID" '$5==id {print $1, $2, $3, $4, $5; exit}' "$SAMPLELIST"
)

# safety check
if [[ -z "${geneID:-}" ]]; then
  echo "ERROR: geneID '$TARGET_GENE_ID' not found in $SAMPLELIST" >&2
  exit 1
fi

echo "processing: ${geneID} (${geneName})"
echo "region: ${seqName}:${start}-${stop}"

#list of genes in specific geneic window
FLANKING_DIR="${WORKDIR}/step1_prep/outputFiles/10_flanking_list"
FLANKING_FILE="${FLANKING_DIR}/flanking_list_for_${geneID}.tsv"

#actual gene file
ACTUAL_GENE_FILE="$WORKDIR/biscutellaVariaData/geneCoord.tsv"

#naming of output files
OUTPUT_FILE_NAME="${geneID}_${geneName}" #to pass into R
FASTA_OUTPUTFILE="${SEQUENCE_OUTPUTDIR}/${OUTPUT_FILE_NAME}_window.fa"
BLAST_OUTPUTFILE="${BLAST_ALIGNMENT_OUTPUTDIR}/${OUTPUT_FILE_NAME}_self_alignment.tsv"


BEDTOOLS_CONTAINER="/containers/apptainer/bedtools_2.31.1.sif"
BLAST_CONTAINER="/containers/apptainer/blast_2.15.0.sif"


#load R module - packages used in Rscript should be available on CRAN
module load R-bundle-CRAN/2023.11-foss-2021a
#with visualization of script
RSCRIPT="${WORKDIR}/step1_prep/scripts/11_blast_alignment_plot.R"

echo "files loaded"

#EXTRACT SEQUENCES - BEDTOOLS
cd $BEDTOOLS_OUTPUTDIR
echo "extracting sequences of geneic windows..."

BEDFILE="${BEDTOOLS_OUTPUTDIR}/${OUTPUT_FILE_NAME}.bed"
echo -e "${seqName}\t${start}\t${stop}" > ${BEDFILE}

apptainer exec --bind /data/ ${BEDTOOLS_CONTAINER} \
  bedtools getfasta \
  -fi ${FASTA} \
  -bed ${BEDFILE} \
  -split | fold -w 60 > ${FASTA_OUTPUTFILE}


if [ ! -f "${FASTA_OUTPUTFILE}" ]; then
    echo "ERROR: Failed to extract sequence"
    exit 1
fi


echo "sequence extract done"
echo "bed file in ${BEDTOOLS_OUTPUTDIR}/${OUTPUT_FILE_NAME}"
echo "sequence extract: ${FASTA_OUTPUTFILE}"


echo "performing alignment..."

cd ${BLAST_ALIGNMENT_OUTPUTDIR}
#performing alignment
apptainer exec --bind /data/ ${BLAST_CONTAINER} blastn \
  -query ${FASTA_OUTPUTFILE} \
  -db ${BLAST_DB} \
  -outfmt "6 qstart qend sstart send length pident evalue bitscore qseqid sseqid" \
  -out ${BLAST_OUTPUTFILE}\
  -evalue 1e-10 \
  -word_size 11 \
  -reward 2 \
  -penalty -3 \
  -gapopen 5 \
  -gapextend 2 \
  -dust no \
  -soft_masking false \
  -num_threads ${SLURM_CPUS_PER_TASK}

if [ ! -f "${BLAST_OUTPUTFILE}" ]; then
    echo "ERROR: BLAST alignment failed"
    exit 1
fi

echo "blast self-alignment complete: ${BLAST_OUTPUTFILE}"


#PLOTTING
Rscript "${RSCRIPT}" \
    "${FASTA_OUTPUTFILE}" \
    "${BLAST_OUTPUTFILE}" \
    "${FLANKING_FILE}" \
    "${ACTUAL_GENE_FILE}" \
    "${geneID}" \
    "${geneName}" \
    "${R_OUTDIR}"


# fasta_file <- "5gene_window.fasta"
# blast_file <- "self_alignment_results/self_alignment.tsv"
# gene_info_file <- "flanking_list_for_VARIA_0007106.tsv"  # Your gene coordinates
# focal_gene_info <- "actual_gene_file.tsv"  # Focal gene info
# output_dir <- "self_alignment_results"

#extract fasta seq of gene window
#run blast on the sequence
#plot result
#parallelize?