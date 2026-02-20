#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --array=1-20
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=sortGFF
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/6_filter_sorted_gff_for_flanking_genes/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/6_filter_sorted_gff_for_flanking_genes/error_%A_%a.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#filters sorted gff file for varia_id gene with 2 left flanking genes, query gene, and 2 right flanking genes
#returns new file with those 5 varia_id genes
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of proteins
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

# Input files
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTDIR="$WORKDIR/step1_prep/outputFiles/6_filter_sorted_gff_for_flanking_genes"

#gff file
GFF="$WORKDIR/step1_prep/outputFiles/2_sortFilteredGFF/sorted_filteredOutput.gff"

#actual gene coords
SAMPLE="$WORKDIR/biscutellaVariaData/geneCoord.tsv"

# OUT="flanks_2Genes.tsv"

mkdir -p $OUTDIR

cd $OUTDIR

#extract gene info based off of current task
seqName=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLE`
start=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLE`
stop=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLE`
geneName=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $4; exit}' $SAMPLE`
variaID=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $5; exit}' $SAMPLE`

#define individual file for each variaID
OUTFILE=$"filteredGFFfor_${variaID}.gff"

#get the line number of the gene of interest
line_number=$(grep -n "$variaID" ${GFF} | cut -d: -f1)

#get line numbers for the flanking genes
left_flank_1stgene=$((line_number-2))
left_flank_2ndgene=$((line_number-1))
right_flank_1stgene=$((line_number+1))
right_flank_2ndgene=$((line_number+2))

#extract info from GFF and add to respective file
sed -n "${left_flank_1stgene}p" $GFF >> $OUTFILE
sed -n "${left_flank_2ndgene}p" $GFF >> $OUTFILE
sed -n "${line_number}p" $GFF >> $OUTFILE
sed -n "${right_flank_1stgene}p" $GFF >> $OUTFILE
sed -n "${right_flank_2ndgene}p" $GFF >> $OUTFILE