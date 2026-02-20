#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=coords
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/7_flanking_genes_gff/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/7_flanking_genes_gff/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#generates gff file of 5-gene window
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

# Input files
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTDIR="$WORKDIR/step1_prep/outputFiles/7_flanking_genes_gff"
COORD_INFO_FOLDER="$WORKDIR/step1_prep/outputFiles/6_filter_sorted_gff_for_flanking_genes"

#define output file 
OUTFILE="gene_w_flanks_sample_list.gff"

mkdir -p $OUTDIR

cd $OUTDIR

for coordFile in $COORD_INFO_FOLDER/filteredGFFfor_VARIA*; do
    scaffold=`awk -v line=3 'NR==line{print $1; exit}' $coordFile`
    source=`awk -v line=3 'NR==line{print $2; exit}' $coordFile`
    feature=`awk -v line=3 'NR==line{print $3; exit}' $coordFile`
    start_coord=`awk -v line=1 'NR==line{print $4; exit}' $coordFile`
    stop_coord=`awk -v line=5 'NR==line{print $5; exit}' $coordFile`
    score=`awk -v line=3 'NR==line{print $6; exit}' $coordFile`
    strand=`awk -v line=3 'NR==line{print $7; exit}' $coordFile`
    phase=`awk -v line=3 'NR==line{print $8; exit}' $coordFile`
    id_info=`awk -v line=3 'NR==line{print $9; exit}' $coordFile`

    #append to same file
    echo -e "${scaffold}\t${source}\t${feature}\t${start_coord}\t${stop_coord}\t${score}\t${strand}\t${phase}\t${id_info}" >> $OUTFILE

    echo "Done: written info from $coordFile to $OUTFILE"
done

echo "Done: written to $OUTFILE"

#copy file into additional directory
cp $OUTFILE $SECONDARY_OUTDIR/$OUTFILE