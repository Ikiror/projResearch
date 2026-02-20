#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=4
#SBATCH --job-name=blastDB
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/9_blast_db/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/9_blast_db/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#job to generate blast database 
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

echo "loading files..."

#file paths
WORKDIR="/data/users/aikiror/researchProject/pipeline"

#fasta file
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"

#output directories
OUTPUTDIR="${WORKDIR}/step1_prep/outputFiles/9_blast_db"

#blast output directories
BLAST_OUTPUTDIR="${OUTPUTDIR}/blast_output" #location of /windowsdb
mkdir -p $BLAST_OUTPUTDIR

#blast container
BLAST_CONTAINER="/containers/apptainer/blast_2.15.0.sif"

#BLAST - blastn
cd $BLAST_OUTPUTDIR
echo "creating blast database..."
#create blast database
apptainer exec --bind /data/ ${BLAST_CONTAINER} \
  makeblastdb \
  -in ${FASTA} \
  -dbtype nucl \
  -out ${BLAST_OUTPUTDIR}/genome_db \
  -title "Varia_genome" \
  -parse_seqids

echo "database in ${BLAST_OUTPUTDIR}"