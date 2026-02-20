#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --mem=512gb
#SBATCH --cpus-per-task=32
#SBATCH --job-name=VCF0foldExtract
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/05_VCF_0fold_extract/output_%J.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step3_popGenStats/logFiles/05_VCF_0fold_extract/error_%J.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#filters vcf file by 0-fold degenerate bed file
#alter SBATCH --output,error, workdir to match file paths
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,

# Define variables
WORKDIR="/data/users/aikiror/researchProject/pipeline"
OUTPUTDIR="${WORKDIR}/step3_popGenStats/outputFiles/05_VCF_0fold_extract"

#container
CONTAINER="/containers/apptainer/gatk4_4.6.2.0_samtools_1.21.sif"

#0fold sites from .bed file
BED0FOLD="${WORKDIR}/step3_popGenStats/outputFiles/02_extract_0fold/0fold_degenerate_sites.bed"

#fasta file
FASTA="${WORKDIR}/biscutellaVariaData/Varia_review2.FINAL.Syri.fa"

VCF="${WORKDIR}/step3_popGenStats/outputFiles/01_VCFFileFiltering/missingness_filtered.vcf.gz"
OUTFILENAME="${OUTPUTDIR}/0_fold_degenerate_filtered.vcf.gz"

#make outputdir if it doesnt already exist
mkdir -p $OUTPUTDIR

#-R => fasta; -O => output path; -V => VCF file w/ variants; -L => genomic interval to operate over; 
#--select-type => select only a certain type of variants from the input file to include
#--restrict-alleles-to => select only variants of a particular allelicity

#filter vcf file by 0-fold degenerate bed file
apptainer exec --bind /data/ ${CONTAINER} gatk SelectVariants -R ${FASTA} \
    -V ${VCF} \
    -O ${OUTFILENAME} \
    -L ${BED0FOLD} \
    --select-type SNP \
    --restrict-alleles-to BIALLELIC


 