#!/bin/bash
#SBATCH --array=1-20
#SBATCH --time=00:20:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sampleListGenerator
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/8_flanking_gene_coord_list/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/8_flanking_gene_coord_list/error_%A_%a.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#generates a gene coord list of the 5-gene window
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of genes
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

#workdir
WORKDIR="/data/users/aikiror/researchProject/pipeline"

# Define the output folder and file where the list will be saved
OUTPUTDIR="${WORKDIR}/biscutellaVariaData"
mkdir -p ${OUTPUTDIR}

OUTPUT_FILE="${OUTPUTDIR}/flanking_2genes_gene_coords.tsv"

cd $OUTPUTDIR

SAMPLELIST="${WORKDIR}/step1_prep/outputFiles/7_flanking_genes_gff/gene_w_flanks_sample_list.gff"
GENE_INFO="${WORKDIR}/biscutellaVariaData/geneCoord.tsv"

#extract info directly from columns
scaffold=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
start=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $4; exit}' $SAMPLELIST`
stop=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $5; exit}' $SAMPLELIST`

#extract varia_id from .gff formating
varia_id=$(awk -v line="$SLURM_ARRAY_TASK_ID" '
    NR == line {
        if (match($9, /ID=([^;]+)/, arr)) {
            print arr[1];
        }
        exit
    }' "$SAMPLELIST")
#find geneName based off of current varia_id
geneName=$(grep -P "$varia_id" $GENE_INFO | awk '{print $4}')

echo -e "${scaffold}\t${start}\t${stop}\t${geneName}\t${varia_id}" >> ${OUTPUT_FILE}
