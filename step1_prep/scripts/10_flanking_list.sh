#!/bin/bash
#SBATCH --array=1-20
#SBATCH --time=00:40:00
#SBATCH --mem=1g
#SBATCH --cpus-per-task=1
#SBATCH --job-name=filter5genewindowGFF
#SBATCH --output=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/10_flanking_list/output_%A_%a.o
#SBATCH --error=/data/users/aikiror/researchProject/pipeline/step1_prep/logFiles/10_flanking_list/error_%A_%a.e
#SBATCH --partition=pibu_el8

#to run: sbatch path/to/file
#filters the 5-gene window gff for contig geneStart geneStop variaID per gene
#alter SBATCH --output,error, workdir to match file paths
#alter SBATCH --array to match number of proteins
#option to add #SBATCH --mail-user=youremailaddress, SBATCH --mail-type=begin,end,fail

WORKDIR="/data/users/aikiror/researchProject/pipeline"
# Define the output folder and file where the list will be saved

OUTPUTDIR="${WORKDIR}/step1_prep/outputFiles/10_flanking_list"
mkdir -p ${OUTPUTDIR}

FOCAL_GENE_LIST="/biscutellaVariaData/geneCoord.tsv"
FLANKING_GFF_FOLDER="$WORKDIR/step1_prep/outputFiles/6_filter_sorted_gff_for_flanking_genes"
FLANKING_GFF_FILE=$(ls ${FLANKING_GFF_FOLDER}/filteredGFFfor_VARIA_*.gff | sed -n "${SLURM_ARRAY_TASK_ID}p")

#get focal varia_id per file
base=$(basename ${FLANKING_GFF_FILE})
focal_id="${base#filteredGFFfor_}"
focal_id="${focal_id%.gff}"

OUTPUT_FILE="${OUTPUTDIR}/flanking_list_for_${focal_id}.tsv"
MASTER_FILE="${OUTPUTDIR}/master_flanking_list.tsv"
cd $OUTPUTDIR

#match focal varia id to focal gene info in actual gene list - will be a placeholder for flanking varia_ids 
geneName=$(grep -P "$focal_id" $FOCAL_GENE_LIST | awk '{print $4}')

LINES=$(wc -l < "$FLANKING_GFF_FILE") #number of lines to iterate over

for ((line=1; line<=LINES; line++)); do
    scaffold=$(awk -v line="$line" 'NR==line{print $1; exit}' "$FLANKING_GFF_FILE")
    start_coord=$(awk -v line="$line" 'NR==line{print $4; exit}' "$FLANKING_GFF_FILE")
    stop_coord=$(awk -v line="$line" 'NR==line{print $5; exit}' "$FLANKING_GFF_FILE")
    varia_id=$(awk -v line="$line" '
        NR == line {
            if (match($9, /ID=([^;]+)/, arr)) {
                print arr[1];
            }
            exit
        }' "$FLANKING_GFF_FILE")

    #append to same file
    echo -e "${scaffold}\t${start_coord}\t${stop_coord}\t${geneName}\t${varia_id}" >> ${OUTPUT_FILE}
    echo -e "${scaffold}\t${start_coord}\t${stop_coord}\t${geneName}\t${varia_id}" >> ${MASTER_FILE}

    echo "Done: written info from $FLANKING_GFF_FILE to $OUTPUT_FILE and $MASTER_FILE"
done

