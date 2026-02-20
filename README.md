Files needed
Data
    -annotation feature map
    -all gene names 
    -bam file
    fasta file
    -gff file
    -sample list of proteinIDs of genes of interest in a file: geneName proteinID (e.g FLC Q9S7Q7)
    -sample list of the gene coordinates of the genes of interest: contig geneStart geneEnd gene geneID (e.g: Bv4 65406076 65408129 SEP3 VARIA_0009872)
Degeneracy info
    -bedfile of degeneracy at all sites

Step 1 - annotation Pre-processing
This step prepares genome annotation files for downstream population genomic analyses

1_filterGFF.sh
Purpose : Filters a full MAKER-generated GFF file to retain only entries annotated as "gene" in column 3.
Input -> gff file
Output -> filtered gff file

2_filterForGenes.sh
Purpose: Filters an annotation feature map to extract entries corresponding to specific candidate genes.
This script runs as a SLURM array job. Each array task reads one line from a two-colum sample list and extracts the gene and protein_id, then filters the annot feat map for that protein ID
Input -> annotation feature map; proteinIDs and genes sample list
Output -> filteredForGeneProteinID.tsv per protein
