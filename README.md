# Biscutella laevigata Population Genomics Pipeline

Population genomics analysis of flowering-time genes across elevational gradients in alpine *Biscutella laevigata*, comparing 0-fold vs 4-fold degenerate sites to detect signatures of natural selection.

## Quick Start

```bash
# Clone repository
git clone https://github.com/Ikiror/projResearch
cd biscutella-popgen

This project uses two tools: gawk v5.0.0 and piawka v0.8.11 (https://github.com/novikovalab/piawka)

Changes have been made to piawka that deviate slightly form the orginial to better handle division by 0.

Original: (piawka/scripts/piawka)

"function printOutput( locus, nSites, pop1, pop2, nUsed, metric, numerator, denominator, nGeno, nMiss ) {
  out=locus"\t"nSites"\t"pop1"\t"pop2"\t"nUsed"\t"metric"\t"numerator/denominator"\t"numerator"\t"denominator"\t"nGeno"\t"nMiss
  say("", 1) # overwrite error messages
  if ( args["jobs"] > 1 ) {
    print out > outbuffer
  } else {
    print out
  }
}" 

Changed:
function printOutput( locus, nSites, pop1, pop2, nUsed, metric, numerator, denominator, nGeno, nMiss ) {
  value = (denominator == 0 ? "nan" : numerator / denominator)
  out = locus"\t"nSites"\t"pop1"\t"pop2"\t"nUsed"\t"metric"\t"value"\t"numerator"\t"denominator"\t"nGeno"\t"nMiss
  #out=locus"\t"nSites"\t"pop1"\t"pop2"\t"nUsed"\t"metric"\t"numerator/denominator"\t"numerator"\t"denominator"\t"nGeno"\t"nMiss
  say("", 1) # overwrite error messages
  if ( args["jobs"] > 1 ) {
    print out > outbuffer
  } else {
    print out
  }
}


pipeline/
├── step1_prep/                          # Data preparation pipeline
│   └── scripts/
│       ├── 1_filterGFF.sh               # Filter GFF for gene features
│       ├── 2_sortFilteredGFF.sh         # Sort by contig and coordinates
│       ├── 3_filterForGenes.sh          # Extract genes by protein ID
│       ├── 4_actual_gene_coord_list.sh  # Generate gene coordinate list
│       ├── 5_fasta_indexing.sh          # Index reference genome
│       ├── 6_filter_sorted_gff_for_flanking_genes.sh  # Extract 5-gene windows
│       ├── 7_flanking_genes_gff.sh      # Create flanking gene GFF
│       ├── 8_flanking_gene_coord_list.sh  # Generate flanking coordinates
│       ├── 9_blast_db.sh                # Build BLAST database
│       ├── 10_flanking_list.sh          # Create per-gene flanking lists
│       ├── 11_bedtools_blast_alignment_plot.sh  # BLAST alignment & plot
│       ├── 11_blast_alignment_plot.R    # R script for alignment visualization
│       └── 11.5_single_run_blast_alignment.sh   # Single-gene BLAST run
│
├── step2_calcCoverageDepth/             # Coverage analysis
│   └── scripts/
│       ├── 0.5_bamIndexing.sh           # Index BAM file
│       └── 1_samtools_depth_per_region.sh  # Calculate depth per region
│
├── step3_popGenStats/                   # Population genetics statistics
│   └── scripts/
│       ├── 01_VCFFileFiltering.sh       # Filter VCF files
│       ├── 01pt5_VCF_indexing.sh        # Index VCF files
│       ├── 02_extract_0fold.sh          # Extract 0-fold sites
│       ├── 03_extract_4fold.sh          # Extract 4-fold sites
│       ├── 04_dictFileCreation.sh       # Create sequence dictionary
│       ├── 05_VCF_0fold_extract.sh      # VCF extraction (0-fold)
│       ├── 06_VCF_4fold_extract.sh      # VCF extraction (4-fold)
│       ├── 07_generate0FoldBedByContig.sh   # Generate 0-fold BED by contig
│       ├── 08_VCF_0fold_extract_by_contig.sh  # Extract 0-fold by contig
│       ├── 08.5_individual_VCF_0fold_extract_by_contig.sh  # Individual 0-fold
│       ├── 09_generate4FoldBedByContig.sh   # Generate 4-fold BED by contig
│       ├── 10_VCF_4fold_extract_by_contig.sh  # Extract 4-fold by contig
│       ├── 10.5_individual_VCF_4fold_extract_by_contig.sh  # Individual 4-fold
│       ├── 11_combineVCF.sh             # Combine VCF files
│       ├── 11.5_alt_combine_bed_file.sh  # Alternative BED combine
│       ├── 12_piawka_prep.sh            # Prepare for pixy analysis
│       ├── 14.95_piawka_stats_tajimasD_full_bedfile_0fold.sh  # Pixy stats (0-fold)
│       └── 14.97_piawka_stats_tajimasD_full_bedfile_4fold.sh  # Pixy stats (4-fold)
│
├── step4_plotting/                      # Visualization & analysis
│   └── scripts/
│       ├── 1_piawka_R_visualization_w_boxplots.R  # Main R analysis script
│       ├── 1_single_run_1_piawka_R_visualization_w_boxplots.sh  # Single gene run
│       └── 1_multirun_piawka_R_visualization_w_boxplots.sh  # Batch processing
│
├── biscutellaVariaData/                 # Input data directory
│   ├── Varia_review2.FINAL.Syri.fa      # Reference genome
│   ├── Varia.all.maker.renamed.gff      # Gene annotations
│   ├── annot_feature_map_VARIA.txt      # Annotation feature map
│   ├── ProteinID_sampleList.tsv         # Gene-to-protein ID mapping
│   ├── geneCoord.tsv                    # Gene coordinates
│   ├── flanking_2genes_gene_coords.tsv  # 5-gene window coordinates
│   ├── master_flanking_list.tsv         # All gene windows
│   └── VARI_short_reads.bam             # Sequencing reads
│
├── degeneracyInfo/                      # Degeneracy classification data
├── vcf_file/                            # VCF files (variant calls)
├── tools/                               # Additional tools/utilities
└── README.md


Pipeline Workflow
Step 1: Data Preparation (step1_prep/)
Purpose: Extract gene coordinates, build reference databases, identify gene windows

Filter & Sort GFF (1-2_*.sh)

Extract gene features from annotation
Sort by scaffold and genomic position


Gene Identification (3-4_*.sh)

Match genes to protein IDs
Generate coordinate lists


Reference Setup (5_*.sh, 9_*.sh)

Index FASTA file
Build BLAST database


Gene Window Definition (6-8_*.sh, 10_*.sh)

Extract focal gene ±2 flanking genes (5-gene windows)
Create coordinate lists for local genomic context


Sequence Alignment (11_*.sh, 11_*.R)

BLAST alignment within gene windows
Visualization of genomic context



Step 2: Coverage Analysis (step2_calcCoverageDepth/)
Purpose: Calculate sequencing depth across gene regions

0.5: Index BAM files
1: Calculate depth per region using samtools

Step 3: Population Genetics Statistics (step3_popGenStats/)
Purpose: Extract variants, classify degeneracy, calculate statistics
Workflow:

VCF Preparation (01-01pt5_*.sh)

Filter and index VCF files


Degeneracy Classification (02-03_*.sh)

Extract 0-fold (nonsynonymous) sites
Extract 4-fold (synonymous) sites


Reference Files (04_*.sh)

Create sequence dictionary


VCF Extraction (05-06_*.sh)

Extract variants at 0-fold sites
Extract variants at 4-fold sites


Contig-Level Processing (07-10.5_*.sh)

Generate BED files by contig
Extract VCF data by contig for both degeneracy classes


File Consolidation (11-11.5_*.sh)

Combine VCF files
Alternative BED file combination


Pixy Preparation & Statistics (12_*.sh, 14.95-14.97_*.sh)

Prepare files for pixy analysis
Calculate Fst, π, and Tajima's D for 0-fold and 4-fold sites



Step 4: Visualization & Analysis (step4_plotting/)
Purpose: Generate comprehensive figures and calculate π₀/π₄ ratios
Main Script: 1_piawka_R_visualization_w_boxplots.R
Execution Scripts:

1_single_run_1_piawka_R_visualization_w_boxplots.sh - Single gene
1_multirun_piawka_R_visualization_w_boxplots.sh - Batch processing

Metrics Calculated:

Fst: Population differentiation (low vs high elevation)
π (nucleotide diversity): Genetic variation within populations
Tajima's D: Test for selection (±2 thresholds)
π₀/π₄ ratio: Functional vs neutral diversity

Output: 12-panel comparison figures + CSV tables

Input Files (description and format)
Required Input Data (in biscutellaVariaData/)
FileDescriptionFormat
Reference genome -> FASTA
Gene annotations -> GFF3
Annotation feature map -> TSV
Gene-to-proteinID mapping -> TSV: gene, uniProtID
Sequencing reads (aligned) -> BAM


-Ikiror
-UniBe
