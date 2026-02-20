#!/usr/bin/env Rscript

################################################################################
# SELF-ALIGNMENT DOT PLOT FOR STRUCTURAL VALIDATION
# 5-gene window with focal gene shading
################################################################################

args <- commandArgs(trailingOnly = TRUE)

if(length(args) < 7) {
  stop("Usage: Rscript script.R <fasta_file> <blast_file> <gene_info_file> <focal_gene_info> <geneID> <geneName> <output_dir>")
}

##############################################
# ARGUMENTS
##############################################

fasta_file <- args[1]           # Extracted FASTA sequence
blast_file <- args[2]           # BLAST alignment results
gene_info_file <- args[3]       # Flanking list
focal_gene_info <- args[4]      # Actual gene file
focal_id <- args[5]             # e.g., "VARIA_0007106"
focal_gene_name <- args[6]      # e.g., "FUL"
output_dir <- args[7]           # Output directory

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message("=================================================================")
message("BLAST SELF-ALIGNMENT VISUALIZATION")
message("=================================================================")
message("Gene ID: ", focal_id)
message("Gene name: ", focal_gene_name)
message("FASTA: ", fasta_file)
message("BLAST: ", blast_file)
message("Output: ", output_dir)
message("=================================================================")

##############################################
# LOAD REQUIRED LIBRARIES
##############################################

library(seqinr)


##############################################
# PARAMETERS
##############################################

# Filtering thresholds
min_identity <- 90  # Minimum percent identity to plot
min_length <- 100    # Minimum alignment length (bp)

# Color scheme
col_forward <- "blue"
col_reverse <- "red"
col_focal <- rgb(1, 0.8, 0.8, 0.3)  # Light red shading for focal gene

##############################################
# EXTRACT GENOMIC COORDINATES FROM FASTA HEADER
##############################################

message("\nExtracting genomic coordinates from FASTA header...")

# Read FASTA file
fasta <- read.fasta(fasta_file)

# Get the header (first sequence)
fasta_header <- names(fasta)[1]
message("FASTA header: ", fasta_header)

# Parse header to extract genomic position
# Expected format: >seqName:start-end or >seqName:start-end(+/-)
# Example: >Bv4:8726467-8846758 or >scaffold_31:765477-820426

# Extract contig and coordinates
if(grepl(":", fasta_header)) {
  # Split by colon
  parts <- strsplit(fasta_header, ":")[[1]]
  window_contig <- parts[1]
  
  # Extract start-end (remove any trailing characters like +/- in parentheses)
  # Keep the dash between coordinates!
  coords <- parts[2]
  # Remove parentheses and their contents: (anything)
  coords <- gsub("\\(.*\\)", "", coords)
  # Now split by dash
  coord_parts <- strsplit(coords, "-")[[1]]
  
  if(length(coord_parts) != 2) {
    stop("ERROR: Cannot parse coordinates from header: ", fasta_header)
  }

  window_start <- as.numeric(coord_parts[1])
  window_end <- as.numeric(coord_parts[2])
  
  # Validate
  if(is.na(window_start) || is.na(window_end)) {
    stop("ERROR: Could not convert coordinates to numbers. Header: ", fasta_header)
  }
  
} else {
  stop("ERROR: Cannot parse FASTA header. Expected format: >contig:start-end\nGot: ", fasta_header)
}

message("Window contig: ", window_contig)
message("Window start: ", window_start)
message("Window end: ", window_end)

# Get sequence length
seq_length <- length(fasta[[1]])
message("Sequence length: ", format(seq_length, big.mark = ","), " bp")

# Verify that coordinates match
expected_length <- window_end - window_start + 1
message("Expected length from coordinates: ", format(expected_length, big.mark = ","), " bp")

if(abs(seq_length - expected_length) > 10) {  # Allow small differences due to extraction
  warning("Sequence length (", seq_length, ") doesn't match coordinates (", expected_length, ")")
  message("Difference: ", abs(seq_length - expected_length), " bp")
} else {
  message("Sequence length matches coordinates")
}

##############################################
# LOAD GENE COORDINATES
##############################################

message("\nLoading gene coordinates...")

gene_coord_colnames <- c("contig", "gene_start", "gene_end", "gene", "varia_id")

# Load flanking genes
genes <- read.table(gene_info_file, header=FALSE, sep="\t", 
                   stringsAsFactors=FALSE, col.names=gene_coord_colnames)

# Load focal gene
focal_gene <- read.table(focal_gene_info, header=FALSE, sep="\t",
                        stringsAsFactors=FALSE, col.names=gene_coord_colnames)

# Get focal gene boundaries (absolute genomic coordinates)
focal_gene_row <- focal_gene[focal_gene$varia_id == focal_id, ]
if(nrow(focal_gene_row) == 0) {
  stop("ERROR: Focal gene ", focal_id, " not found in ", focal_gene_info)
}

focal_start_genomic <- focal_gene_row$gene_start
focal_end_genomic <- focal_gene_row$gene_end

message("Focal gene: ", focal_gene_name, " (", focal_id, ")")
message("Genomic position: ", focal_start_genomic, "-", focal_end_genomic)

# Filter genes to only those in this window
genes_in_window <- genes[
  genes$contig == window_contig &
  genes$gene_start >= window_start &
  genes$gene_end <= window_end,
]

message("Found ", nrow(genes_in_window), " genes in window")

##############################################
# LOAD BLAST RESULTS
##############################################

message("\nLoading BLAST results...")

blast_cols <- c("qstart", "qend", "sstart", "send", "length", "pident", 
                "evalue", "bitscore", "qseqid", "sseqid")
blast <- read.table(blast_file, header=FALSE, sep="\t", 
                   stringsAsFactors=FALSE, col.names=blast_cols)

message("Loaded ", nrow(blast), " alignments")

##############################################
# CONVERT BLAST COORDINATES TO GENOMIC COORDINATES
##############################################

message("\nConverting BLAST coordinates to genomic positions...")

# BLAST coordinates are relative to the extracted sequence (1-based)
# We need to convert them to absolute genomic coordinates

# For query coordinates (qstart, qend):
blast$qstart_genomic <- window_start + blast$qstart - 1
blast$qend_genomic <- window_start + blast$qend - 1

# For subject coordinates (sstart, send):
# Subject is also the same window, so same conversion
blast$sstart_genomic <- window_start + blast$sstart - 1
blast$send_genomic <- window_start + blast$send - 1

message("Coordinate conversion complete")

##############################################
# FILTER ALIGNMENTS
##############################################

message("\nFiltering alignments...")

# Remove self-alignments (perfect diagonal)
blast_filtered <- blast[
  blast$qstart != blast$sstart | blast$qend != blast$send,
]

# Apply identity and length filters
blast_filtered <- blast_filtered[
  blast_filtered$pident >= min_identity & 
  blast_filtered$length >= min_length,
]

# Determine strand (forward or reverse)
blast_filtered$strand <- ifelse(
  blast_filtered$sstart < blast_filtered$send, 
  "forward", 
  "reverse"
)

message("After filtering: ", nrow(blast_filtered), " alignments")
message("  Forward: ", sum(blast_filtered$strand == "forward"))
message("  Reverse: ", sum(blast_filtered$strand == "reverse"))

##############################################
# CREATE DOT PLOT
##############################################

message("\nCreating dot plot...")

png_file <- file.path(output_dir, paste0("self_alignment_dotplot_", focal_id, "_", focal_gene_name, ".png"))

png(png_file, width = 1800, height = 1800, res = 150, type = "cairo")

# Set up plot with genomic coordinates
par(mar = c(5.5, 10, 4.5, 2))
#par(mar = c(5.5, 12, 4.5, 2), mgp = c(3.5, 1, 0))

par(mar = c(5.5, 10, 4.5, 2))

plot(0, 0, type = "n",
     xlim = c(window_start, window_end),
     ylim = c(window_start, window_end),
     xlab = "Genomic Position (bp)",
     ylab = "",
     main = paste0("Self-alignment: 5-gene window\nFocal gene: ", 
                  focal_gene_name, " (", focal_id, ")"),
     cex.lab = 1.4,
     cex.main = 1.5,
     cex.axis = 1.2,
     xaxt = "n",
     yaxt = "n",
     asp = 1)  # Keep square

# Custom axes with proper formatting
x_ticks <- seq(window_start, window_end, length.out = 5)
y_ticks <- seq(window_start, window_end, length.out = 5)

axis(1, at = x_ticks, labels = format(x_ticks, scientific = FALSE, big.mark = ","), 
     cex.axis = 1.2)
axis(2, at = y_ticks, labels = format(y_ticks, scientific = FALSE, big.mark = ","), 
     cex.axis = 1.2, las = 1)

# Add y-axis label with proper spacing
mtext("Genomic Position (bp)", side = 2, line = 7.5, cex = 1.4)

# Add grid
grid(col = "gray90", lty = 1)

# Shade focal gene region (both axes)
rect(focal_start_genomic, window_start, focal_end_genomic, window_end, 
     col = col_focal, border = "darkred", lwd = 2)
rect(window_start, focal_start_genomic, window_end, focal_end_genomic, 
     col = col_focal, border = "darkred", lwd = 2)

# Add focal gene label
text(focal_start_genomic + (focal_end_genomic - focal_start_genomic)/2, 
     window_end * 0.98,
     labels = "Focal Gene", col = "darkred", cex = 1.2, font = 2)

# Plot alignments
# Forward strand (blue)
forward <- blast_filtered[blast_filtered$strand == "forward", ]
if(nrow(forward) > 0) {
  segments(forward$qstart_genomic, forward$sstart_genomic, 
           forward$qend_genomic, forward$send_genomic,
           col = adjustcolor(col_forward, alpha.f = 0.4), 
           lwd = 2)  # *** INCREASED FROM 0.5 to 2 ***
}

# Reverse strand (red)
reverse <- blast_filtered[blast_filtered$strand == "reverse", ]
if(nrow(reverse) > 0) {
  segments(reverse$qstart_genomic, reverse$sstart_genomic, 
           reverse$qend_genomic, reverse$send_genomic,
           col = adjustcolor(col_reverse, alpha.f = 0.4), 
           lwd = 2)  # *** INCREASED FROM 0.5 to 2 ***
}

# Plot main diagonal (self-alignment)
abline(a = window_start - window_start, b = 1, col = "black", lwd = 2.5, lty = 2)

# Add gene boundaries as vertical/horizontal lines
if(nrow(genes_in_window) > 0) {
  for(i in 1:nrow(genes_in_window)) {
    abline(v = genes_in_window$gene_start[i], col = "gray60", lty = 3, lwd = 0.8)
    abline(v = genes_in_window$gene_end[i], col = "gray60", lty = 3, lwd = 0.8)
    abline(h = genes_in_window$gene_start[i], col = "gray60", lty = 3, lwd = 0.8)
    abline(h = genes_in_window$gene_end[i], col = "gray60", lty = 3, lwd = 0.8)
  }
}

# Add gene annotations at bottom and left
# if(nrow(genes_in_window) > 0) {
#   for(i in 1:nrow(genes_in_window)) {
#     gene_mid <- (genes_in_window$gene_start[i] + genes_in_window$gene_end[i]) / 2
    
#     # Gene name color
#     gene_col <- if(genes_in_window$varia_id[i] == focal_id) "darkred" else "black"
#     gene_font <- if(genes_in_window$varia_id[i] == focal_id) 2 else 1
    
#     # Bottom annotation
#     text(gene_mid, window_start - (window_end - window_start) * 0.02, 
#          labels = genes_in_window$gene[i], 
#          srt = 45, adj = 1, cex = 0.9, col = gene_col, font = gene_font, xpd = TRUE)
    
#     # Left annotation
#     text(window_start - (window_end - window_start) * 0.02, gene_mid,
#          labels = genes_in_window$gene[i], 
#          srt = 45, adj = 1, cex = 0.9, col = gene_col, font = gene_font, xpd = TRUE)
#   }
# }

# Enhanced legend
legend("topleft",
       legend = c(paste0("Forward (", sum(blast_filtered$strand == "forward"), ")"),
                 paste0("Reverse (", sum(blast_filtered$strand == "reverse"), ")"),
                 "Focal gene region",
                 "Perfect self-match"),
       col = c(col_forward, col_reverse, "darkred", "black"),
       lty = c(1, 1, NA, 2),
       lwd = c(2, 2, NA, 2.5),
       fill = c(NA, NA, col_focal, NA),
       border = c(NA, NA, "darkred", NA),
       cex = 1.2,
       bg = "white",
       box.lwd = 1.5)

# Statistics box
legend("bottomright",
       legend = c(
         paste0("Alignments: ", nrow(blast_filtered)),
         paste0("Min identity: ", min_identity, "%"),
         paste0("Min length: ", min_length, " bp"),
         paste0("Window: ", format(window_start, big.mark = ","), "-", 
                format(window_end, big.mark = ","))
       ),
       cex = 1.0,
       bg = "white",
       box.lwd = 1.5)

dev.off()

message("\nPlot saved: ", png_file)

##############################################
# SUMMARY STATISTICS
##############################################

message("\n", strrep("=", 70))
message("SUMMARY STATISTICS")
message(strrep("=", 70))

# Count alignments in focal gene region
focal_alignments <- blast_filtered[
  (blast_filtered$qstart_genomic >= focal_start_genomic & 
   blast_filtered$qstart_genomic <= focal_end_genomic) |
  (blast_filtered$qend_genomic >= focal_start_genomic & 
   blast_filtered$qend_genomic <= focal_end_genomic) |
  (blast_filtered$sstart_genomic >= focal_start_genomic & 
   blast_filtered$sstart_genomic <= focal_end_genomic) |
  (blast_filtered$send_genomic >= focal_start_genomic & 
   blast_filtered$send_genomic <= focal_end_genomic),
]

message("Total alignments: ", nrow(blast_filtered))
message("  In/near focal gene: ", nrow(focal_alignments))
message("  Outside focal gene: ", nrow(blast_filtered) - nrow(focal_alignments))
message("\nMean identity: ", round(mean(blast_filtered$pident), 2), "%")
message("Mean length: ", round(mean(blast_filtered$length), 0), " bp")

# Detect potential duplications
high_identity <- blast_filtered[blast_filtered$pident > 98, ]
if(nrow(high_identity) > 0) {
  message("\nWARNING: ", nrow(high_identity), " high-identity alignments (>98%) detected!")
  message("This may indicate:")
  message("  - Tandem duplications")
  message("  - Segmental duplications")
  message("  - Assembly artifacts")
}

# Check for inversions in focal gene
focal_inversions <- focal_alignments[focal_alignments$strand == "reverse", ]
if(nrow(focal_inversions) > 0) {
  message("\nWARNING: ", nrow(focal_inversions), " reverse alignments in focal gene region!")
  message("This may indicate inversions or misassemblies.")
}

message("\n", strrep("=", 70))
message("Analysis complete!")
message(strrep("=", 70))