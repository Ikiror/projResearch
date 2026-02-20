#!/usr/bin/env Rscript

################################################################################
# COMPLETE POPULATION GENOMICS ANALYSIS
# 0-fold vs 4-fold degenerate sites with Coverage
# Elevation adaptation study
################################################################################

##############################################
# CLUSTER-FRIENDLY SETUP - !!!!! if running locally, type=cairo may cause errors
##############################################
#options(bitmapType='cairo')
#options(warn=-1)

args <- commandArgs(trailingOnly = TRUE)
##############################################
# FILE PATHS
##############################################
message("Loading data files...")
# 0-fold data
piawka_0fold_stats <- args[1] #piawka stats for the 0fold deg sites
combined_bed_0fold <- args[2] #bedfile of the 0fold deg sites

# 4-fold data
piawka_4fold_stats <- args[3] #as above but for 4fold deg sites
combined_bed_4fold <- args[4]

# Gene info (shared)
focal_gene_info <- args[5] #info about the focal gene being looked into -> contig gene_start gene_end gene_name varia_id
gene_file_w_flanks <- args[6] #same as focal gene info; gene_start=window_start gene_end=window_end; window = +/-2genes flanking focal
selected_gene_of_interest_window <- args[7] #like focal gene info but w/ actual gene info of each gene in window rows: -2 -1 focal +1 +2
master_flanking_windows <- args[8] #master list of all selected_gene_of_interest_window's

# Coverage files directory (optional)
coverage_dir <- args[9]  # Directory containing coverage files (focal_id,"_", focal_contig, "_", focal_gene_name,"_coverage.tsv"); coverage files: contig, position, depth

# Output directory
outdir <- args[10]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

message("Data files loaded")

message(strrep("=", 70))
message("STARTING POPULATION GENOMICS ANALYSIS")
message(strrep("=", 70))

################################################################################
# FUNCTION: ADD PANEL LABEL
################################################################################

add_panel_label <- function(label, x = 0.98, y = 0.98, cex = 1.5, font = 4) {
  # Add panel label (A, B, C, etc.) to top-right of plot
  u <- par("usr")
  x_pos <- u[1] + (u[2] - u[1]) * x
  y_pos <- u[3] + (u[4] - u[3]) * y
  text(x_pos, y_pos, labels = label, cex = cex, font = font, adj = c(1, 1))
}

################################################################################
# FUNCTION: LOAD AND PROCESS COVERAGE DATA
################################################################################

load_coverage_data <- function(coverage_file, focal_gene_start, focal_gene_end, 
                               window_start, window_end) {
  
  message("Loading coverage data from: ", coverage_file)
  
  # Check if file exists
  if(!file.exists(coverage_file)) {
    message("WARNING: Coverage file not found: ", coverage_file)
    return(NULL)
  }
  
  # Load coverage file (contig, position, depth)
  coverage <- read.table(coverage_file, header=FALSE, sep="\t", 
                         stringsAsFactors=FALSE,
                         col.names=c("contig", "pos", "depth"))
  
  message("Loaded ", nrow(coverage), " positions with coverage data")
  
  # Filter for the geneic window region
  coverage_window <- coverage[
    coverage$pos >= window_start & 
      coverage$pos <= window_end,
  ]
  
  message("Filtered to ", nrow(coverage_window), " positions in window")
  
  # Calculate summary statistics
  mean_depth <- mean(coverage_window$depth, na.rm=TRUE)
  median_depth <- median(coverage_window$depth, na.rm=TRUE)
  sd_depth <- sd(coverage_window$depth, na.rm=TRUE)
  
  message("Mean coverage: ", round(mean_depth, 2), "x")
  message("Median coverage: ", round(median_depth, 2), "x")
  message("SD coverage: ", round(sd_depth, 2), "x")
  
  # Optional: Bin coverage for smoother visualization if data is dense
  # Only bin if we have more than 10,000 positions
  if(nrow(coverage_window) > 10000) {
    message("Binning coverage data for visualization (", nrow(coverage_window), " positions)...")
    bin_size <- 1000  # 1kb bins
    
    coverage_window$bin <- floor(coverage_window$pos / bin_size) * bin_size
    
    coverage_binned <- aggregate(
      depth ~ bin,
      data = coverage_window,
      FUN = function(x) c(mean = mean(x), min = min(x), max = max(x))
    )
    
    coverage_binned <- do.call(data.frame, coverage_binned)
    names(coverage_binned)[2:4] <- c("mean_depth", "min_depth", "max_depth")
    
    return(list(
      raw = coverage_window,
      binned = coverage_binned,
      mean_depth = mean_depth,
      median_depth = median_depth,
      sd_depth = sd_depth,
      use_binned = TRUE
    ))
  } else {
    # For smaller regions, use raw data
    return(list(
      raw = coverage_window,
      binned = NULL,
      mean_depth = mean_depth,
      median_depth = median_depth,
      sd_depth = sd_depth,
      use_binned = FALSE
    ))
  }
}

################################################################################
# FUNCTION: PLOT COVERAGE PANEL
################################################################################

plot_coverage_panel <- function(coverage_data, focal_gene_start, focal_gene_end,
                                focal_gene_name, focal_id) {
  
  if(is.null(coverage_data)) {
    # Create empty plot with message
    plot.new()
    text(0.5, 0.5, "Coverage data not available", cex = 1.5, col = "gray50")
    add_panel_label("A")
    return()
  }
  
  # Determine which data to plot
  if(coverage_data$use_binned && !is.null(coverage_data$binned)) {
    plot_data <- coverage_data$binned
    x_col <- "bin"
    y_col <- "mean_depth"
  } else {
    plot_data <- coverage_data$raw
    x_col <- "pos"
    y_col <- "depth"
  }
  
  # Calculate y-axis limit
  y_max <- max(plot_data[[y_col]], na.rm=TRUE) * 1.3
  
  # Create plot
  plot(plot_data[[x_col]], plot_data[[y_col]],
       type = "n",
       ylim = c(0, y_max),
       main = paste0("Sequencing Coverage ", focal_gene_name, " region: ", focal_id),
       ylab = "Coverage Depth (reads)",
       xlab = "Genomic Position",
       xaxt = "n",
       cex.main = 1.2,
       cex.lab = 1.1,
       cex.axis = 1.0)
  
  # Custom x-axis without scientific notation
  axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","))
  
  # Shade focal gene region
  rect(focal_gene_start, 0, focal_gene_end, y_max,
       col = rgb(1, 0.8, 0.8, 0.4), border = "darkred", lwd = 1.5)
  
  # Plot coverage as filled area
  if(coverage_data$use_binned && !is.null(coverage_data$binned)) {
    # For binned data, show mean with min/max range
    polygon(c(plot_data$bin, rev(plot_data$bin)),
            c(plot_data$min_depth, rev(plot_data$max_depth)),
            col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
    
    lines(plot_data$bin, plot_data$mean_depth, 
          col = "darkblue", lwd = 2)
  } else {
    # For raw data, plot as filled polygon
    polygon(c(plot_data[[x_col]], rev(plot_data[[x_col]])),
            c(plot_data[[y_col]], rep(0, nrow(plot_data))),
            col = rgb(0.2, 0.4, 0.8, 0.5), border = NA)
    
    lines(plot_data[[x_col]], plot_data[[y_col]], 
          col = "darkblue", lwd = 1.5)
  }
  
  # Add horizontal line for mean coverage
  abline(h = coverage_data$mean_depth, lty = 2, col = "red", lwd = 1.5)
  
  # Add text label for mean coverage - on plot
  text(min(plot_data[[x_col]]), coverage_data$mean_depth,
       labels = paste0("Mean: ", round(coverage_data$mean_depth, 1), "x"),
       pos = 3, cex = 1, col = "red", offset = 0.3)
  
  # Add median line
  abline(h = coverage_data$median_depth, lty = 2, col = "darkgreen", lwd = 1)
  
  # Add text label for median coverage - on plot
  text(max(plot_data[[x_col]]), coverage_data$median_depth,
       labels = paste0("Median: ", round(coverage_data$median_depth, 1), "x"),
       pos = 3, cex = 1, col = "darkgreen", offset = 0.3)
  
  # Add legend
  legend("topleft",
         legend = c(paste0("Coverage depth"), 
                    paste0("Mean (", round(coverage_data$mean_depth, 1), "x)"),
                    paste0("Median (", round(coverage_data$median_depth, 1), "x)")),
         col = c("darkblue", "red", "darkgreen"),
         lty = c(1, 2, 2),
         lwd = c(1.5, 1.5, 1),
         cex = 1,
         bg = "white")
  
  # Add panel label
  add_panel_label("A")
}

################################################################################
# FUNCTION DEFINITIONS
################################################################################

analyze_degeneracy <- function(piawka_file, bed_file, degeneracy_label, coverage_file = NULL) {
  
  message("\n", strrep("=", 70))
  message("ANALYZING ", degeneracy_label, " SITES")
  message(strrep("=", 70))
  
  ##############################################
  # LOAD DATA
  ##############################################
  
  message("Loading piawka data...")
  piawka <- read.table(piawka_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  message("Loading bed file...")
  bed <- read.table(
    bed_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("contig_bed", "start0", "end0", "locus", "degeneracy", 
                  "ref_base", "aa_ref", "alt_map")
  )
  bed$pos1 <- bed$start0 + 1 #bedfile is zero based, account for that when incorporating it with other data
  
  ##############################################
  # PROCESS PIAWKA DATA
  ##############################################
  
  message("Processing piawka data...")
  piawka$varia_id <- sub("-.*$", "", piawka$locus) #VARIA_XXXXX-RA:229 -> VARIA_XXXXXX
  
  relevant_varia_ids <- unique(flank_master_list$varia_id) #get unique varia ids from master list in case of overlap
  filtered_piawka <- piawka[piawka$varia_id %in% relevant_varia_ids, ] #filter piawka based off of master list
  
  gene_info_unique <- flank_master_list[, c("varia_id", "contig", "gene_start", "gene_end", "gene")]
  gene_info_unique <- gene_info_unique[!duplicated(gene_info_unique$varia_id), ] #remove duplicates
  
  annotated_filtered_piawka <- merge( #adding info on gene coord etc to piawka; join by varia id
    filtered_piawka,
    gene_info_unique,
    by = "varia_id",
    all.x = TRUE
  )
  
  bed_positions <- bed[, c("locus", "contig_bed", "pos1")] #adding bedfile info to piawka, join by locus
  annotated_filtered_piawka <- merge(
    annotated_filtered_piawka,
    bed_positions,
    by = "locus",
    all.x = TRUE
  )
  names(annotated_filtered_piawka)[names(annotated_filtered_piawka) == "contig_bed"] <- "contig_site" 
  
  selected_geneic_window_annotated_filtered_piawka <- annotated_filtered_piawka[ #filter annot piawka by specific window using varia id
    annotated_filtered_piawka$varia_id %in% selected_geneic_window$varia_id,
  ]
  
  message("Filtered to ", nrow(selected_geneic_window_annotated_filtered_piawka), " sites in window")
  
  
  ##############################################
  # LOAD COVERAGE DATA
  ##############################################
  
  coverage_data <- NULL
  
  if(!is.null(coverage_file) && file.exists(coverage_file)) {
    # Get window boundaries from selected_geneic_window
    window_start <- min(selected_geneic_window$gene_start)
    window_end <- max(selected_geneic_window$gene_end)
    
    coverage_data <- load_coverage_data(
      coverage_file = coverage_file,
      focal_gene_start = focal_gene_start,
      focal_gene_end = focal_gene_end,
      window_start = window_start,
      window_end = window_end
    )
  } else {
    message("No coverage file provided or file not found")
  }
  
  ##############################################
  # EXTRACT STATISTICS FOR DIPLOID
  ##############################################
  
  message("Extracting diploid statistics...")
  
  diploid_data <- selected_geneic_window_annotated_filtered_piawka[
    (selected_geneic_window_annotated_filtered_piawka$pop1 == "low_diploid" & 
       selected_geneic_window_annotated_filtered_piawka$pop2 == "high_diploid") |
      (selected_geneic_window_annotated_filtered_piawka$pop1 == "high_diploid" & 
         selected_geneic_window_annotated_filtered_piawka$pop2 == "low_diploid"),
  ]
  
  # # *** COUNT SAMPLES IN DIPLOID COMPARISONS ***
  # n_samples_dip_low <- length(unique(piawka$sample[piawka$pop1 == "low_diploid" | piawka$pop2 == "."]))
  # n_samples_dip_high <- length(unique(piawka$sample[piawka$pop1 == "high_diploid" | piawka$pop2 == "."]))
  
  # message("Sample counts:")
  # message("  Diploid low: ", n_samples_dip_low)
  # message("  Diploid high: ", n_samples_dip_high)


  
  diploid_fst <- diploid_data[diploid_data$metric == "Fst_WC", 
                              c("locus", "value", "pos1", "gene", "varia_id")]
  names(diploid_fst)[names(diploid_fst) == "value"] <- "fst"
  
  diploid_pi_low <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "low_diploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "pi",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(diploid_pi_low)[names(diploid_pi_low) == "value"] <- "pi_low"
  
  diploid_pi_high <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "high_diploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "pi",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(diploid_pi_high)[names(diploid_pi_high) == "value"] <- "pi_high"
  
  diploid_tajd_low <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "low_diploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "TajD",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(diploid_tajd_low)[names(diploid_tajd_low) == "value"] <- "tajd_low"
  
  diploid_tajd_high <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "high_diploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "TajD",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(diploid_tajd_high)[names(diploid_tajd_high) == "value"] <- "tajd_high"
  
  diploid_stats <- diploid_fst
  diploid_stats$fst[diploid_stats$fst < 0] <- 0
  diploid_stats$start_pos <- diploid_stats$pos1
  diploid_stats <- diploid_stats[order(diploid_stats$start_pos), ]
  
  
  ##############################################
  # EXTRACT STATISTICS FOR TETRAPLOID
  ##############################################
  
  message("Extracting tetraploid statistics...")
  
  tetraploid_data <- selected_geneic_window_annotated_filtered_piawka[
    (selected_geneic_window_annotated_filtered_piawka$pop1 == "low_tetraploid" & 
       selected_geneic_window_annotated_filtered_piawka$pop2 == "high_tetraploid") |
      (selected_geneic_window_annotated_filtered_piawka$pop1 == "high_tetraploid" & 
         selected_geneic_window_annotated_filtered_piawka$pop2 == "low_tetraploid"),
  ]
  
  # # *** COUNT SAMPLES IN TETRAPLOID COMPARISONS ***
  # n_samples_tet_low <- length(unique(piawka$sample[piawka$pop1 == "low_tetraploid" | piawka$pop2 == "low_tetraploid"]))
  # n_samples_tet_high <- length(unique(piawka$sample[piawka$pop1 == "high_tetraploid" | piawka$pop2 == "high_tetraploid"]))
  # message("Sample counts:")
  # message("  Tetraploid low: ", n_samples_tet_low)
  # message("  Tetraploid high: ", n_samples_tet_high)
  
  tetraploid_fst <- tetraploid_data[tetraploid_data$metric == "Fst_WC", 
                                    c("locus", "value", "pos1", "gene", "varia_id")]
  names(tetraploid_fst)[names(tetraploid_fst) == "value"] <- "fst"
  
  tetraploid_pi_low <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "low_tetraploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "pi",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(tetraploid_pi_low)[names(tetraploid_pi_low) == "value"] <- "pi_low"
  
  tetraploid_pi_high <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "high_tetraploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "pi",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(tetraploid_pi_high)[names(tetraploid_pi_high) == "value"] <- "pi_high"
  
  tetraploid_tajd_low <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "low_tetraploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "TajD",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(tetraploid_tajd_low)[names(tetraploid_tajd_low) == "value"] <- "tajd_low"
  
  tetraploid_tajd_high <- selected_geneic_window_annotated_filtered_piawka[
    selected_geneic_window_annotated_filtered_piawka$pop1 == "high_tetraploid" & 
      selected_geneic_window_annotated_filtered_piawka$pop2 == "." &
      selected_geneic_window_annotated_filtered_piawka$metric == "TajD",
    c("locus", "value", "pos1", "gene", "varia_id")
  ]
  names(tetraploid_tajd_high)[names(tetraploid_tajd_high) == "value"] <- "tajd_high"
  
  tetraploid_stats <- tetraploid_fst
  tetraploid_stats$fst[tetraploid_stats$fst < 0] <- 0
  tetraploid_stats$start_pos <- tetraploid_stats$pos1
  tetraploid_stats <- tetraploid_stats[order(tetraploid_stats$start_pos), ]
  
#   ##############################################
#   # EXTRACT BACKGROUND FST
#   ##############################################
  
#   message("Extracting background Fst...")
  
#   diploid_background_fst <- selected_geneic_window_annotated_filtered_piawka[
#     selected_geneic_window_annotated_filtered_piawka$metric == "Fst_WC" &
#       !(
#         (selected_geneic_window_annotated_filtered_piawka$pop1 == "low_diploid" & 
#            selected_geneic_window_annotated_filtered_piawka$pop2 == "high_diploid") |
#           (selected_geneic_window_annotated_filtered_piawka$pop1 == "high_diploid" & 
#              selected_geneic_window_annotated_filtered_piawka$pop2 == "low_diploid")
#       ),
#     c("locus", "value", "pos1", "pop1", "pop2")
#   ]
#   names(diploid_background_fst)[names(diploid_background_fst) == "value"] <- "fst"
#   diploid_background_fst$fst[diploid_background_fst$fst < 0] <- 0
  
#   tetraploid_background_fst <- selected_geneic_window_annotated_filtered_piawka[
#     selected_geneic_window_annotated_filtered_piawka$metric == "Fst_WC" &
#       !(
#         (selected_geneic_window_annotated_filtered_piawka$pop1 == "low_tetraploid" & 
#            selected_geneic_window_annotated_filtered_piawka$pop2 == "high_tetraploid") |
#           (selected_geneic_window_annotated_filtered_piawka$pop1 == "high_tetraploid" & 
#              selected_geneic_window_annotated_filtered_piawka$pop2 == "low_tetraploid")
#       ),
#     c("locus", "value", "pos1", "pop1", "pop2")
#   ]
#   names(tetraploid_background_fst)[names(tetraploid_background_fst) == "value"] <- "fst"
#   tetraploid_background_fst$fst[tetraploid_background_fst$fst < 0] <- 0
  
  ##############################################
  # CALCULATE OUTLIERS (95th percentile of low vs high elevation)
  ##############################################
  
  message("Calculating outliers...")
  
  if(nrow(diploid_stats) == 0) {
    message("WARNING: No diploid Fst data for ", degeneracy_label)
    fst_threshold_dip <- NA
    outliers_dip <- logical(0)
    outliers_in_focal_dip <- logical(0)
  } else {
    fst_threshold_dip <- quantile(diploid_stats$fst, 0.95, na.rm=TRUE)
    outliers_dip <- diploid_stats$fst > fst_threshold_dip
    
    outliers_in_focal_dip <- outliers_dip & 
      diploid_stats$start_pos >= focal_gene_start & 
      diploid_stats$start_pos <= focal_gene_end
  }
  
  if(nrow(tetraploid_stats) == 0) {
    message("WARNING: No tetraploid Fst data for ", degeneracy_label)
    fst_threshold_tet <- NA
    outliers_tet <- logical(0)
    outliers_in_focal_tet <- logical(0)
  } else {
    fst_threshold_tet <- quantile(tetraploid_stats$fst, 0.95, na.rm=TRUE)
    outliers_tet <- tetraploid_stats$fst > fst_threshold_tet
    
    outliers_in_focal_tet <- outliers_tet & 
      tetraploid_stats$start_pos >= focal_gene_start & 
      tetraploid_stats$start_pos <= focal_gene_end
  }
  
  message("Diploid outliers: ", sum(outliers_dip, na.rm=TRUE), " (", sum(outliers_in_focal_dip, na.rm=TRUE), " in focal gene)")
  message("Tetraploid outliers: ", sum(outliers_tet, na.rm=TRUE), " (", sum(outliers_in_focal_tet, na.rm=TRUE), " in focal gene)")
  if(!is.na(fst_threshold_dip)) message("Diploid Fst 95th percentile threshold: ", round(fst_threshold_dip, 4))
  if(!is.na(fst_threshold_tet)) message("Tetraploid Fst 95th percentile threshold: ", round(fst_threshold_tet, 4))
  
  ##############################################
  # EXPORT DATA (should handle EMPTY DATAFRAMES)
  ##############################################
  
  message("Exporting statistics to CSV...")
  
  # Fst statistics
  if(nrow(diploid_stats) > 0) {
    diploid_all_stats <- diploid_stats
    diploid_all_stats$ploidy <- "diploid"
    diploid_all_stats$is_outlier <- outliers_dip
    diploid_all_stats$in_focal_gene <- diploid_all_stats$start_pos >= focal_gene_start & 
      diploid_all_stats$start_pos <= focal_gene_end
    # diploid_all_stats$n_samples_low <- n_samples_dip_low
    # diploid_all_stats$n_samples_high <- n_samples_dip_high
    diploid_all_stats$fst_95th_threshold <- fst_threshold_dip
  } else {
    diploid_all_stats <- data.frame()
  }
  
  if(nrow(tetraploid_stats) > 0) {
    tetraploid_all_stats <- tetraploid_stats
    tetraploid_all_stats$ploidy <- "tetraploid"
    tetraploid_all_stats$is_outlier <- outliers_tet
    tetraploid_all_stats$in_focal_gene <- tetraploid_all_stats$start_pos >= focal_gene_start & 
      tetraploid_all_stats$start_pos <= focal_gene_end
    # tetraploid_all_stats$n_samples_low <- n_samples_tet_low
    # tetraploid_all_stats$n_samples_high <- n_samples_tet_high
    tetraploid_all_stats$fst_95th_threshold <- fst_threshold_tet
  } else {
    tetraploid_all_stats <- data.frame()
  }
  
  # Only rbind if we have data
  if(nrow(diploid_all_stats) > 0 && nrow(tetraploid_all_stats) > 0) {
    all_fst_data <- rbind(
      diploid_all_stats[, c("locus", "start_pos", "gene", "varia_id", "fst", "ploidy", "is_outlier", "in_focal_gene", "fst_95th_threshold")],
      tetraploid_all_stats[, c("locus", "start_pos", "gene", "varia_id", "fst", "ploidy", "is_outlier", "in_focal_gene", "fst_95th_threshold")]
    )
  } else if(nrow(diploid_all_stats) > 0) {
    all_fst_data <- diploid_all_stats[, c("locus", "start_pos", "gene", "varia_id", "fst", "ploidy", "is_outlier", "in_focal_gene", "fst_95th_threshold")]
  } else if(nrow(tetraploid_all_stats) > 0) {
    all_fst_data <- tetraploid_all_stats[, c("locus", "start_pos", "gene", "varia_id", "fst", "ploidy", "is_outlier", "in_focal_gene", "fst_95th_threshold")]
  } else {
    message("WARNING: No Fst data to export!")
    all_fst_data <- data.frame()
  }
  
  if(nrow(all_fst_data) > 0) {
    write.csv(all_fst_data, 
              file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_all_fst_statistics_", degeneracy_label, ".csv")), 
              row.names = FALSE)
  }
  
  # Pi statistics (NO DXY) - FIXED FOR EMPTY DATAFRAMES
  all_pi_data_list <- list()
  
  # Helper function to process pi data
  process_pi_data <- function(df, df_name, ploidy_val, elevation_val) {
    if(nrow(df) == 0) {
      message("Warning: No data for ", ploidy_val, " ", elevation_val, " pi")
      return(NULL)
    }
    
    # Make a copy
    df_export <- df
    
    # Rename value column
    value_col <- grep("^pi_", names(df_export), value = TRUE)
    if(length(value_col) > 0) {
      names(df_export)[names(df_export) == value_col] <- "value"
    }
    
    # Add metadata
    df_export$ploidy <- ploidy_val
    df_export$elevation <- elevation_val
    df_export$metric <- "pi"
    df_export$in_focal_gene <- df_export$pos1 >= focal_gene_start & 
      df_export$pos1 <= focal_gene_end
    
    return(df_export[, c("locus", "pos1", "gene", "varia_id", "value", "ploidy", "elevation", "metric", "in_focal_gene")])
  }
  
  # Process each pi dataset
  all_pi_data_list[[1]] <- process_pi_data(diploid_pi_low, "diploid_pi_low", "diploid", "low")
  all_pi_data_list[[2]] <- process_pi_data(diploid_pi_high, "diploid_pi_high", "diploid", "high")
  all_pi_data_list[[3]] <- process_pi_data(tetraploid_pi_low, "tetraploid_pi_low", "tetraploid", "low")
  all_pi_data_list[[4]] <- process_pi_data(tetraploid_pi_high, "tetraploid_pi_high", "tetraploid", "high")
  
  # Remove NULL entries and combine
  all_pi_data_list <- all_pi_data_list[!sapply(all_pi_data_list, is.null)]
  
  if(length(all_pi_data_list) > 0) {
    all_pi_data <- do.call(rbind, all_pi_data_list)
    write.csv(all_pi_data, 
              file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_all_pi_statistics_", degeneracy_label, ".csv")), 
              row.names = FALSE)
  } else {
    message("Warning: No pi data to export")
  }
  

  # Tajima's D statistics - FIXED FOR EMPTY DATAFRAMES
  all_tajd_data_list <- list()
  
  # Helper function to process tajd data
  process_tajd_data <- function(df, df_name, ploidy_val, elevation_val) {
    if(nrow(df) == 0) {
      message("Warning: No data for ", ploidy_val, " ", elevation_val, " TajD")
      return(NULL)
    }
    
    # Make a copy
    df_export <- df
    
    # Rename value column
    value_col <- grep("^tajd_", names(df_export), value = TRUE)
    if(length(value_col) > 0) {
      names(df_export)[names(df_export) == value_col] <- "value"
    }
    
    # Add metadata
    df_export$ploidy <- ploidy_val
    df_export$elevation <- elevation_val
    df_export$metric <- "TajD"
    df_export$in_focal_gene <- df_export$pos1 >= focal_gene_start & 
      df_export$pos1 <= focal_gene_end
    
    return(df_export[, c("locus", "pos1", "gene", "varia_id", "value", "ploidy", "elevation", "metric", "in_focal_gene")])
  }
  
  # Process each tajd dataset
  all_tajd_data_list[[1]] <- process_tajd_data(diploid_tajd_low, "diploid_tajd_low", "diploid", "low")
  all_tajd_data_list[[2]] <- process_tajd_data(diploid_tajd_high, "diploid_tajd_high", "diploid", "high")
  all_tajd_data_list[[3]] <- process_tajd_data(tetraploid_tajd_low, "tetraploid_tajd_low", "tetraploid", "low")
  all_tajd_data_list[[4]] <- process_tajd_data(tetraploid_tajd_high, "tetraploid_tajd_high", "tetraploid", "high")
  
  # Remove NULL entries and combine
  all_tajd_data_list <- all_tajd_data_list[!sapply(all_tajd_data_list, is.null)]
  
  if(length(all_tajd_data_list) > 0) {
    all_tajd_data <- do.call(rbind, all_tajd_data_list)
    write.csv(all_tajd_data, 
              file.path(outdir, paste0(focal_gene_name, "_", focal_id,"_all_tajd_statistics_", degeneracy_label, ".csv")), 
              row.names = FALSE)
  } else {
    message("Warning: No Tajima's D data to export")
  }
  
  # Export outliers
  if(sum(outliers_dip, na.rm=TRUE) > 0) {
    diploid_outliers <- diploid_stats[outliers_dip, c("locus", "start_pos", "gene", "varia_id", "fst")]
    diploid_outliers$in_focal_gene <- diploid_outliers$varia_id == focal_id
    diploid_outliers$ploidy <- "diploid"
    diploid_outliers$fst_95th_threshold <- fst_threshold_dip
    # diploid_outliers$n_samples_low <- n_samples_dip_low
    # diploid_outliers$n_samples_high <- n_samples_dip_high
    write.csv(diploid_outliers, 
              file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_diploid_fst_outliers_", degeneracy_label, ".csv")), 
              row.names = FALSE)
  }
  
  if(sum(outliers_tet, na.rm=TRUE) > 0) {
    tetraploid_outliers <- tetraploid_stats[outliers_tet, c("locus", "start_pos", "gene", "varia_id", "fst")]
    tetraploid_outliers$in_focal_gene <- tetraploid_outliers$varia_id == focal_id
    tetraploid_outliers$ploidy <- "tetraploid"
    tetraploid_outliers$fst_95th_threshold <- fst_threshold_tet
    # tetraploid_outliers$n_samples_low <- n_samples_tet_low
    # tetraploid_outliers$n_samples_high <- n_samples_tet_high
    write.csv(tetraploid_outliers, 
              file.path(outdir, paste0(focal_gene_name, "_", focal_id,"_tetraploid_fst_outliers_", degeneracy_label, ".csv")), 
              row.names = FALSE)
  }
 
  
  # Return all data for plotting
  return(list(
    diploid_stats = diploid_stats,
    tetraploid_stats = tetraploid_stats,
    diploid_pi_low = diploid_pi_low,
    diploid_pi_high = diploid_pi_high,
    tetraploid_pi_low = tetraploid_pi_low,
    tetraploid_pi_high = tetraploid_pi_high,
    diploid_tajd_low = diploid_tajd_low,
    diploid_tajd_high = diploid_tajd_high,
    tetraploid_tajd_low = tetraploid_tajd_low,
    tetraploid_tajd_high = tetraploid_tajd_high,
    fst_threshold_dip = fst_threshold_dip,
    fst_threshold_tet = fst_threshold_tet,
    outliers_dip = outliers_dip,
    outliers_tet = outliers_tet,
    all_fst_data = all_fst_data,
    all_pi_data = if(exists("all_pi_data")) all_pi_data else NULL,
    all_tajd_data = if(exists("all_tajd_data")) all_tajd_data else NULL,
    coverage_data = coverage_data
    # n_samples_dip_low = n_samples_dip_low,
    # n_samples_dip_high = n_samples_dip_high,
    # n_samples_tet_low = n_samples_tet_low,
    # n_samples_tet_high = n_samples_tet_high
  ))
}

################################################################################
# PLOTTING FUNCTIONS
################################################################################

##############################################
# COLOR SCHEMES
##############################################

#fst colors
col_fst_dip <- "dodgerblue"       
col_fst_tet <- "#f039bf"                   
col_outlier <- "red"              

#pi colors
col_pi_dip_low <- "#4682B4"        
col_pi_dip_high <- "#191970" 
col_pi_tet_low <- "#ff3e9b"       
col_pi_tet_high <- "#980ec2"     

#tajd colors
col_tajd_dip_low <- "#4682B4"     
col_tajd_dip_high <- "#191970"    
col_tajd_tet_low <- "#ff3e9b"      
col_tajd_tet_high <- "#980ec2" 

#pi0/pi4 ratio colors
col_focal_ratio="#000000"
col_window_ratio="gray30"

##############################################
# FUNCTION: PLOT SINGLE DEGENERACY
##############################################

plot_degeneracy <- function(results, degeneracy_label) {
  
  message("\nCreating plots for ", degeneracy_label, " sites...")
  
  # Extract results
  diploid_stats <- results$diploid_stats
  tetraploid_stats <- results$tetraploid_stats
  diploid_pi_low <- results$diploid_pi_low
  diploid_pi_high <- results$diploid_pi_high
  tetraploid_pi_low <- results$tetraploid_pi_low
  tetraploid_pi_high <- results$tetraploid_pi_high
  diploid_tajd_low <- results$diploid_tajd_low
  diploid_tajd_high <- results$diploid_tajd_high
  tetraploid_tajd_low <- results$tetraploid_tajd_low
  tetraploid_tajd_high <- results$tetraploid_tajd_high
  fst_threshold_dip <- results$fst_threshold_dip
  fst_threshold_tet <- results$fst_threshold_tet
  outliers_dip <- results$outliers_dip
  outliers_tet <- results$outliers_tet
  coverage_data <- results$coverage_data
  # n_samples_dip_low <- results$n_samples_dip_low
  # n_samples_dip_high <- results$n_samples_dip_high
  # n_samples_tet_low <- results$n_samples_tet_low
  # n_samples_tet_high <- results$n_samples_tet_high
  
  # Check if we have enough data to plot
  if(nrow(diploid_stats) == 0 && nrow(tetraploid_stats) == 0) {
    message("WARNING: No Fst data to plot for ", degeneracy_label)
    return(NULL)
  }

  ##############################################
  # CALCULATE GLOBAL X-AXIS LIMITS (FOR CONSISTENT SHADING)
  ##############################################
  
  message("Calculating global x-axis limits...")
  
  # Collect all genomic positions across all datasets
  all_positions <- c()
  
  if(!is.null(coverage_data)) {
    if(coverage_data$use_binned && !is.null(coverage_data$binned)) {
      all_positions <- c(all_positions, coverage_data$binned$bin)
    } else {
      all_positions <- c(all_positions, coverage_data$raw$pos)
    }
  }
  
  if(nrow(diploid_stats) > 0) all_positions <- c(all_positions, diploid_stats$start_pos)
  if(nrow(tetraploid_stats) > 0) all_positions <- c(all_positions, tetraploid_stats$start_pos)
  if(nrow(diploid_pi_low) > 0) all_positions <- c(all_positions, diploid_pi_low$pos1)
  if(nrow(diploid_pi_high) > 0) all_positions <- c(all_positions, diploid_pi_high$pos1)
  if(nrow(tetraploid_pi_low) > 0) all_positions <- c(all_positions, tetraploid_pi_low$pos1)
  if(nrow(tetraploid_pi_high) > 0) all_positions <- c(all_positions, tetraploid_pi_high$pos1)
  if(nrow(diploid_tajd_low) > 0) all_positions <- c(all_positions, diploid_tajd_low$pos1)
  if(nrow(diploid_tajd_high) > 0) all_positions <- c(all_positions, diploid_tajd_high$pos1)
  if(nrow(tetraploid_tajd_low) > 0) all_positions <- c(all_positions, tetraploid_tajd_low$pos1)
  if(nrow(tetraploid_tajd_high) > 0) all_positions <- c(all_positions, tetraploid_tajd_high$pos1)
  
  # Calculate global x-axis limits
  global_xlim <- c(min(all_positions, na.rm = TRUE), max(all_positions, na.rm = TRUE))
  
  message("Global x-axis range: ", format(global_xlim[1], big.mark=","), " - ", format(global_xlim[2], big.mark=","))
  message("Focal gene range: ", format(focal_gene_start, big.mark=","), " - ", format(focal_gene_end, big.mark=","))

  
  ##############################################
  # IDENTIFY OUTLIER POSITIONS AND LOCI (ONLY IN FOCAL GENE)
  ##############################################
  
  # Get outlier positions for diploid (ONLY in focal gene)
  outlier_positions_dip <- NULL
  outlier_loci_dip <- NULL
  if(length(outliers_dip) > 0 && sum(outliers_dip, na.rm=TRUE) > 0) {
    # Filter for focal gene only
    outliers_in_focal_dip <- outliers_dip & 
      diploid_stats$start_pos >= focal_gene_start & 
      diploid_stats$start_pos <= focal_gene_end
    
    if(sum(outliers_in_focal_dip) > 0) {
      outlier_positions_dip <- diploid_stats$start_pos[outliers_in_focal_dip]
      outlier_loci_dip <- diploid_stats$locus[outliers_in_focal_dip]
    }
  }
  
  # Get outlier positions for tetraploid (ONLY in focal gene)
  outlier_positions_tet <- NULL
  outlier_loci_tet <- NULL
  if(length(outliers_tet) > 0 && sum(outliers_tet, na.rm=TRUE) > 0) {
    # Filter for focal gene only
    outliers_in_focal_tet <- outliers_tet & 
      tetraploid_stats$start_pos >= focal_gene_start & 
      tetraploid_stats$start_pos <= focal_gene_end
    
    if(sum(outliers_in_focal_tet) > 0) {
      outlier_positions_tet <- tetraploid_stats$start_pos[outliers_in_focal_tet]
      outlier_loci_tet <- tetraploid_stats$locus[outliers_in_focal_tet]
    }
  }
  
  # Combine all outlier positions (only focal gene outliers now)
  all_outlier_positions <- c(outlier_positions_dip, outlier_positions_tet)
  all_outlier_loci <- c(outlier_loci_dip, outlier_loci_tet)
  
  # Remove duplicates (if same position is outlier in both)
  if(length(all_outlier_positions) > 0) {
    unique_idx <- !duplicated(all_outlier_positions)
    all_outlier_positions <- all_outlier_positions[unique_idx]
    all_outlier_loci <- all_outlier_loci[unique_idx]
  }
  
  message("Found ", length(all_outlier_positions), " outlier positions in focal gene")
  
  # Calculate y-axis limits safely
  fst_values <- c(diploid_stats$fst, tetraploid_stats$fst)
  fst_values <- fst_values[!is.na(fst_values) & is.finite(fst_values)]
  
  if(length(fst_values) == 0) {
    message("WARNING: No valid Fst values to plot")
    return(NULL)
  }
  
  y_max_fst_combined <- max(fst_values) * 1.1
  
  # Pi values
  pi_values <- c(diploid_pi_low$pi_low, diploid_pi_high$pi_high,
                 tetraploid_pi_low$pi_low, tetraploid_pi_high$pi_high)
  pi_values <- pi_values[!is.na(pi_values) & is.finite(pi_values)]
  
  y_max_pi_combined <- if(length(pi_values) > 0) max(pi_values) * 1.1 else 1
  
  # Tajima's D values
  tajd_values <- c(diploid_tajd_low$tajd_low, diploid_tajd_high$tajd_high,
                   tetraploid_tajd_low$tajd_low, tetraploid_tajd_high$tajd_high)
  tajd_values <- tajd_values[!is.na(tajd_values) & is.finite(tajd_values)]
  
  if(length(tajd_values) > 0) {
    y_min_tajd_combined <- min(tajd_values)
    y_max_tajd_combined <- max(tajd_values)
    
    # Ensure we include -2 and +2 thresholds in the plot
    y_min_tajd_combined <- min(y_min_tajd_combined, -2.2)
    y_max_tajd_combined <- max(y_max_tajd_combined, 2.2)
    
    if(abs(y_max_tajd_combined - y_min_tajd_combined) < 1) {
      y_min_tajd_combined <- -2.5
      y_max_tajd_combined <- 2.5
    } else {
      y_min_tajd_combined <- y_min_tajd_combined * 1.1
      y_max_tajd_combined <- y_max_tajd_combined * 1.1
    }
  } else {
    y_min_tajd_combined <- -2.5
    y_max_tajd_combined <- 2.5
  }
  
  ##############################################
  # HELPER FUNCTION: ADD VERTICAL LINES FOR OUTLIERS
  ##############################################
  
  add_outlier_lines <- function(positions, y_min, y_max, lty = 3, col = "red", alpha = 0.3) {
    if(length(positions) > 0) {
      for(pos in positions) {
        abline(v = pos, lty = lty, col = adjustcolor(col, alpha.f = alpha), lwd = 1)
      }
    }
  }
  
  ##############################################
  # HELPER FUNCTION: SMART TEXT LABEL PLACEMENT (NO OVERLAP)
  ##############################################
  
  smart_text_labels <- function(x_positions, y_positions, labels, pos = 3, cex = 0.7, col = "red", offset = 0.5) {
    if(length(x_positions) == 0) return()
    
    # Sort by x position
    order_idx <- order(x_positions)
    x_positions <- x_positions[order_idx]
    y_positions <- y_positions[order_idx]
    labels <- labels[order_idx]
    
    # If only one label, place it normally
    if(length(x_positions) == 1) {
      text(x_positions, y_positions, labels = labels, pos = pos, cex = cex, col = col, offset = offset)
      return()
    }
    
    # For multiple labels, alternate positions to avoid overlap
    for(i in seq_along(x_positions)) {
      # Alternate between top (pos=3) and top-right (pos=4)
      current_pos <- if(i %% 2 == 1) 3 else 4
      # Slightly adjust offset for alternating positions
      current_offset <- if(i %% 2 == 1) offset else offset * 0.7
      
      text(x_positions[i], y_positions[i], 
           labels = labels[i], 
           pos = current_pos, 
           cex = cex, 
           col = col, 
           offset = current_offset)
    }
  }
  
  # Determine number of panels
  n_panels <- if(!is.null(coverage_data)) 4 else 3
  
  # *** USE PNG WITH CAIRO FOR CLUSTER COMPATIBILITY ***
  png_file <- file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_plot_combined_", degeneracy_label, ".png"))
  
  # Use tryCatch to ensure dev.off() is called even if plotting fails
  tryCatch({
    
    # *** PNG WITH CAIRO (CLUSTER-SAFE) ***
    png(png_file, width = 1800, height = n_panels * 750, res = 150, type = "cairo")
    
    # *** INCREASED FONT SIZES ***
    par(mfrow = c(n_panels, 1), mar = c(6, 6, 4.5, 2))
    
    ##############################################
    # PANEL A: COVERAGE
    ##############################################
    
    if(!is.null(coverage_data)) {
      # Determine which data to plot
      if(coverage_data$use_binned && !is.null(coverage_data$binned)) {
        plot_data <- coverage_data$binned
        x_col <- "bin"
        y_col <- "mean_depth"
      } else {
        plot_data <- coverage_data$raw
        x_col <- "pos"
        y_col <- "depth"
      }
      
      # Calculate y-axis limit
      y_max <- max(plot_data[[y_col]], na.rm=TRUE) * 1.3
      
      # Create plot with GLOBAL x-axis limits
      plot(plot_data[[x_col]], plot_data[[y_col]],
           type = "n",
           xlim = global_xlim,  # *** GLOBAL X-AXIS ***
           ylim = c(0, y_max),
           main = paste0("Sequencing Coverage ", focal_gene_name, " region: ", focal_id),
           ylab = "Coverage Depth (reads)",
           xlab = "Genomic Position",
           xaxt = "n",
           cex.main = 1.6,   # *** INCREASED ***
           cex.lab = 1.5,    # *** INCREASED ***
           cex.axis = 1.3)   # *** INCREASED ***
      
      # Custom x-axis without scientific notation
      axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
           cex.axis = 1.3)
      
      # Shade focal gene region
      rect(focal_gene_start, 0, focal_gene_end, y_max,
           col = rgb(1, 0.8, 0.8, 0.4), border = "darkred", lwd = 1.5)
      
      # Plot coverage as filled area
      if(coverage_data$use_binned && !is.null(coverage_data$binned)) {
        polygon(c(plot_data$bin, rev(plot_data$bin)),
                c(plot_data$min_depth, rev(plot_data$max_depth)),
                col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
        
        lines(plot_data$bin, plot_data$mean_depth, 
              col = "darkblue", lwd = 2)
      } else {
        polygon(c(plot_data[[x_col]], rev(plot_data[[x_col]])),
                c(plot_data[[y_col]], rep(0, nrow(plot_data))),
                col = rgb(0.2, 0.4, 0.8, 0.5), border = NA)
        
        lines(plot_data[[x_col]], plot_data[[y_col]], 
              col = "darkblue", lwd = 1.5)
      }
      
      # Add horizontal line for mean coverage
      abline(h = coverage_data$mean_depth, lty = 2, col = "red", lwd = 1.5)
      
      # Add text label for mean coverage
      text(global_xlim[1] + (global_xlim[2] - global_xlim[1]) * 0.02, 
           coverage_data$mean_depth,
           labels = paste0("Mean: ", round(coverage_data$mean_depth, 1), "x"),
           pos = 3, cex = 1.2, col = "red", offset = 0.3)  # *** INCREASED ***
      
      # Add median line
      abline(h = coverage_data$median_depth, lty = 2, col = "darkgreen", lwd = 1)
      
      # Add text label for median coverage
      text(global_xlim[2] - (global_xlim[2] - global_xlim[1]) * 0.02, 
           coverage_data$median_depth,
           labels = paste0("Median: ", round(coverage_data$median_depth, 1), "x"),
           pos = 3, cex = 1.2, col = "darkgreen", offset = 0.3)  # *** INCREASED ***
      
      # Add legend
      legend("topleft",
             legend = c(paste0("Coverage depth"), 
                        paste0("Mean (", round(coverage_data$mean_depth, 1), "x)"),
                        paste0("Median (", round(coverage_data$median_depth, 1), "x)")),
             col = c("darkblue", "red", "darkgreen"),
             lty = c(1, 2, 2),
             lwd = c(1.5, 1.5, 1),
             cex = 1.2,  # *** INCREASED ***
             bg = "white")
      
      # Add outlier lines to coverage (only if outliers exist in focal gene)
      if(length(all_outlier_positions) > 0) {
        add_outlier_lines(all_outlier_positions, 0, y_max, col = "red", alpha = 0.4)
      }
      
      add_panel_label("A")
    }
    
    
    ##############################################
    # PANEL B: COMBINED FST
    ##############################################
    
    plot(diploid_stats$start_pos, diploid_stats$fst, type = "n",
         xlim = global_xlim,  # *** GLOBAL X-AXIS ***
         ylim = c(0, y_max_fst_combined),
         main = paste0("Fst (Low vs High) ", degeneracy_label, " sites (", focal_gene_name, ": ", focal_id, ")"),
         ylab = "Fst_WC", 
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.4,  
         cex.lab = 1.5,   
         cex.axis = 1.3)  
    
    # Custom x-axis without scientific notation
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    # Add outlier vertical lines (only focal gene outliers)
    add_outlier_lines(all_outlier_positions, 0, y_max_fst_combined, col = "red", alpha = 0.3)
    
    abline(h = 0, col = "gray30", lwd = 0.5)
    # if(!is.na(fst_threshold_dip)) abline(h = fst_threshold_dip, lty = 2, col = col_fst_dip, lwd = 1.5)
    # if(!is.na(fst_threshold_tet)) abline(h = fst_threshold_tet, lty = 2, col = col_fst_tet, lwd = 1.5)

    ###
    #lines and labels of thresholds
    if(!is.na(fst_threshold_dip)) {abline(h = fst_threshold_dip, lty = 2, col = col_fst_dip, lwd = 1.5)
  
      # Add text label for diploid threshold
      text(global_xlim[2] - (global_xlim[2] - global_xlim[1]) * 0.02, 
          fst_threshold_dip,
          labels = paste0(" Dip 95th: ", round(fst_threshold_dip, 3)),
          pos = 3, cex = 1.0, col = col_fst_dip, offset = 0.3)
    }
    if(!is.na(fst_threshold_tet)) {abline(h = fst_threshold_tet, lty = 2, col = col_fst_tet, lwd = 1.5)
  
      # Add text label for diploid threshold
      text(global_xlim[1] + (global_xlim[2] - global_xlim[1]) * 0.02, 
          fst_threshold_tet,
          labels = paste0(" Tet 95th: ", round(fst_threshold_tet, 3)),
          pos = 3, cex = 1.0, col = col_fst_tet, offset = 0.3)
    }
    ###
    
    rect(focal_gene_start, 0, focal_gene_end, y_max_fst_combined,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_stats) > 0 && length(outliers_dip) > 0) {
      points(diploid_stats$start_pos[!outliers_dip], diploid_stats$fst[!outliers_dip], 
             pch = 19, cex = 1.4, col = col_fst_dip)  # *** INCREASED ***
    }
    if(nrow(tetraploid_stats) > 0 && length(outliers_tet) > 0) {
      points(tetraploid_stats$start_pos[!outliers_tet], tetraploid_stats$fst[!outliers_tet], 
             pch = 17, cex = 1.4, col = col_fst_tet)  # *** INCREASED ***
    }
    
    # Plot outliers with smart label placement
    if(length(outliers_dip) > 0 && sum(outliers_dip, na.rm=TRUE) > 0) {
      points(diploid_stats$start_pos[outliers_dip], diploid_stats$fst[outliers_dip], 
             pch = 21, cex = 1.6, col = col_outlier, bg = col_outlier, lwd = 2)  # *** INCREASED ***
      
      smart_text_labels(
        x_positions = diploid_stats$start_pos[outliers_dip],
        y_positions = diploid_stats$fst[outliers_dip],
        labels = diploid_stats$locus[outliers_dip],
        pos = 3,
        cex = 0.7,  # *** INCREASED ***
        col = col_outlier,
        offset = 0.5
      )
    }
    
    if(length(outliers_tet) > 0 && sum(outliers_tet, na.rm=TRUE) > 0) {
      points(tetraploid_stats$start_pos[outliers_tet], tetraploid_stats$fst[outliers_tet], 
             pch = 24, cex = 1.6, col = col_outlier, bg = col_outlier, lwd = 2)  # *** INCREASED ***
      
      smart_text_labels(
        x_positions = tetraploid_stats$start_pos[outliers_tet],
        y_positions = tetraploid_stats$fst[outliers_tet],
        labels = tetraploid_stats$locus[outliers_tet],
        pos = 3,
        cex = 0.7,  # *** INCREASED ***
        col = col_outlier,
        offset = 0.5
      )
    }
    
    # Build legend text with threshold values
    legend_text <- c("Diploid", "Tetraploid", "Background", "Outlier (>95th)")
    # if(!is.na(fst_threshold_dip)) {
    #   legend_text[4] <- paste0("Outlier (>", round(fst_threshold_dip, 3), ")")
    # }
    
    legend("topleft",
           legend = legend_text,
           col = c(col_fst_dip, col_fst_tet, col_outlier),
           pch = c(19, 17, 21),
           pt.cex = c(1.4, 1.4, 1.5),  # *** INCREASED ***
           pt.bg = c(NA, NA, col_outlier),
           cex = 1.2,  # *** INCREASED ***
           bg = "white")
    
    add_panel_label("B")


    ##############################################
    # PANEL C: COMBINED NUCLEOTIDE DIVERSITY (NO LABELS)
    ##############################################
    
    if(length(pi_values) > 0) {
      plot(diploid_pi_low$pos1, diploid_pi_low$pi_low, type = "n",
           xlim = global_xlim,  # *** GLOBAL X-AXIS ***
           ylim = c(0, y_max_pi_combined),
           main = paste0("Nucleotide Diversity ", degeneracy_label, " sites (", focal_gene_name, ": ", focal_id, ")"),
           ylab = "(Nucleotide Diversity)", 
           xlab = "Genomic Position",
           xaxt = "n",
           cex.main = 1.4,   
           cex.lab = 1.5,    
           cex.axis = 1.3)  
      
      # Custom x-axis without scientific notation 
      axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
           cex.axis = 1.3)
      
      # Add outlier vertical lines (only focal gene outliers)
      add_outlier_lines(all_outlier_positions, 0, y_max_pi_combined, col = "red", alpha = 0.3)
      
      abline(h = 0, col = "gray30", lwd = 0.5)
      
      rect(focal_gene_start, 0, focal_gene_end, y_max_pi_combined,
           col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
      
      if(nrow(diploid_pi_low) > 0) {
        points(diploid_pi_low$pos1, diploid_pi_low$pi_low, pch=19, cex=1.4, col=col_pi_dip_low)
      }
      if(nrow(diploid_pi_high) > 0) {
        points(diploid_pi_high$pos1, diploid_pi_high$pi_high, pch=19, cex=1.4, col=col_pi_dip_high)
      }
      if(nrow(tetraploid_pi_low) > 0) {
        points(tetraploid_pi_low$pos1, tetraploid_pi_low$pi_low, pch=17, cex=1.4, col=col_pi_tet_low)
      }
      if(nrow(tetraploid_pi_high) > 0) {
        points(tetraploid_pi_high$pos1, tetraploid_pi_high$pi_high, pch=17, cex=1.4, col=col_pi_tet_high)
      }
      
      # Highlight outlier positions with circles ONLY (NO LABELS)
      if(length(all_outlier_positions) > 0) {
        
        if(nrow(diploid_pi_low) > 0) {
          matches <- diploid_pi_low$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(diploid_pi_low$pos1[matches], 
                   diploid_pi_low$pi_low[matches],
                   pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
        
        if(nrow(diploid_pi_high) > 0) {
          matches <- diploid_pi_high$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(diploid_pi_high$pos1[matches], 
                   diploid_pi_high$pi_high[matches],
                   pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
        
        if(nrow(tetraploid_pi_low) > 0) {
          matches <- tetraploid_pi_low$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(tetraploid_pi_low$pos1[matches], 
                   tetraploid_pi_low$pi_low[matches],
                   pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
        
        if(nrow(tetraploid_pi_high) > 0) {
          matches <- tetraploid_pi_high$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(tetraploid_pi_high$pos1[matches], 
                   tetraploid_pi_high$pi_high[matches],
                   pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
      }
      
      legend("topleft",
             legend = c("Dip Low", "Dip High", "Tet Low", "Tet High", "Fst outlier"), 
             col = c(col_pi_dip_low, col_pi_dip_high, col_pi_tet_low, col_pi_tet_high, col_outlier),
             pch = c(19, 19, 17, 17, 1),
             pt.cex = c(1.4, 1.4, 1.4, 1.4, 1.5),  # *** INCREASED ***
             lwd = c(NA, NA, NA, NA, 2.5),
             cex = 1.2,  # *** INCREASED ***
             bg = "white")
    } else {
      plot.new()
      text(0.5, 0.5, "No nucleotide diversity data available", cex = 1.5, col = "gray50")
    }
    
    add_panel_label("C")
    
    ##############################################
    # PANEL D: COMBINED TAJIMA'S D (NO LABELS)
    ##############################################
    
    if(length(tajd_values) > 0) {
      plot(diploid_tajd_low$pos1, diploid_tajd_low$tajd_low, type = "n",
           xlim = global_xlim,  # *** GLOBAL X-AXIS ***
           ylim = c(y_min_tajd_combined, y_max_tajd_combined),
           main = paste0("Tajima's D ", degeneracy_label, " sites (", focal_gene_name, ": ", focal_id, ")"),
           ylab = "Tajima's D", 
           xlab = "Genomic Position",
           xaxt = "n",
           cex.main = 1.4,   
           cex.lab = 1.5,   
           cex.axis = 1.3)  
      
      # Custom x-axis without scientific notation
      axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
           cex.axis = 1.3)
      
      # Add outlier vertical lines (only focal gene outliers)
      add_outlier_lines(all_outlier_positions, y_min_tajd_combined, y_max_tajd_combined, 
                        col = "red", alpha = 0.3)
      
      abline(h = 0, col = "black", lwd = 1.5)
      abline(h = -2, lty = 2, col = "blue", lwd = 1.5)  # *** -2 THRESHOLD ***
      abline(h = 2, lty = 2, col = "blue", lwd = 1.5)   # *** +2 THRESHOLD ***
      
      rect(focal_gene_start, y_min_tajd_combined, focal_gene_end, y_max_tajd_combined,
           col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
      
      if(nrow(diploid_tajd_low) > 0) {
        points(diploid_tajd_low$pos1, diploid_tajd_low$tajd_low, pch=19, cex=1.4, col=col_tajd_dip_low)
      }
      if(nrow(diploid_tajd_high) > 0) {
        points(diploid_tajd_high$pos1, diploid_tajd_high$tajd_high, pch=19, cex=1.4, col=col_tajd_dip_high)
      }
      if(nrow(tetraploid_tajd_low) > 0) {
        points(tetraploid_tajd_low$pos1, tetraploid_tajd_low$tajd_low, pch=17, cex=1.4, col=col_tajd_tet_low)
      }
      if(nrow(tetraploid_tajd_high) > 0) {
        points(tetraploid_tajd_high$pos1, tetraploid_tajd_high$tajd_high, pch=17, cex=1.4, col=col_tajd_tet_high)
      }
      
      # Highlight outlier positions with circles ONLY (NO LABELS)
      if(length(all_outlier_positions) > 0) {
        
        if(nrow(diploid_tajd_low) > 0) {
          matches <- diploid_tajd_low$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(diploid_tajd_low$pos1[matches], 
                   diploid_tajd_low$tajd_low[matches],
                   pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
        
        if(nrow(diploid_tajd_high) > 0) {
          matches <- diploid_tajd_high$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(diploid_tajd_high$pos1[matches], 
                   diploid_tajd_high$tajd_high[matches],
                   pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
        
        if(nrow(tetraploid_tajd_low) > 0) {
          matches <- tetraploid_tajd_low$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(tetraploid_tajd_low$pos1[matches], 
                   tetraploid_tajd_low$tajd_low[matches],
                   pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
        
        if(nrow(tetraploid_tajd_high) > 0) {
          matches <- tetraploid_tajd_high$pos1 %in% all_outlier_positions
          if(sum(matches) > 0) {
            points(tetraploid_tajd_high$pos1[matches], 
                   tetraploid_tajd_high$tajd_high[matches],
                   pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)  # *** INCREASED ***
          }
        }
      }
      
      legend("topleft",
             legend = c("Neutral", "±2 threshold", "Fst outlier"),
             col = c("black", "blue", col_outlier), 
             pch = c(NA, NA, 1),
             lty = c(1, 2, NA),
             pt.cex = c(NA, NA, 2.0),  # *** INCREASED ***
             lwd = c(1.5, 1.5, 2.5),
             cex = 1.2,  # *** INCREASED ***
             bg = "white")
    } else {
      plot.new()
      text(0.5, 0.5, "No Tajima's D data available", cex = 1.5, col = "gray50")
    }
    
    add_panel_label("D")
    
    # Close device
    dev.off()
    
    message("Plot saved: ", png_file)
    
  }, error = function(e) {
    # If there's an error, make sure to close the device
    if(dev.cur() > 1) dev.off()
    message("ERROR creating plot: ", e$message)
    return(NULL)
  })
}


#plot both side by side
################################################################################
# FUNCTION: PLOT SIDE-BY-SIDE COMPARISON (0-FOLD VS 4-FOLD)
# Exact replica of plot_degeneracy for both degeneracies side-by-side
################################################################################

plot_comparison_side_by_side <- function(results_0fold, results_4fold) {
  
  message("\nCreating side-by-side comparison plot (0-fold vs 4-fold)...")
  
  ##############################################
  # EXTRACT RESULTS FROM BOTH DATASETS
  ##############################################
  
  # 0-fold data
  diploid_stats_0 <- results_0fold$diploid_stats
  tetraploid_stats_0 <- results_0fold$tetraploid_stats
  diploid_pi_low_0 <- results_0fold$diploid_pi_low
  diploid_pi_high_0 <- results_0fold$diploid_pi_high
  tetraploid_pi_low_0 <- results_0fold$tetraploid_pi_low
  tetraploid_pi_high_0 <- results_0fold$tetraploid_pi_high
  diploid_tajd_low_0 <- results_0fold$diploid_tajd_low
  diploid_tajd_high_0 <- results_0fold$diploid_tajd_high
  tetraploid_tajd_low_0 <- results_0fold$tetraploid_tajd_low
  tetraploid_tajd_high_0 <- results_0fold$tetraploid_tajd_high
  fst_threshold_dip_0 <- results_0fold$fst_threshold_dip
  fst_threshold_tet_0 <- results_0fold$fst_threshold_tet
  outliers_dip_0 <- results_0fold$outliers_dip
  outliers_tet_0 <- results_0fold$outliers_tet
  coverage_data_0 <- results_0fold$coverage_data
  
  # 4-fold data
  diploid_stats_4 <- results_4fold$diploid_stats
  tetraploid_stats_4 <- results_4fold$tetraploid_stats
  diploid_pi_low_4 <- results_4fold$diploid_pi_low
  diploid_pi_high_4 <- results_4fold$diploid_pi_high
  tetraploid_pi_low_4 <- results_4fold$tetraploid_pi_low
  tetraploid_pi_high_4 <- results_4fold$tetraploid_pi_high
  diploid_tajd_low_4 <- results_4fold$diploid_tajd_low
  diploid_tajd_high_4 <- results_4fold$diploid_tajd_high
  tetraploid_tajd_low_4 <- results_4fold$tetraploid_tajd_low
  tetraploid_tajd_high_4 <- results_4fold$tetraploid_tajd_high
  fst_threshold_dip_4 <- results_4fold$fst_threshold_dip
  fst_threshold_tet_4 <- results_4fold$fst_threshold_tet
  outliers_dip_4 <- results_4fold$outliers_dip
  outliers_tet_4 <- results_4fold$outliers_tet
  coverage_data_4 <- results_4fold$coverage_data
  
  # # Sample counts
  # n_samples_dip_low <- results_0fold$n_samples_dip_low
  # n_samples_dip_high <- results_0fold$n_samples_dip_high
  # n_samples_tet_low <- results_0fold$n_samples_tet_low
  # n_samples_tet_high <- results_0fold$n_samples_tet_high
  
  ##############################################
  # CALCULATE GLOBAL X-AXIS LIMITS
  ##############################################
  
  all_positions <- c()
  
  # Collect from 0-fold
  if(!is.null(coverage_data_0)) {
    if(coverage_data_0$use_binned && !is.null(coverage_data_0$binned)) {
      all_positions <- c(all_positions, coverage_data_0$binned$bin)
    } else {
      all_positions <- c(all_positions, coverage_data_0$raw$pos)
    }
  }
  if(nrow(diploid_stats_0) > 0) all_positions <- c(all_positions, diploid_stats_0$start_pos)
  if(nrow(tetraploid_stats_0) > 0) all_positions <- c(all_positions, tetraploid_stats_0$start_pos)
  all_positions <- c(all_positions, diploid_pi_low_0$pos1, diploid_pi_high_0$pos1, 
                     tetraploid_pi_low_0$pos1, tetraploid_pi_high_0$pos1,
                     diploid_tajd_low_0$pos1, diploid_tajd_high_0$pos1,
                     tetraploid_tajd_low_0$pos1, tetraploid_tajd_high_0$pos1)
  
  # Collect from 4-fold
  if(!is.null(coverage_data_4)) {
    if(coverage_data_4$use_binned && !is.null(coverage_data_4$binned)) {
      all_positions <- c(all_positions, coverage_data_4$binned$bin)
    } else {
      all_positions <- c(all_positions, coverage_data_4$raw$pos)
    }
  }
  if(nrow(diploid_stats_4) > 0) all_positions <- c(all_positions, diploid_stats_4$start_pos)
  if(nrow(tetraploid_stats_4) > 0) all_positions <- c(all_positions, tetraploid_stats_4$start_pos)
  all_positions <- c(all_positions, diploid_pi_low_4$pos1, diploid_pi_high_4$pos1, 
                     tetraploid_pi_low_4$pos1, tetraploid_pi_high_4$pos1,
                     diploid_tajd_low_4$pos1, diploid_tajd_high_4$pos1,
                     tetraploid_tajd_low_4$pos1, tetraploid_tajd_high_4$pos1)
  
  global_xlim <- c(min(all_positions, na.rm = TRUE), max(all_positions, na.rm = TRUE))
  
  ##############################################
  # GET FOCAL GENE OUTLIER POSITIONS
  ##############################################
  
  get_focal_outliers <- function(stats, outliers) {
    if(nrow(stats) == 0 || length(outliers) == 0) return(c())
    outliers_in_focal <- outliers & 
      stats$start_pos >= focal_gene_start & 
      stats$start_pos <= focal_gene_end
    if(sum(outliers_in_focal) > 0) {
      return(stats$start_pos[outliers_in_focal])
    }
    return(c())
  }
  
  # Get outlier positions for 0-fold
  outlier_positions_0 <- unique(c(
    get_focal_outliers(diploid_stats_0, outliers_dip_0),
    get_focal_outliers(tetraploid_stats_0, outliers_tet_0)
  ))
  
  # Get outlier positions for 4-fold
  outlier_positions_4 <- unique(c(
    get_focal_outliers(diploid_stats_4, outliers_dip_4),
    get_focal_outliers(tetraploid_stats_4, outliers_tet_4)
  ))
  
  message("0-fold outliers in focal gene: ", length(outlier_positions_0))
  message("4-fold outliers in focal gene: ", length(outlier_positions_4))
  
  ##############################################
  # CALCULATE Y-AXIS LIMITS FOR EACH METRIC
  ##############################################
  
  # Coverage y-limits
  if(!is.null(coverage_data_0) && !is.null(coverage_data_4)) {
    if(coverage_data_0$use_binned) {
      max_cov_0 <- max(coverage_data_0$binned$mean_depth, na.rm=TRUE)
    } else {
      max_cov_0 <- max(coverage_data_0$raw$depth, na.rm=TRUE)
    }
    if(coverage_data_4$use_binned) {
      max_cov_4 <- max(coverage_data_4$binned$mean_depth, na.rm=TRUE)
    } else {
      max_cov_4 <- max(coverage_data_4$raw$depth, na.rm=TRUE)
    }
    y_max_cov <- max(max_cov_0, max_cov_4, na.rm=TRUE) * 1.3
  } else if(!is.null(coverage_data_0)) {
    if(coverage_data_0$use_binned) {
      y_max_cov <- max(coverage_data_0$binned$mean_depth, na.rm=TRUE) * 1.3
    } else {
      y_max_cov <- max(coverage_data_0$raw$depth, na.rm=TRUE) * 1.3
    }
  } else if(!is.null(coverage_data_4)) {
    if(coverage_data_4$use_binned) {
      y_max_cov <- max(coverage_data_4$binned$mean_depth, na.rm=TRUE) * 1.3
    } else {
      y_max_cov <- max(coverage_data_4$raw$depth, na.rm=TRUE) * 1.3
    }
  } else {
    y_max_cov <- 100
  }
  
  # Fst y-limits
  fst_values <- c(diploid_stats_0$fst, tetraploid_stats_0$fst,
                  diploid_stats_4$fst, tetraploid_stats_4$fst)
  fst_values <- fst_values[!is.na(fst_values) & is.finite(fst_values)]
  y_max_fst <- if(length(fst_values) > 0) max(fst_values) * 1.1 else 1
  
  # Pi y-limits
  pi_values <- c(diploid_pi_low_0$pi_low, diploid_pi_high_0$pi_high,
                 tetraploid_pi_low_0$pi_low, tetraploid_pi_high_0$pi_high,
                 diploid_pi_low_4$pi_low, diploid_pi_high_4$pi_high,
                 tetraploid_pi_low_4$pi_low, tetraploid_pi_high_4$pi_high)
  pi_values <- pi_values[!is.na(pi_values) & is.finite(pi_values)]
  y_max_pi <- if(length(pi_values) > 0) max(pi_values) * 1.1 else 1
  
  # Tajima's D y-limits
  tajd_values <- c(diploid_tajd_low_0$tajd_low, diploid_tajd_high_0$tajd_high,
                   tetraploid_tajd_low_0$tajd_low, tetraploid_tajd_high_0$tajd_high,
                   diploid_tajd_low_4$tajd_low, diploid_tajd_high_4$tajd_high,
                   tetraploid_tajd_low_4$tajd_low, tetraploid_tajd_high_4$tajd_high)
  tajd_values <- tajd_values[!is.na(tajd_values) & is.finite(tajd_values)]
  
  if(length(tajd_values) > 0) {
    y_min_tajd <- min(tajd_values, -2.2)
    y_max_tajd <- max(tajd_values, 2.2)
    y_min_tajd <- y_min_tajd * 1.1
    y_max_tajd <- y_max_tajd * 1.1
  } else {
    y_min_tajd <- -2.5
    y_max_tajd <- 2.5
  }
  
  ##############################################
  # HELPER FUNCTIONS
  ##############################################
  
  add_outlier_lines <- function(positions, y_min, y_max) {
    if(length(positions) > 0) {
      for(pos in positions) {
        abline(v = pos, lty = 3, col = adjustcolor("red", alpha.f = 0.3), lwd = 1)
      }
    }
  }
  
  smart_text_labels <- function(x_positions, y_positions, labels, pos = 3, cex = 0.7, col = "red", offset = 0.5) {
    if(length(x_positions) == 0) return()
    
    order_idx <- order(x_positions)
    x_positions <- x_positions[order_idx]
    y_positions <- y_positions[order_idx]
    labels <- labels[order_idx]
    
    if(length(x_positions) == 1) {
      text(x_positions, y_positions, labels = labels, pos = pos, cex = cex, col = col, offset = offset)
      return()
    }
    
    for(i in seq_along(x_positions)) {
      current_pos <- if(i %% 2 == 1) 3 else 4
      current_offset <- if(i %% 2 == 1) offset else offset * 0.7
      
      text(x_positions[i], y_positions[i], 
           labels = labels[i], 
           pos = current_pos, 
           cex = cex, 
           col = col, 
           offset = current_offset)
    }
  }
  
  ##############################################
  # CREATE PLOT
  ##############################################
  
  png_file <- file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_comparison_0fold_vs_4fold.png"))
  
  tryCatch({
    
    png(png_file, width = 2400, height = 2400, res = 150, type = "cairo")
    
    # 6 rows x 2 columns layout
    par(mfrow = c(6, 2), mar = c(4.5, 5.5, 3.5, 1.5), oma = c(0, 0, 3, 0))
    #suggested layout changes if my alterations dont work:
    # layout(matrix(c(1,2,
    #             3,4,
    #             5,6,
    #             7,8,
    #             9,10,
    #             11,12),  # NEW ROW
    #           nrow=6, ncol=2, byrow=TRUE),
    #    heights=c(1,1,1,1,1.2,1.2))  # Added height for row 6
    
    ##############################################
    # ROW 1: COVERAGE
    ##############################################
    
    # === PANEL A: 0-FOLD COVERAGE ===
    if(!is.null(coverage_data_0)) {
      if(coverage_data_0$use_binned && !is.null(coverage_data_0$binned)) {
        plot_data <- coverage_data_0$binned
        x_col <- "bin"
        y_col <- "mean_depth"
      } else {
        plot_data <- coverage_data_0$raw
        x_col <- "pos"
        y_col <- "depth"
      }
      
      plot(plot_data[[x_col]], plot_data[[y_col]],
           type = "n",
           xlim = global_xlim,
           ylim = c(0, y_max_cov),
           main = "Coverage Depth Across Gene Window",
           ylab = "Coverage Depth (reads)",
           xlab = "Genomic Position",
           xaxt = "n",
           cex.main = 1.6,
           cex.lab = 1.5,
           cex.axis = 1.3)
      
      axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
           cex.axis = 1.3)
      
      rect(focal_gene_start, 0, focal_gene_end, y_max_cov,
           col = rgb(1, 0.8, 0.8, 0.4), border = "darkred", lwd = 1.5)
      
      if(coverage_data_0$use_binned && !is.null(coverage_data_0$binned)) {
        polygon(c(plot_data$bin, rev(plot_data$bin)),
                c(plot_data$min_depth, rev(plot_data$max_depth)),
                col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
        lines(plot_data$bin, plot_data$mean_depth, col = "darkblue", lwd = 2)
      } else {
        polygon(c(plot_data[[x_col]], rev(plot_data[[x_col]])),
                c(plot_data[[y_col]], rep(0, nrow(plot_data))),
                col = rgb(0.2, 0.4, 0.8, 0.5), border = NA)
        lines(plot_data[[x_col]], plot_data[[y_col]], col = "darkblue", lwd = 1.5)
      }

      abline(h = coverage_data_0$mean_depth, lty = 2, col = "red", lwd = 1.5)
      #text label mean
      text(global_xlim[1] + (global_xlim[2] - global_xlim[1]) * 0.02, 
           coverage_data_0$mean_depth,
           labels = paste0("Mean: ", round(coverage_data_0$mean_depth, 1), "x"),
           pos = 3, cex = 1.2, col = "red", offset = 0.3)
      abline(h = coverage_data_0$median_depth, lty = 2, col = "darkgreen", lwd = 1)
      #text label median
      text(global_xlim[2] - (global_xlim[2] - global_xlim[1]) * 0.02, 
           coverage_data_0$median_depth,
           labels = paste0("Median: ", round(coverage_data_0$median_depth, 1), "x"),
           pos = 3, cex = 1.2, col = "darkgreen", offset = 0.3)      
      
      add_outlier_lines(outlier_positions_0, 0, y_max_cov)

      # Add legend
      legend("topleft",
             legend = c(paste0("Coverage depth"), 
                        paste0("Mean (", round(coverage_data_0$mean_depth, 1), "x)"),
                        paste0("Median (", round(coverage_data_0$median_depth, 1), "x)")),
             col = c("darkblue", "red", "darkgreen"),
             lty = c(1, 2, 2),
             lwd = c(1.5, 1.5, 1),
             cex = 1.2,  # *** INCREASED ***
             bg = "white")      
      
      # *** PANEL LABEL A ***
      add_panel_label("A")
      
    } else {
      plot.new()
      text(0.5, 0.5, "No coverage data", cex = 1.3, col = "gray50")

      add_panel_label("A")
    }
    
    # === PANEL B: 4-FOLD COVERAGE ===
    if(!is.null(coverage_data_4)) {
      if(coverage_data_4$use_binned && !is.null(coverage_data_4$binned)) {
        plot_data <- coverage_data_4$binned
        x_col <- "bin"
        y_col <- "mean_depth"
      } else {
        plot_data <- coverage_data_4$raw
        x_col <- "pos"
        y_col <- "depth"
      }
      
      plot(plot_data[[x_col]], plot_data[[y_col]],
           type = "n",
           xlim = global_xlim,
           ylim = c(0, y_max_cov),
           main = "Coverage Depth Across Gene Window",
           ylab = "Coverage Depth (reads)",
           xlab = "Genomic Position",
           xaxt = "n",
           cex.main = 1.6,
           cex.lab = 1.5,
           cex.axis = 1.3)
      
      axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
           cex.axis = 1.3)
      
      rect(focal_gene_start, 0, focal_gene_end, y_max_cov,
           col = rgb(1, 0.8, 0.8, 0.4), border = "darkred", lwd = 1.5)
      
      if(coverage_data_4$use_binned && !is.null(coverage_data_4$binned)) {
        polygon(c(plot_data$bin, rev(plot_data$bin)),
                c(plot_data$min_depth, rev(plot_data$max_depth)),
                col = rgb(0.2, 0.4, 0.8, 0.3), border = NA)
        lines(plot_data$bin, plot_data$mean_depth, col = "darkblue", lwd = 2)
      } else {
        polygon(c(plot_data[[x_col]], rev(plot_data[[x_col]])),
                c(plot_data[[y_col]], rep(0, nrow(plot_data))),
                col = rgb(0.2, 0.4, 0.8, 0.5), border = NA)
        lines(plot_data[[x_col]], plot_data[[y_col]], col = "darkblue", lwd = 1.5)
      }
      
      abline(h = coverage_data_4$mean_depth, lty = 2, col = "red", lwd = 1.5)
      #text label mean
      text(global_xlim[1] + (global_xlim[2] - global_xlim[1]) * 0.02, 
           coverage_data_4$mean_depth,
           labels = paste0("Mean: ", round(coverage_data_4$mean_depth, 1), "x"),
           pos = 3, cex = 1.2, col = "red", offset = 0.3)

      abline(h = coverage_data_4$median_depth, lty = 2, col = "darkgreen", lwd = 1)
      #text label median
      text(global_xlim[2] - (global_xlim[2] - global_xlim[1]) * 0.02, 
           coverage_data_4$median_depth,
           labels = paste0("Median: ", round(coverage_data_4$median_depth, 1), "x"),
           pos = 3, cex = 1.2, col = "darkgreen", offset = 0.3) 

      add_outlier_lines(outlier_positions_4, 0, y_max_cov)
      
      # Add legend
      legend("topleft",
             legend = c(paste0("Coverage depth"), 
                        paste0("Mean (", round(coverage_data_4$mean_depth, 1), "x)"),
                        paste0("Median (", round(coverage_data_4$median_depth, 1), "x)")),
             col = c("darkblue", "red", "darkgreen"),
             lty = c(1, 2, 2),
             lwd = c(1.5, 1.5, 1),
             cex = 1.2,  # *** INCREASED ***
             bg = "white")

      # *** PANEL LABEL B ***
      add_panel_label("B")
      
    } else {
      plot.new()
      text(0.5, 0.5, "No coverage data", cex = 1.3, col = "gray50")
      add_panel_label("B")
    }
    
    ##############################################
    # ROW 2: FST
    ##############################################     
    
    
    # === PANEL C: 0-FOLD FST ===
    plot(diploid_stats_0$start_pos, diploid_stats_0$fst, type = "n",
         xlim = global_xlim,
         ylim = c(0, y_max_fst),
         main = paste0("Fst (Low vs High) 0-fold sites"),
         ylab = "Fst_WC",
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.4,
         cex.lab = 1.5,
         cex.axis = 1.3)
    
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    add_outlier_lines(outlier_positions_0, 0, y_max_fst)
    
    abline(h = 0, col = "gray30", lwd = 0.5)

    #lines and labels of thresholds
    if(!is.na(fst_threshold_dip_0)) {abline(h = fst_threshold_dip_0, lty = 2, col = col_fst_dip, lwd = 1.5)
  
      # Add text label for diploid threshold
      text(global_xlim[2] - (global_xlim[2] - global_xlim[1]) * 0.02, 
          fst_threshold_dip_0,
          labels = paste0(" Dip 95th: ", round(fst_threshold_dip_0, 3)),
          pos = 3, cex = 1.0, col = col_fst_dip, offset = 0.3)
    }
    if(!is.na(fst_threshold_tet_0)) {abline(h = fst_threshold_tet_0, lty = 2, col = col_fst_tet, lwd = 1.5)
  
      # Add text label for diploid threshold
      text(global_xlim[1] + (global_xlim[2] - global_xlim[1]) * 0.02, 
          fst_threshold_tet_0,
          labels = paste0(" Tet 95th: ", round(fst_threshold_tet_0, 3)),
          pos = 3, cex = 1.0, col = col_fst_tet, offset = 0.3)
    }   

    #just line no label
    # if(!is.na(fst_threshold_dip_0)) abline(h = fst_threshold_dip_0, lty = 2, col = col_fst_dip, lwd = 1.5)
    # if(!is.na(fst_threshold_tet_0)) abline(h = fst_threshold_tet_0, lty = 2, col = col_fst_tet, lwd = 1.5)
    
    rect(focal_gene_start, 0, focal_gene_end, y_max_fst,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_stats_0) > 0 && length(outliers_dip_0) > 0) {
      points(diploid_stats_0$start_pos[!outliers_dip_0], diploid_stats_0$fst[!outliers_dip_0], 
             pch = 19, cex = 1.4, col = col_fst_dip)
      
      if(sum(outliers_dip_0) > 0) {
        focal_outliers_idx <- outliers_dip_0 & 
          diploid_stats_0$start_pos >= focal_gene_start & 
          diploid_stats_0$start_pos <= focal_gene_end
        
        if(sum(focal_outliers_idx) > 0) {
          points(diploid_stats_0$start_pos[focal_outliers_idx], 
                 diploid_stats_0$fst[focal_outliers_idx], 
                 pch = 21, cex = 1.6, col = col_outlier, bg = col_outlier, lwd = 2)
          
          smart_text_labels(
            x_positions = diploid_stats_0$start_pos[focal_outliers_idx],
            y_positions = diploid_stats_0$fst[focal_outliers_idx],
            labels = diploid_stats_0$locus[focal_outliers_idx],
            pos = 3, cex = 0.8, col = col_outlier, offset = 0.5
          )
        }
      }
    }
    
    if(nrow(tetraploid_stats_0) > 0 && length(outliers_tet_0) > 0) {
      points(tetraploid_stats_0$start_pos[!outliers_tet_0], tetraploid_stats_0$fst[!outliers_tet_0], 
             pch = 17, cex = 1.4, col = col_fst_tet)
      
      if(sum(outliers_tet_0) > 0) {
        focal_outliers_idx <- outliers_tet_0 & 
          tetraploid_stats_0$start_pos >= focal_gene_start & 
          tetraploid_stats_0$start_pos <= focal_gene_end
        
        if(sum(focal_outliers_idx) > 0) {
          points(tetraploid_stats_0$start_pos[focal_outliers_idx], 
                 tetraploid_stats_0$fst[focal_outliers_idx], 
                 pch = 24, cex = 1.6, col = col_outlier, bg = col_outlier, lwd = 2)
          
          smart_text_labels(
            x_positions = tetraploid_stats_0$start_pos[focal_outliers_idx],
            y_positions = tetraploid_stats_0$fst[focal_outliers_idx],
            labels = tetraploid_stats_0$locus[focal_outliers_idx],
            pos = 3, cex = 0.8, col = col_outlier, offset = 0.5
          )
        }
      }
    }

    # Build legend text with threshold values
    legend_text <- c("Diploid", "Tetraploid", "Background", "Outlier (>95th)")
    # if(!is.na(fst_threshold_dip)) {
    #   legend_text[4] <- paste0("Outlier (>", round(fst_threshold_dip, 3), ")")
    # }
    
    legend("topleft",
           legend = legend_text,
           col = c(col_fst_dip, col_fst_tet, col_outlier),
           pch = c(19, 17, 21),
           pt.cex = c(1.4, 1.4, 1.5),  # *** INCREASED ***
           pt.bg = c(NA, NA, col_outlier),
           cex = 1.2,  # *** INCREASED ***
           bg = "white")
    
    # *** PANEL LABEL C ***
    add_panel_label("C")
    
    # === PANEL D: 4-FOLD FST ===
    plot(diploid_stats_4$start_pos, diploid_stats_4$fst, type = "n",
         xlim = global_xlim,
         ylim = c(0, y_max_fst),
         main = paste0("Fst (Low vs High) 4-fold sites"),
         ylab = "Fst_WC",
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.4,
         cex.lab = 1.5,
         cex.axis = 1.3)
    
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    add_outlier_lines(outlier_positions_4, 0, y_max_fst)
    
    abline(h = 0, col = "gray30", lwd = 0.5)
    # if(!is.na(fst_threshold_dip_4)) abline(h = fst_threshold_dip_4, lty = 2, col = col_fst_dip, lwd = 1.5)
    # if(!is.na(fst_threshold_tet_4)) abline(h = fst_threshold_tet_4, lty = 2, col = col_fst_tet, lwd = 1.5)

    ####
    #lines and labels of thresholds
    if(!is.na(fst_threshold_dip_4)) {abline(h = fst_threshold_dip_4, lty = 2, col = col_fst_dip, lwd = 1.5)
  
      # Add text label for diploid threshold
      text(global_xlim[2] - (global_xlim[2] - global_xlim[1]) * 0.02, 
          fst_threshold_dip_4,
          labels = paste0(" Dip 95th: ", round(fst_threshold_dip_4, 3)),
          pos = 3, cex = 1.0, col = col_fst_dip, offset = 0.3)
    }
    if(!is.na(fst_threshold_tet_4)) {abline(h = fst_threshold_tet_4, lty = 2, col = col_fst_tet, lwd = 1.5)
  
      # Add text label for diploid threshold
      text(global_xlim[1] + (global_xlim[2] - global_xlim[1]) * 0.02, 
          fst_threshold_tet_4,
          labels = paste0(" Tet 95th: ", round(fst_threshold_tet_4, 3)),
          pos = 3, cex = 1.0, col = col_fst_tet, offset = 0.3)
    }  
    ####
    
    rect(focal_gene_start, 0, focal_gene_end, y_max_fst,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_stats_4) > 0 && length(outliers_dip_4) > 0) {
      points(diploid_stats_4$start_pos[!outliers_dip_4], diploid_stats_4$fst[!outliers_dip_4], 
             pch = 19, cex = 1.4, col = col_fst_dip)
      
      if(sum(outliers_dip_4) > 0) {
        focal_outliers_idx <- outliers_dip_4 & 
          diploid_stats_4$start_pos >= focal_gene_start & 
          diploid_stats_4$start_pos <= focal_gene_end
        
        if(sum(focal_outliers_idx) > 0) {
          points(diploid_stats_4$start_pos[focal_outliers_idx], 
                 diploid_stats_4$fst[focal_outliers_idx], 
                 pch = 21, cex = 1.6, col = col_outlier, bg = col_outlier, lwd = 2)
          
          smart_text_labels(
            x_positions = diploid_stats_4$start_pos[focal_outliers_idx],
            y_positions = diploid_stats_4$fst[focal_outliers_idx],
            labels = diploid_stats_4$locus[focal_outliers_idx],
            pos = 3, cex = 0.8, col = col_outlier, offset = 0.5
          )
        }
      }
    }
    
    if(nrow(tetraploid_stats_4) > 0 && length(outliers_tet_4) > 0) {
      points(tetraploid_stats_4$start_pos[!outliers_tet_4], tetraploid_stats_4$fst[!outliers_tet_4], 
             pch = 17, cex = 1.4, col = col_fst_tet)
      
      if(sum(outliers_tet_4) > 0) {
        focal_outliers_idx <- outliers_tet_4 & 
          tetraploid_stats_4$start_pos >= focal_gene_start & 
          tetraploid_stats_4$start_pos <= focal_gene_end
        
        if(sum(focal_outliers_idx) > 0) {
          points(tetraploid_stats_4$start_pos[focal_outliers_idx], 
                 tetraploid_stats_4$fst[focal_outliers_idx], 
                 pch = 24, cex = 1.6, col = col_outlier, bg = col_outlier, lwd = 2)
          
          smart_text_labels(
            x_positions = tetraploid_stats_4$start_pos[focal_outliers_idx],
            y_positions = tetraploid_stats_4$fst[focal_outliers_idx],
            labels = tetraploid_stats_4$locus[focal_outliers_idx],
            pos = 3, cex = 0.8, col = col_outlier, offset = 0.5
          )
        }
      }
    }

    # Build legend text with threshold values
    legend_text <- c("Diploid", "Tetraploid", "Background", "Outlier (>95th)")
    # if(!is.na(fst_threshold_dip)) {
    #   legend_text[4] <- paste0("Outlier (>", round(fst_threshold_dip, 3), ")")
    # }
    
    legend("topleft",
           legend = legend_text,
           col = c(col_fst_dip, col_fst_tet, col_outlier),
           pch = c(19, 17, 21),
           pt.cex = c(1.4, 1.4,  1.5),  # *** INCREASED ***
           pt.bg = c(NA, NA, col_outlier),
           cex = 1.2,  # *** INCREASED ***
           bg = "white")
    
    # *** PANEL LABEL D ***
    add_panel_label("D")

    
    ##############################################
    # ROW 3: NUCLEOTIDE DIVERSITY
    ##############################################
    
    # === PANEL E: 0-FOLD PI ===
    plot(diploid_pi_low_0$pos1, diploid_pi_low_0$pi_low, type = "n",
         xlim = global_xlim,
         ylim = c(0, y_max_pi),
         main = "Nucleotide Diversity 0-fold sites",
         ylab = "Nucleotide Diversity",
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.6,
         cex.lab = 1.5,
         cex.axis = 1.3)
    
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    add_outlier_lines(outlier_positions_0, 0, y_max_pi)
    
    abline(h = 0, col = "gray30", lwd = 0.5)
    
    rect(focal_gene_start, 0, focal_gene_end, y_max_pi,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_pi_low_0) > 0) points(diploid_pi_low_0$pos1, diploid_pi_low_0$pi_low, pch=19, cex=1.4, col=col_pi_dip_low)
    if(nrow(diploid_pi_high_0) > 0) points(diploid_pi_high_0$pos1, diploid_pi_high_0$pi_high, pch=19, cex=1.4, col=col_pi_dip_high)
    if(nrow(tetraploid_pi_low_0) > 0) points(tetraploid_pi_low_0$pos1, tetraploid_pi_low_0$pi_low, pch=17, cex=1.4, col=col_pi_tet_low)
    if(nrow(tetraploid_pi_high_0) > 0) points(tetraploid_pi_high_0$pos1, tetraploid_pi_high_0$pi_high, pch=17, cex=1.4, col=col_pi_tet_high)
    
    # Highlight outliers - diploids
    if(length(outlier_positions_0) > 0) {
      for(df in list(diploid_pi_low_0, diploid_pi_high_0)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_0
          if(sum(matches) > 0) {
            pi_col <- grep("^pi_", names(df), value = TRUE)
            if(length(pi_col) > 0) {
              points(df$pos1[matches], df[[pi_col]][matches],
                     pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
      
      # Highlight outliers - tetraploids
      for(df in list(tetraploid_pi_low_0, tetraploid_pi_high_0)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_0
          if(sum(matches) > 0) {
            pi_col <- grep("^pi_", names(df), value = TRUE)
            if(length(pi_col) > 0) {
              points(df$pos1[matches], df[[pi_col]][matches],
                     pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
    }

    legend("topleft",
          legend = c("Dip Low", "Dip High", "Tet Low", "Tet High", "Fst outlier"), 
          col = c(col_pi_dip_low, col_pi_dip_high, col_pi_tet_low, col_pi_tet_high, col_outlier),
          pch = c(19, 19, 17, 17, 1),
          pt.cex = c(1.4, 1.4, 1.4, 1.4, 1.5),  # *** INCREASED ***
          lwd = c(NA, NA, NA, NA, 2.5),
          cex = 1.2,  # *** INCREASED ***
          bg = "white")

    # *** PANEL LABEL E ***
    add_panel_label("E")
    
    # === PANEL F: 4-FOLD PI ===
    plot(diploid_pi_low_4$pos1, diploid_pi_low_4$pi_low, type = "n",
         xlim = global_xlim,
         ylim = c(0, y_max_pi),
         main = "Nucleotide Diversity 4-fold sites",
         ylab = "Nucleotide Diversity",
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.6,
         cex.lab = 1.5,
         cex.axis = 1.3)
    
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    add_outlier_lines(outlier_positions_4, 0, y_max_pi)
    
    abline(h = 0, col = "gray30", lwd = 0.5)
    
    rect(focal_gene_start, 0, focal_gene_end, y_max_pi,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_pi_low_4) > 0) points(diploid_pi_low_4$pos1, diploid_pi_low_4$pi_low, pch=19, cex=1.4, col=col_pi_dip_low)
    if(nrow(diploid_pi_high_4) > 0) points(diploid_pi_high_4$pos1, diploid_pi_high_4$pi_high, pch=19, cex=1.4, col=col_pi_dip_high)
    if(nrow(tetraploid_pi_low_4) > 0) points(tetraploid_pi_low_4$pos1, tetraploid_pi_low_4$pi_low, pch=17, cex=1.4, col=col_pi_tet_low)
    if(nrow(tetraploid_pi_high_4) > 0) points(tetraploid_pi_high_4$pos1, tetraploid_pi_high_4$pi_high, pch=17, cex=1.4, col=col_pi_tet_high)
    
    # Highlight outliers - diploids
    if(length(outlier_positions_4) > 0) {
      for(df in list(diploid_pi_low_4, diploid_pi_high_4)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_4
          if(sum(matches) > 0) {
            pi_col <- grep("^pi_", names(df), value = TRUE)
            if(length(pi_col) > 0) {
              points(df$pos1[matches], df[[pi_col]][matches],
                     pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
      
      # Highlight outliers - tetraploids
      for(df in list(tetraploid_pi_low_4, tetraploid_pi_high_4)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_4
          if(sum(matches) > 0) {
            pi_col <- grep("^pi_", names(df), value = TRUE)
            if(length(pi_col) > 0) {
              points(df$pos1[matches], df[[pi_col]][matches],
                     pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
    }
    
    legend("topleft",
          legend = c("Dip Low", "Dip High", "Tet Low", "Tet High", "Fst outlier"), 
          col = c(col_pi_dip_low, col_pi_dip_high, col_pi_tet_low, col_pi_tet_high, col_outlier),
          pch = c(19, 19, 17, 17, 1),
          pt.cex = c(1.4, 1.4, 1.4, 1.4, 1.5),  # *** INCREASED ***
          lwd = c(NA, NA, NA, NA, 2.5),
          cex = 1.2,  # *** INCREASED ***
          bg = "white") 

    # *** PANEL LABEL F ***
    add_panel_label("F")

    
    ##############################################
    # ROW 4: TAJIMA'S D
    ##############################################
    
    # === PANEL G: 0-FOLD TAJIMA'S D ===
    plot(diploid_tajd_low_0$pos1, diploid_tajd_low_0$tajd_low, type = "n",
         xlim = global_xlim,
         ylim = c(y_min_tajd, y_max_tajd),
         main = "Tajima's D 0-fold sites",
         ylab = "Tajima's D",
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.6,
         cex.lab = 1.5,
         cex.axis = 1.3)
    
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    add_outlier_lines(outlier_positions_0, y_min_tajd, y_max_tajd)
    
    abline(h = 0, col = "black", lwd = 1.5)
    abline(h = -2, lty = 2, col = "blue", lwd = 1.5)
    abline(h = 2, lty = 2, col = "blue", lwd = 1.5)
    
    rect(focal_gene_start, y_min_tajd, focal_gene_end, y_max_tajd,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_tajd_low_0) > 0) points(diploid_tajd_low_0$pos1, diploid_tajd_low_0$tajd_low, pch=19, cex=1.4, col=col_tajd_dip_low)
    if(nrow(diploid_tajd_high_0) > 0) points(diploid_tajd_high_0$pos1, diploid_tajd_high_0$tajd_high, pch=19, cex=1.4, col=col_tajd_dip_high)
    if(nrow(tetraploid_tajd_low_0) > 0) points(tetraploid_tajd_low_0$pos1, tetraploid_tajd_low_0$tajd_low, pch=17, cex=1.4, col=col_tajd_tet_low)
    if(nrow(tetraploid_tajd_high_0) > 0) points(tetraploid_tajd_high_0$pos1, tetraploid_tajd_high_0$tajd_high, pch=17, cex=1.4, col=col_tajd_tet_high)
    
    # Highlight outliers - diploids
    if(length(outlier_positions_0) > 0) {
      for(df in list(diploid_tajd_low_0, diploid_tajd_high_0)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_0
          if(sum(matches) > 0) {
            tajd_col <- grep("^tajd_", names(df), value = TRUE)
            if(length(tajd_col) > 0) {
              points(df$pos1[matches], df[[tajd_col]][matches],
                     pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
      
      # Highlight outliers - tetraploids
      for(df in list(tetraploid_tajd_low_0, tetraploid_tajd_high_0)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_0
          if(sum(matches) > 0) {
            tajd_col <- grep("^tajd_", names(df), value = TRUE)
            if(length(tajd_col) > 0) {
              points(df$pos1[matches], df[[tajd_col]][matches],
                     pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
    }

    legend("topleft",
          legend = c("Neutral", "±2 threshold", "Fst outlier"),
          col = c("black", "blue", col_outlier), 
          pch = c(NA, NA, 1),
          lty = c(1, 2, NA),
          pt.cex = c(NA, NA, 2.0),  # *** INCREASED ***
          lwd = c(1.5, 1.5, 2.5),
          cex = 1.2,  # *** INCREASED ***
          bg = "white")
    
    
    # *** PANEL LABEL G ***
    add_panel_label("G")

    
    # === PANEL H: 4-FOLD TAJIMA'S D ===
    plot(diploid_tajd_low_4$pos1, diploid_tajd_low_4$tajd_low, type = "n",
         xlim = global_xlim,
         ylim = c(y_min_tajd, y_max_tajd),
         main = "Tajima's D 4-fold sites",
         ylab = "Tajima's D",
         xlab = "Genomic Position",
         xaxt = "n",
         cex.main = 1.6,
         cex.lab = 1.5,
         cex.axis = 1.3)
    
    axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = FALSE, big.mark = ","), 
         cex.axis = 1.3)
    
    add_outlier_lines(outlier_positions_4, y_min_tajd, y_max_tajd)
    
    abline(h = 0, col = "black", lwd = 1.5)
    abline(h = -2, lty = 2, col = "blue", lwd = 1.5)
    abline(h = 2, lty = 2, col = "blue", lwd = 1.5)
    
    rect(focal_gene_start, y_min_tajd, focal_gene_end, y_max_tajd,
         col = rgb(0.7, 0.7, 0.7, 0.3), border = NA)
    
    if(nrow(diploid_tajd_low_4) > 0) points(diploid_tajd_low_4$pos1, diploid_tajd_low_4$tajd_low, pch=19, cex=1.4, col=col_tajd_dip_low)
    if(nrow(diploid_tajd_high_4) > 0) points(diploid_tajd_high_4$pos1, diploid_tajd_high_4$tajd_high, pch=19, cex=1.4, col=col_tajd_dip_high)
    if(nrow(tetraploid_tajd_low_4) > 0) points(tetraploid_tajd_low_4$pos1, tetraploid_tajd_low_4$tajd_low, pch=17, cex=1.4, col=col_tajd_tet_low)
    if(nrow(tetraploid_tajd_high_4) > 0) points(tetraploid_tajd_high_4$pos1, tetraploid_tajd_high_4$tajd_high, pch=17, cex=1.4, col=col_tajd_tet_high)
    
    # Highlight outliers - diploids
    if(length(outlier_positions_4) > 0) {
      for(df in list(diploid_tajd_low_4, diploid_tajd_high_4)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_4
          if(sum(matches) > 0) {
            tajd_col <- grep("^tajd_", names(df), value = TRUE)
            if(length(tajd_col) > 0) {
              points(df$pos1[matches], df[[tajd_col]][matches],
                     pch = 1, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
      
      # Highlight outliers - tetraploids
      for(df in list(tetraploid_tajd_low_4, tetraploid_tajd_high_4)) {
        if(nrow(df) > 0) {
          matches <- df$pos1 %in% outlier_positions_4
          if(sum(matches) > 0) {
            tajd_col <- grep("^tajd_", names(df), value = TRUE)
            if(length(tajd_col) > 0) {
              points(df$pos1[matches], df[[tajd_col]][matches],
                     pch = 2, cex = 1.6, col = col_outlier, lwd = 2.5)
            }
          }
        }
      }
    }

    legend("topleft",
          legend = c("Neutral", "±2 threshold", "Fst outlier"),
          col = c("black", "blue", col_outlier), 
          pch = c(NA, NA, 1),
          lty = c(1, 2, NA),
          pt.cex = c(NA, NA, 2.0),  # *** INCREASED ***
          lwd = c(1.5, 1.5, 2.5),
          cex = 1.2,  # *** INCREASED ***
          bg = "white")
    
    
    # *** PANEL LABEL H ***
    add_panel_label("H")

    ##############################################
    # ROW 5: π₀/π₄ RATIOS
    ##############################################

    # Calculate all ratios
    #DIPLOID
    #low elevation
    focal_pi_0_dip_low <- mean(diploid_pi_low_0[diploid_pi_low_0$pos1 >= focal_gene_start & 
                                        diploid_pi_low_0$pos1 <= focal_gene_end, ]$pi_low, na.rm = TRUE)
    focal_pi_4_dip_low <- mean(diploid_pi_low_4[diploid_pi_low_4$pos1 >= focal_gene_start & 
                                        diploid_pi_low_4$pos1 <= focal_gene_end, ]$pi_low, na.rm = TRUE)
    focal_ratio_dip_low <- focal_pi_0_dip_low / focal_pi_4_dip_low

    #high elevation
    focal_pi_0_dip_high <- mean(diploid_pi_high_0[diploid_pi_high_0$pos1 >= focal_gene_start & 
                                            diploid_pi_high_0$pos1 <= focal_gene_end, ]$pi_high, na.rm = TRUE)
    focal_pi_4_dip_high <- mean(diploid_pi_high_4[diploid_pi_high_4$pos1 >= focal_gene_start & 
                                            diploid_pi_high_4$pos1 <= focal_gene_end, ]$pi_high, na.rm = TRUE)
    focal_ratio_dip_high <- focal_pi_0_dip_high / focal_pi_4_dip_high

    window_pi_0_dip_low <- mean(diploid_pi_low_0$pi_low, na.rm = TRUE)
    window_pi_4_dip_low <- mean(diploid_pi_low_4$pi_low, na.rm = TRUE)
    window_ratio_dip_low <- window_pi_0_dip_low / window_pi_4_dip_low

    window_pi_0_dip_high <- mean(diploid_pi_high_0$pi_high, na.rm = TRUE)
    window_pi_4_dip_high <- mean(diploid_pi_high_4$pi_high, na.rm = TRUE)
    window_ratio_dip_high <- window_pi_0_dip_high / window_pi_4_dip_high

    #TETRAPLOID
    # Calculate all ratios
    #low elevations
    focal_pi_0_tet_low <- mean(tetraploid_pi_low_0[tetraploid_pi_low_0$pos1 >= focal_gene_start & 
                                            tetraploid_pi_low_0$pos1 <= focal_gene_end, ]$pi_low, na.rm = TRUE)
    focal_pi_4_tet_low <- mean(tetraploid_pi_low_4[tetraploid_pi_low_4$pos1 >= focal_gene_start & 
                                            tetraploid_pi_low_4$pos1 <= focal_gene_end, ]$pi_low, na.rm = TRUE)
    focal_ratio_tet_low <- focal_pi_0_tet_low / focal_pi_4_tet_low

    #high elevations
    focal_pi_0_tet_high <- mean(tetraploid_pi_high_0[tetraploid_pi_high_0$pos1 >= focal_gene_start & 
                                            tetraploid_pi_high_0$pos1 <= focal_gene_end, ]$pi_high, na.rm = TRUE)
    focal_pi_4_tet_high <- mean(tetraploid_pi_high_4[tetraploid_pi_high_4$pos1 >= focal_gene_start & 
                                            tetraploid_pi_high_4$pos1 <= focal_gene_end, ]$pi_high, na.rm = TRUE)
    focal_ratio_tet_high <- focal_pi_0_tet_high / focal_pi_4_tet_high

    window_pi_0_tet_low <- mean(tetraploid_pi_low_0$pi_low, na.rm = TRUE)
    window_pi_4_tet_low <- mean(tetraploid_pi_low_4$pi_low, na.rm = TRUE)
    window_ratio_tet_low <- window_pi_0_tet_low / window_pi_4_tet_low

    window_pi_0_tet_high <- mean(tetraploid_pi_high_0$pi_high, na.rm = TRUE)
    window_pi_4_tet_high <- mean(tetraploid_pi_high_4$pi_high, na.rm = TRUE)
    window_ratio_tet_high <- window_pi_0_tet_high / window_pi_4_tet_high

    ##############################################
    # EXPORT π₀/π₄ RATIOS TO CSV
    ##############################################

    # Create summary data frame
    pi_ratio_summary <- data.frame(
      Gene = focal_gene_name,
      VARIA_ID = focal_id,
      
      # Diploid - Focal gene
      Dip_Low_Focal_pi0 = focal_pi_0_dip_low,
      Dip_Low_Focal_pi4 = focal_pi_4_dip_low,
      Dip_Low_Focal_ratio = focal_ratio_dip_low,
      
      Dip_High_Focal_pi0 = focal_pi_0_dip_high,
      Dip_High_Focal_pi4 = focal_pi_4_dip_high,
      Dip_High_Focal_ratio = focal_ratio_dip_high,
      
      # Diploid - Window
      Dip_Low_Window_pi0 = window_pi_0_dip_low,
      Dip_Low_Window_pi4 = window_pi_4_dip_low,
      Dip_Low_Window_ratio = window_ratio_dip_low,
      
      Dip_High_Window_pi0 = window_pi_0_dip_high,
      Dip_High_Window_pi4 = window_pi_4_dip_high,
      Dip_High_Window_ratio = window_ratio_dip_high,
      
      # Tetraploid - Focal gene
      Tet_Low_Focal_pi0 = focal_pi_0_tet_low,
      Tet_Low_Focal_pi4 = focal_pi_4_tet_low,
      Tet_Low_Focal_ratio = focal_ratio_tet_low,
      
      Tet_High_Focal_pi0 = focal_pi_0_tet_high,
      Tet_High_Focal_pi4 = focal_pi_4_tet_high,
      Tet_High_Focal_ratio = focal_ratio_tet_high,
      
      # Tetraploid - Window
      Tet_Low_Window_pi0 = window_pi_0_tet_low,
      Tet_Low_Window_pi4 = window_pi_4_tet_low,
      Tet_Low_Window_ratio = window_ratio_tet_low,
      
      Tet_High_Window_pi0 = window_pi_0_tet_high,
      Tet_High_Window_pi4 = window_pi_4_tet_high,
      Tet_High_Window_ratio = window_ratio_tet_high,
      
      stringsAsFactors = FALSE
    )

    # Save to CSV
    write.csv(pi_ratio_summary, 
              file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_pi_ratios.csv")),
              row.names = FALSE)

    message("π_0/π_4 ratios exported to: ", file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_pi_ratios.csv")))
    #diagnostics
    message("\n=== π_0/π_4 Ratio Diagnostics ===")
    message("Diploid Low - Focal: π_0=", round(focal_pi_0_dip_low, 4), " π_4=", round(focal_pi_4_dip_low, 4), " ratio=", round(focal_ratio_dip_low, 4))
    message("Diploid High - Focal: π_0=", round(focal_pi_0_dip_high, 4), " π_4=", round(focal_pi_4_dip_high, 4), " ratio=", round(focal_ratio_dip_high, 4))
    message("Tetraploid Low - Focal: π_0=", round(focal_pi_0_tet_low, 4), " π_4=", round(focal_pi_4_tet_low, 4), " ratio=", round(focal_ratio_tet_low, 4))
    message("Tetraploid High - Focal: π_0=", round(focal_pi_0_tet_high, 4), " π_4=", round(focal_pi_4_tet_high, 4), " ratio=", round(focal_ratio_tet_high, 4))


    #plotting check
    all_ratios <- c(focal_ratio_dip_low, window_ratio_dip_low, focal_ratio_dip_high, window_ratio_dip_high,
                focal_ratio_tet_low, window_ratio_tet_low, focal_ratio_tet_high, window_ratio_tet_high)

    if(all(is.na(all_ratios) | is.infinite(all_ratios) | is.nan(all_ratios))) {
      message("WARNING: All π₀/π₄ ratios are NA, Inf, or NaN - skipping panel")
      
      # Create empty panels with message
      plot.new()
      text(0.5, 0.5, "π_0/π_4 ratios unavailable\n(insufficient data)", cex = 1.5, col = "gray50")
      add_panel_label("I")
      
      plot.new()
      text(0.5, 0.5, "π_0/π_4 ratios unavailable\n(insufficient data)", cex = 1.5, col = "gray50")
      add_panel_label("J")
      
    } else {
      
      # Replace Inf/NaN with NA for plotting
      focal_ratio_dip_low <- ifelse(is.finite(focal_ratio_dip_low), focal_ratio_dip_low, NA)
      focal_ratio_dip_high <- ifelse(is.finite(focal_ratio_dip_high), focal_ratio_dip_high, NA)
      window_ratio_dip_low <- ifelse(is.finite(window_ratio_dip_low), window_ratio_dip_low, NA)
      window_ratio_dip_high <- ifelse(is.finite(window_ratio_dip_high), window_ratio_dip_high, NA)
      
      focal_ratio_tet_low <- ifelse(is.finite(focal_ratio_tet_low), focal_ratio_tet_low, NA)
      focal_ratio_tet_high <- ifelse(is.finite(focal_ratio_tet_high), focal_ratio_tet_high, NA)
      window_ratio_tet_low <- ifelse(is.finite(window_ratio_tet_low), window_ratio_tet_low, NA)
      window_ratio_tet_high <- ifelse(is.finite(window_ratio_tet_high), window_ratio_tet_high, NA)
      
      

    # === PANEL I: DIPLOID RATIOS ===
    # Create matrix for grouped barplot
    ratio_matrix_dip <- matrix(
    c(focal_ratio_dip_low, window_ratio_dip_low,      # Column 1: Low elevation
        focal_ratio_dip_high, window_ratio_dip_high),   # Column 2: High elevation
    nrow = 2,
    byrow = FALSE
    )
    colnames(ratio_matrix_dip) <- c("Low Elev", "High Elev")

    # Calculate y-limit, handling NA values
    y_max_dip <- max(ratio_matrix_dip, 1, na.rm = TRUE)
    if(!is.finite(y_max_dip)) y_max_dip <- 1.5

    # Create grouped bar plot
    bp_dip <- barplot(ratio_matrix_dip,
            beside = TRUE,
            ylim = c(0, y_max_dip * 1.3),
            main = "π_0/π_4 Ratio - Diploid",
            ylab = "π_0/π_4",
            col = c(col_focal_ratio, col_window_ratio,    # Low elev colors
                    col_focal_ratio, col_window_ratio), # High elev colors
            cex.main = 1.6,
            cex.lab = 1.5,
            cex.axis = 1.3,
            cex.names = 1.2,
            legend.text = c(paste0(focal_gene_name, " (focal)"), "5-gene window"),
            args.legend = list(x = "topleft", cex = 1.0, bty = "n"))

    abline(h = 1.0, lty = 2, col = "gray40", lwd = 2)

    # Add text labels on bars
    # text(bp_dip, ratio_matrix_dip + 0.02, 
    #     labels = round(ratio_matrix_dip, 3), 
    #     pos = 3, cex = 0.9, font = 2)
      # Add text labels (only for finite values)
    for(i in 1:ncol(ratio_matrix_dip)) {
      for(j in 1:nrow(ratio_matrix_dip)) {
        if(is.finite(ratio_matrix_dip[j, i])) {
          text(bp_dip[j, i], ratio_matrix_dip[j, i] + 0.02, 
              labels = round(ratio_matrix_dip[j, i], 3), 
              pos = 3, cex = 0.9, font = 2)
        } else {
          text(bp_dip[j, i], 0.05, labels = "NA", pos = 3, cex = 0.9, col = "gray50")
        }
      }
    }

    add_panel_label("I")

    # === PANEL J: TETRAPLOID RATIOS ===

    # Create matrix
    ratio_matrix_tet <- matrix(
    c(focal_ratio_tet_low, window_ratio_tet_low,
        focal_ratio_tet_high, window_ratio_tet_high),
    nrow = 2, byrow = FALSE
    )
    colnames(ratio_matrix_tet) <- c("Low Elev", "High Elev")

    y_max_tet <- max(ratio_matrix_tet, 1, na.rm = TRUE)
    if(!is.finite(y_max_tet)) y_max_tet <- 1.5

    bp_tet <- barplot(ratio_matrix_tet,
            beside = TRUE,
            ylim = c(0, y_max_tet * 1.3),
            main = "π_0/π_4 Ratio - Tetraploid",
            ylab = "π_0/π_4",
            col = c(col_focal_ratio, col_window_ratio,    # Low elev colors
                    col_focal_ratio, col_window_ratio), # High elev colors
            cex.main = 1.6,
            cex.lab = 1.5,
            cex.axis = 1.3,
            cex.names = 1.2,
            legend.text = c(paste0(focal_gene_name, " (focal)"), "5-gene window"),
            args.legend = list(x = "topleft", cex = 1.0, bty = "n"))

    abline(h = 1.0, lty = 2, col = "gray40", lwd = 2)

    # text(bp_tet, ratio_matrix_tet + 0.02, 
    #     labels = round(ratio_matrix_tet, 3), 
    #     pos = 3, cex = 0.9, font = 2)
    
    for(i in 1:ncol(ratio_matrix_tet)) {
      for(j in 1:nrow(ratio_matrix_tet)) {
        if(is.finite(ratio_matrix_tet[j, i])) {
          text(bp_tet[j, i], ratio_matrix_tet[j, i] + 0.02, 
              labels = round(ratio_matrix_tet[j, i], 3), 
              pos = 3, cex = 0.9, font = 2)
        } else {
          text(bp_tet[j, i], 0.05, labels = "NA", pos = 3, cex = 0.9, col = "gray50")
        }
      }
    }

    add_panel_label("J")

    #ROW 6 - pi0 and pi4 boxplots
    ##############################################
    # DATA PREPARATION FOR PANELS K & L (Boxplots)
    ##############################################

    # Check that focal gene coordinates are defined
    if(!exists("focal_gene_start") || !exists("focal_gene_end")) {
      message("ERROR: focal_gene_start or focal_gene_end not found")
      plot.new()
      add_panel_label("K")
      plot.new()
      add_panel_label("L")
    } else {

    message("\n=== Preparing boxplot data ===")
    message("Focal gene: ", focal_gene_start, " to ", focal_gene_end)

    # Extract all π values (not just means) for each site class and region

    # === DIPLOIDS ===

    # Focal gene - Low elevation - 0-fold 
    pi0_dip_low_focal_values <- diploid_pi_low_0$pi_low[
      diploid_pi_low_0$pos1 >= focal_gene_start & 
      diploid_pi_low_0$pos1 <= focal_gene_end
    ]
    pi0_dip_low_focal_values <- pi0_dip_low_focal_values[!is.na(pi0_dip_low_focal_values)]

    # Focal gene - Low elevation - 4-fold
    pi4_dip_low_focal_values <- diploid_pi_low_4$pi_low[
      diploid_pi_low_4$pos1 >= focal_gene_start & 
      diploid_pi_low_4$pos1 <= focal_gene_end
    ]
    pi4_dip_low_focal_values <- pi4_dip_low_focal_values[!is.na(pi4_dip_low_focal_values)]

    # Gene window - Low elevation - 0-fold (entire window)
    pi0_dip_low_window_values <- diploid_pi_low_0$pi_low
    pi0_dip_low_window_values <- pi0_dip_low_window_values[!is.na(pi0_dip_low_window_values)]

    # Gene window - Low elevation - 4-fold
    pi4_dip_low_window_values <- diploid_pi_low_4$pi_low
    pi4_dip_low_window_values <- pi4_dip_low_window_values[!is.na(pi4_dip_low_window_values)]

    # Focal gene - High elevation - 0-fold
    pi0_dip_high_focal_values <- diploid_pi_high_0$pi_high[
      diploid_pi_high_0$pos1 >= focal_gene_start & 
      diploid_pi_high_0$pos1 <= focal_gene_end
    ]
    pi0_dip_high_focal_values <- pi0_dip_high_focal_values[!is.na(pi0_dip_high_focal_values)]

    # Focal gene - High elevation - 4-fold
    pi4_dip_high_focal_values <- diploid_pi_high_4$pi_high[
      diploid_pi_high_4$pos1 >= focal_gene_start & 
      diploid_pi_high_4$pos1 <= focal_gene_end
    ]
    pi4_dip_high_focal_values <- pi4_dip_high_focal_values[!is.na(pi4_dip_high_focal_values)]

    # Gene window - High elevation - 0-fold
    pi0_dip_high_window_values <- diploid_pi_high_0$pi_high
    pi0_dip_high_window_values <- pi0_dip_high_window_values[!is.na(pi0_dip_high_window_values)]

    # Gene window - High elevation - 4-fold
    pi4_dip_high_window_values <- diploid_pi_high_4$pi_high
    pi4_dip_high_window_values <- pi4_dip_high_window_values[!is.na(pi4_dip_high_window_values)]


    # === TETRAPLOIDS ===

    pi0_tet_low_focal_values <- tetraploid_pi_low_0$pi_low[
      tetraploid_pi_low_0$pos1 >= focal_gene_start & 
      tetraploid_pi_low_0$pos1 <= focal_gene_end
    ]
    pi0_tet_low_focal_values <- pi0_tet_low_focal_values[!is.na(pi0_tet_low_focal_values)]

    pi4_tet_low_focal_values <- tetraploid_pi_low_4$pi_low[
      tetraploid_pi_low_4$pos1 >= focal_gene_start & 
      tetraploid_pi_low_4$pos1 <= focal_gene_end
    ]
    pi4_tet_low_focal_values <- pi4_tet_low_focal_values[!is.na(pi4_tet_low_focal_values)]

    pi0_tet_low_window_values <- tetraploid_pi_low_0$pi_low
    pi0_tet_low_window_values <- pi0_tet_low_window_values[!is.na(pi0_tet_low_window_values)]

    pi4_tet_low_window_values <- tetraploid_pi_low_4$pi_low
    pi4_tet_low_window_values <- pi4_tet_low_window_values[!is.na(pi4_tet_low_window_values)]

    pi0_tet_high_focal_values <- tetraploid_pi_high_0$pi_high[
      tetraploid_pi_high_0$pos1 >= focal_gene_start & 
      tetraploid_pi_high_0$pos1 <= focal_gene_end
    ]
    pi0_tet_high_focal_values <- pi0_tet_high_focal_values[!is.na(pi0_tet_high_focal_values)]

    pi4_tet_high_focal_values <- tetraploid_pi_high_4$pi_high[
      tetraploid_pi_high_4$pos1 >= focal_gene_start & 
      tetraploid_pi_high_4$pos1 <= focal_gene_end
    ]
    pi4_tet_high_focal_values <- pi4_tet_high_focal_values[!is.na(pi4_tet_high_focal_values)]

    pi0_tet_high_window_values <- tetraploid_pi_high_0$pi_high
    pi0_tet_high_window_values <- pi0_tet_high_window_values[!is.na(pi0_tet_high_window_values)]

    pi4_tet_high_window_values <- tetraploid_pi_high_4$pi_high
    pi4_tet_high_window_values <- pi4_tet_high_window_values[!is.na(pi4_tet_high_window_values)]


    # Diagnostic message
    message("\n=== Boxplot Data Summary ===")
    message("Diploid Low - Focal: ", length(pi0_dip_low_focal_values), " π_0 sites, ",
            length(pi4_dip_low_focal_values), " π_4 sites")
    message("Diploid Low - Window: ", length(pi0_dip_low_window_values), " π_0 sites, ",
            length(pi4_dip_low_window_values), " π_4 sites")
    message("Tetraploid Low - Focal: ", length(pi0_tet_low_focal_values), " π_0 sites, ",
            length(pi4_tet_low_focal_values), " π_4 sites")

    ##############################################
    # PANEL K: DIPLOID π₀ AND π₄ DISTRIBUTIONS
    ##############################################

    if(all(c(length(pi0_dip_low_focal_values), 
            length(pi4_dip_low_focal_values),
            length(pi0_dip_high_focal_values),
            length(pi4_dip_high_focal_values)) > 0)) {
    
      boxplot_data_dip <- list(
        pi0_dip_low_focal_values,
        pi4_dip_low_focal_values,
        pi0_dip_low_window_values,
        pi4_dip_low_window_values,
        pi0_dip_high_focal_values,
        pi4_dip_high_focal_values,
        pi0_dip_high_window_values,
        pi4_dip_high_window_values
      )
      
      box_colors <- rep(c("#E41A1C", "#377EB8"), times = 4)
      
      bp <- boxplot(boxplot_data_dip,
                    col = box_colors,
                    border = "black",
                    at = c(1, 2,  4.5, 5.5,  8, 9,  11.5, 12.5),
                    names = rep("", 8),
                    main = "π_0 and π_4 Values - Diploid",
                    ylab = "Nucleotide Diversity",
                    xlab = "",
                    ylim = c(0, max(unlist(boxplot_data_dip), na.rm=TRUE) * 1.1),
                    las = 1,
                    xaxt = "n",
                    cex.main = 1.6,
                    cex.lab = 1.5,
                    cex.axis = 1.0,
                    outline = TRUE,
                    pch = 20,
                    cex = 0.5)
      
      text(x = c(1, 2, 4.5, 5.5, 8, 9, 11.5, 12.5),
          y = par("usr")[3] - 0.04 * diff(par("usr")[3:4]),
          labels = rep(c("π_0", "π_4"), times = 4),
          xpd = TRUE, cex = 0.9, font = 1)
      
      text(x = c(1.5, 5, 8.5, 12),
          y = par("usr")[3] - 0.08 * diff(par("usr")[3:4]),
          labels = c("Focal", "Window", "Focal", "Window"),
          xpd = TRUE, cex = 0.9, font = 1)
      
      text(x = c(2.75, 10),
          y = par("usr")[3] - 0.12 * diff(par("usr")[3:4]),
          labels = c("Low Elevation", "High Elevation"),
          xpd = TRUE, cex = 1.0, font = 2)
      
      text(x = 1.5, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", focal_ratio_dip_low),
          cex = 0.85, font = 2)
      
      text(x = 5, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", window_ratio_dip_low),
          cex = 0.85, font = 2)
      
      text(x = 8.5, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", focal_ratio_dip_high),
          cex = 0.85, font = 2)
      
      text(x = 12, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", window_ratio_dip_high),
          cex = 0.85, font = 2)
      
      add_panel_label("K")
      
      legend("top",
            legend = c("π_0", "π_4"),
            fill = c("#E41A1C", "#377EB8"),
            border = "black",
            bty = "n",
            cex = 0.9)
      
    } else {
      plot.new()
      text(0.5, 0.5, 
          "Insufficient polymorphism data\nfor boxplot visualization",
          cex = 1.2, col = "gray50")
      add_panel_label("K")
    }

    ##############################################
    # PANEL L: TETRAPLOID π₀ AND π₄ DISTRIBUTIONS
    ##############################################

    if(all(c(length(pi0_tet_low_focal_values), 
            length(pi4_tet_low_focal_values),
            length(pi0_tet_high_focal_values),
            length(pi4_tet_high_focal_values)) > 0)) {
      
      boxplot_data_tet <- list(
        pi0_tet_low_focal_values,
        pi4_tet_low_focal_values,
        pi0_tet_low_window_values,
        pi4_tet_low_window_values,
        pi0_tet_high_focal_values,
        pi4_tet_high_focal_values,
        pi0_tet_high_window_values,
        pi4_tet_high_window_values
      )
      
      box_colors <- rep(c("#E41A1C", "#377EB8"), times = 4)
      
      bp <- boxplot(boxplot_data_tet,
                    col = box_colors,
                    border = "black",
                    at = c(1, 2,  4.5, 5.5,  8, 9,  11.5, 12.5),
                    names = rep("", 8),
                    main = "π_0 and π_4 Values - Tetraploid",
                    ylab = "Nucleotide Diversity",
                    xlab = "",
                    ylim = c(0, max(unlist(boxplot_data_tet), na.rm=TRUE) * 1.1),
                    las = 1,
                    xaxt = "n",
                    cex.lab = 1.5,
                    cex.main = 1.6,
                    cex.axis = 1.0,
                    outline = TRUE,
                    pch = 20,
                    cex = 0.5)
      
      text(x = c(1, 2, 4.5, 5.5, 8, 9, 11.5, 12.5),
          y = par("usr")[3] - 0.04 * diff(par("usr")[3:4]),
          labels = rep(c("π_0", "π_4"), times = 4),
          xpd = TRUE, cex = 0.9, font = 1)
      
      text(x = c(1.5, 5, 8.5, 12),
          y = par("usr")[3] - 0.08 * diff(par("usr")[3:4]),
          labels = c("Focal", "Window", "Focal", "Window"),
          xpd = TRUE, cex = 0.9, font = 1)
      
      text(x = c(2.75, 10),
          y = par("usr")[3] - 0.12 * diff(par("usr")[3:4]),
          labels = c("Low Elevation", "High Elevation"),
          xpd = TRUE, cex = 1.0, font = 2)
      
      
      text(x = 1.5, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", focal_ratio_tet_low),
          cex = 0.85, font = 2)
      
      text(x = 5, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", window_ratio_tet_low),
          cex = 0.85, font = 2)
      
      text(x = 8.5, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", focal_ratio_tet_high),
          cex = 0.85, font = 2)
      
      text(x = 12, y = par("usr")[4] * 0.95,
          labels = sprintf("π_0/π_4 = %.3f", window_ratio_tet_high),
          cex = 0.85, font = 2)
      
      add_panel_label("L")
      
      legend("top",
            legend = c("π_0", "π_4"),
            fill = c("#E41A1C", "#377EB8"),
            border = "black",
            bty = "n",
            cex = 0.9)
      
    } else {
      plot.new()
      text(0.5, 0.5, 
          "Insufficient polymorphism data\nfor boxplot visualization",
          cex = 1.2, col = "gray50")
      add_panel_label("L")
    }
    
    }
 
    ##############################################
    # OVERALL TITLE
    ##############################################
    
    mtext(paste0("0-fold vs 4-fold Comparison: ", focal_gene_name, " (", focal_id, ")"), 
          outer = TRUE, cex = 1.8, font = 2, line = 0.5)
    
    dev.off()
    
    message("Comparison plot saved: ", png_file)
    
  }}, error = function(e) {
    if(dev.cur() > 1) dev.off()
    message("ERROR creating comparison plot: ", e$message)
    return(NULL)
  })
}


################################################################################
# FUNCTION: CREATE COMPREHENSIVE SUMMARY TABLE
# Combines Fst, π, Tajima's D, and outlier status for all sites
################################################################################

create_comprehensive_summary <- function(results_0fold, results_4fold, degeneracy_suffix) {
  
  message("\nCreating comprehensive summary table...")
  
  ##############################################
  # HELPER FUNCTION: MERGE ALL STATISTICS
  ##############################################
  
  create_summary_for_degeneracy <- function(results, degeneracy_label) {
    
    # Extract data
    diploid_stats <- results$diploid_stats
    tetraploid_stats <- results$tetraploid_stats
    diploid_pi_low <- results$diploid_pi_low
    diploid_pi_high <- results$diploid_pi_high
    tetraploid_pi_low <- results$tetraploid_pi_low
    tetraploid_pi_high <- results$tetraploid_pi_high
    diploid_tajd_low <- results$diploid_tajd_low
    diploid_tajd_high <- results$diploid_tajd_high
    tetraploid_tajd_low <- results$tetraploid_tajd_low
    tetraploid_tajd_high <- results$tetraploid_tajd_high
    outliers_dip <- results$outliers_dip
    outliers_tet <- results$outliers_tet
    fst_threshold_dip <- results$fst_threshold_dip
    fst_threshold_tet <- results$fst_threshold_tet
    
    # Initialize empty list to collect rows
    summary_rows <- list()
    
    ##############################################
    # PROCESS DIPLOID DATA
    ##############################################
    
    if(nrow(diploid_stats) > 0) {
      
      for(i in 1:nrow(diploid_stats)) {
        
        locus <- diploid_stats$locus[i]
        start_pos <- diploid_stats$start_pos[i]
        fst <- diploid_stats$fst[i]
        is_outlier <- if(length(outliers_dip) > 0) outliers_dip[i] else FALSE
        in_focal <- start_pos >= focal_gene_start & start_pos <= focal_gene_end
        
        # Get π values for this position
        pi_low <- NA
        pi_high <- NA
        
        if(nrow(diploid_pi_low) > 0) {
          match_idx <- which(diploid_pi_low$pos1 == start_pos)
          if(length(match_idx) > 0) pi_low <- diploid_pi_low$pi_low[match_idx[1]]
        }
        
        if(nrow(diploid_pi_high) > 0) {
          match_idx <- which(diploid_pi_high$pos1 == start_pos)
          if(length(match_idx) > 0) pi_high <- diploid_pi_high$pi_high[match_idx[1]]
        }
        
        # Get Tajima's D values for this position
        tajd_low <- NA
        tajd_high <- NA
        
        if(nrow(diploid_tajd_low) > 0) {
          match_idx <- which(diploid_tajd_low$pos1 == start_pos)
          if(length(match_idx) > 0) tajd_low <- diploid_tajd_low$tajd_low[match_idx[1]]
        }
        
        if(nrow(diploid_tajd_high) > 0) {
          match_idx <- which(diploid_tajd_high$pos1 == start_pos)
          if(length(match_idx) > 0) tajd_high <- diploid_tajd_high$tajd_high[match_idx[1]]
        }
        
        # Create row
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          locus = locus,
          degeneracy = degeneracy_label,
          contig = sub(":.*", "", locus),
          start_pos = start_pos,
          gene = focal_gene_name,
          varia_id = focal_id,
          ploidy = "diploid",
          comparison = "low_vs_high",
          fst = fst,
          fst_95th_threshold = fst_threshold_dip,
          is_outlier = is_outlier,
          in_focal_gene = in_focal,
          pi_low = pi_low,
          pi_high = pi_high,
          pi_ratio = if(!is.na(pi_low) && !is.na(pi_high) && pi_high > 0) pi_low / pi_high else NA,
          tajd_low = tajd_low,
          tajd_high = tajd_high,
          tajd_diff = if(!is.na(tajd_low) && !is.na(tajd_high)) tajd_low - tajd_high else NA,
          stringsAsFactors = FALSE
        )
      }
    }
    
    ##############################################
    # PROCESS TETRAPLOID DATA
    ##############################################
    
    if(nrow(tetraploid_stats) > 0) {
      
      for(i in 1:nrow(tetraploid_stats)) {
        
        locus <- tetraploid_stats$locus[i]
        start_pos <- tetraploid_stats$start_pos[i]
        fst <- tetraploid_stats$fst[i]
        is_outlier <- if(length(outliers_tet) > 0) outliers_tet[i] else FALSE
        in_focal <- start_pos >= focal_gene_start & start_pos <= focal_gene_end
        
        # Get π values
        pi_low <- NA
        pi_high <- NA
        
        if(nrow(tetraploid_pi_low) > 0) {
          match_idx <- which(tetraploid_pi_low$pos1 == start_pos)
          if(length(match_idx) > 0) pi_low <- tetraploid_pi_low$pi_low[match_idx[1]]
        }
        
        if(nrow(tetraploid_pi_high) > 0) {
          match_idx <- which(tetraploid_pi_high$pos1 == start_pos)
          if(length(match_idx) > 0) pi_high <- tetraploid_pi_high$pi_high[match_idx[1]]
        }
        
        # Get Tajima's D values
        tajd_low <- NA
        tajd_high <- NA
        
        if(nrow(tetraploid_tajd_low) > 0) {
          match_idx <- which(tetraploid_tajd_low$pos1 == start_pos)
          if(length(match_idx) > 0) tajd_low <- tetraploid_tajd_low$tajd_low[match_idx[1]]
        }
        
        if(nrow(tetraploid_tajd_high) > 0) {
          match_idx <- which(tetraploid_tajd_high$pos1 == start_pos)
          if(length(match_idx) > 0) tajd_high <- tetraploid_tajd_high$tajd_high[match_idx[1]]
        }
        
        # Create row
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          locus = locus,
          degeneracy = degeneracy_label,
          contig = sub(":.*", "", locus),
          start_pos = start_pos,
          gene = focal_gene_name,
          varia_id = focal_id,
          ploidy = "tetraploid",
          comparison = "low_vs_high",
          fst = fst,
          fst_95th_threshold = fst_threshold_tet,
          is_outlier = is_outlier,
          in_focal_gene = in_focal,
          pi_low = pi_low,
          pi_high = pi_high,
          pi_ratio = if(!is.na(pi_low) && !is.na(pi_high) && pi_high > 0) pi_low / pi_high else NA,
          tajd_low = tajd_low,
          tajd_high = tajd_high,
          tajd_diff = if(!is.na(tajd_low) && !is.na(tajd_high)) tajd_low - tajd_high else NA,
          stringsAsFactors = FALSE
        )
      }
    }
    
    # Combine all rows
    if(length(summary_rows) > 0) {
      summary_df <- do.call(rbind, summary_rows)
      return(summary_df)
    } else {
      return(data.frame())
    }
  }
  
  ##############################################
  # CREATE SUMMARIES FOR BOTH DEGENERACIES
  ##############################################
  
  summary_0fold <- create_summary_for_degeneracy(results_0fold, "0fold")
  summary_4fold <- create_summary_for_degeneracy(results_4fold, "4fold")
  
  # Combine
  comprehensive_summary <- rbind(summary_0fold, summary_4fold)
  
  # Sort by position and degeneracy
  if(nrow(comprehensive_summary) > 0) {
    comprehensive_summary <- comprehensive_summary[order(
      comprehensive_summary$degeneracy,
      comprehensive_summary$ploidy,
      comprehensive_summary$start_pos
    ), ]
    
    # Add row names
    rownames(comprehensive_summary) <- NULL
  }
  
  ##############################################
  # SAVE TO FILE
  ##############################################
  
  output_file <- file.path(outdir, paste0(focal_gene_name, "_", focal_id, "_comprehensive_summary.csv"))
  
  write.csv(comprehensive_summary, output_file, row.names = FALSE, quote = FALSE)
  
  message("Comprehensive summary saved: ", output_file)
  message("  Total sites: ", nrow(comprehensive_summary))
  message("  0-fold sites: ", sum(comprehensive_summary$degeneracy == "0fold"))
  message("  4-fold sites: ", sum(comprehensive_summary$degeneracy == "4fold"))
  message("  Outliers: ", sum(comprehensive_summary$is_outlier, na.rm = TRUE))
  message("  Outliers in focal gene: ", sum(comprehensive_summary$is_outlier & comprehensive_summary$in_focal_gene, na.rm = TRUE))
  
  return(comprehensive_summary)
}


################################################################################
# MAIN ANALYSIS
################################################################################

##############################################
# LOAD GENE INFO (SHARED FOR BOTH ANALYSES)
##############################################

message("Loading gene information...")

gene_coord_colnames <- c("contig", "gene_start", "gene_end", "gene", "varia_id")

focal_gene <- read.table(focal_gene_info, header=FALSE, sep="\t", 
                         stringsAsFactors=FALSE, col.names=gene_coord_colnames)
selected_geneic_window <- read.table(selected_gene_of_interest_window, header=FALSE, sep="\t",
                                     stringsAsFactors=FALSE, col.names=gene_coord_colnames)
flank_master_list <- read.table(master_flanking_windows, header=FALSE, sep="\t",
                                stringsAsFactors=FALSE, col.names=gene_coord_colnames)

focal_row <- selected_geneic_window[selected_geneic_window$varia_id %in% focal_gene$varia_id, ]

if (nrow(focal_row) == 0) {
  stop("No focal gene found in window!")
} else if (nrow(focal_row) > 1) {
  message("Multiple focal genes found, using first one")
  focal_row <- focal_row[1, ]
}

focal_gene_start <- focal_row$gene_start[1]
focal_gene_end <- focal_row$gene_end[1]
focal_id <- focal_row$varia_id[1]
focal_gene_name <- focal_row$gene[1]
focal_contig <- focal_row$contig[1]

message("\nFocal gene: ", focal_gene_name, " (", focal_id, ")")
message("Position: ", focal_gene_start, "-", focal_gene_end)

# Create gene labels for plotting
gene_labels <- selected_geneic_window[order(selected_geneic_window$gene_start), ]

##############################################
# RUN ANALYSES
##############################################

# Define coverage files - same for both 0fold and 4fold
coverage_0fold_file <- file.path(coverage_dir, paste0(focal_id,"_", focal_contig, "_", focal_gene_name,"_coverage.tsv"))
coverage_4fold_file <- file.path(coverage_dir, paste0(focal_id,"_", focal_contig, "_", focal_gene_name,"_coverage.tsv"))

# Analyze 0-fold sites
results_0fold <- analyze_degeneracy(
  piawka_0fold_stats, 
  combined_bed_0fold, 
  "0fold",
  coverage_file = coverage_0fold_file
)

# Analyze 4-fold sites
results_4fold <- analyze_degeneracy(
  piawka_4fold_stats, 
  combined_bed_4fold, 
  "4fold",
  coverage_file = coverage_4fold_file
)

##############################################
# CREATE INDIVIDUAL PLOTS
##############################################

plot_degeneracy(results_0fold, "0fold")
plot_degeneracy(results_4fold, "4fold")

##############################################
# CREATE SIDE-BY-SIDE COMPARISON PLOT
##############################################

plot_comparison_side_by_side(results_0fold, results_4fold)

#CREAT COMBINED SUMMARY
comprehensive_summary <- create_comprehensive_summary(results_0fold, results_4fold, "combined")

################################################################################
# SUMMARY STATISTICS TABLE
################################################################################

message("\nCreating summary statistics table...")

create_summary <- function(data, label, ploidy_label) {
  if(nrow(data) == 0) {
    return(data.frame(
      Dataset = paste(ploidy_label, label),
      N_sites = 0,
      Mean_Fst = NA,
      Median_Fst = NA,
      SD_Fst = NA,
      Max_Fst = NA,
      N_outliers = 0,
      Fst_95th_threshold = NA,
      # N_samples_low = NA,
      # N_samples_high = NA,
      Focal_mean = NA,
      Flanking_mean = NA
    ))
  }
  
  data.frame(
    Dataset = paste(ploidy_label, label),
    N_sites = nrow(data),
    Mean_Fst = mean(data$fst, na.rm=TRUE),
    Median_Fst = median(data$fst, na.rm=TRUE),
    SD_Fst = sd(data$fst, na.rm=TRUE),
    Max_Fst = max(data$fst, na.rm=TRUE),
    N_outliers = sum(data$is_outlier, na.rm=TRUE),
    Fst_95th_threshold = unique(data$fst_95th_threshold)[1],
    # N_samples_low = unique(data$n_samples_low)[1],
    # N_samples_high = unique(data$n_samples_high)[1],
    Focal_mean = mean(data$fst[data$in_focal_gene], na.rm=TRUE),
    Flanking_mean = mean(data$fst[!data$in_focal_gene], na.rm=TRUE)
  )
}

# Create summary for each dataset
summary_list <- list()

# 0-fold diploid
if(!is.null(results_0fold$all_fst_data) && nrow(results_0fold$all_fst_data) > 0) {
  summary_list[[1]] <- create_summary(
    results_0fold$all_fst_data[results_0fold$all_fst_data$ploidy == "diploid", ], 
    "0-fold", "Diploid"
  )
} else {
  summary_list[[1]] <- data.frame(
    Dataset = "Diploid 0-fold",
    N_sites = 0, Mean_Fst = NA, Median_Fst = NA, SD_Fst = NA, Max_Fst = NA,
    N_outliers = 0, Fst_95th_threshold = NA, 
    #N_samples_low = NA, N_samples_high = NA,
    Focal_mean = NA, Flanking_mean = NA
  )
}

# 4-fold diploid
if(!is.null(results_4fold$all_fst_data) && nrow(results_4fold$all_fst_data) > 0) {
  summary_list[[2]] <- create_summary(
    results_4fold$all_fst_data[results_4fold$all_fst_data$ploidy == "diploid", ], 
    "4-fold", "Diploid"
  )
} else {
  summary_list[[2]] <- data.frame(
    Dataset = "Diploid 4-fold",
    N_sites = 0, Mean_Fst = NA, Median_Fst = NA, SD_Fst = NA, Max_Fst = NA,
    N_outliers = 0, Fst_95th_threshold = NA, 
    #N_samples_low = NA, N_samples_high = NA,
    Focal_mean = NA, Flanking_mean = NA
  )
}

# 0-fold tetraploid
if(!is.null(results_0fold$all_fst_data) && nrow(results_0fold$all_fst_data) > 0) {
  summary_list[[3]] <- create_summary(
    results_0fold$all_fst_data[results_0fold$all_fst_data$ploidy == "tetraploid", ], 
    "0-fold", "Tetraploid"
  )
} else {
  summary_list[[3]] <- data.frame(
    Dataset = "Tetraploid 0-fold",
    N_sites = 0, Mean_Fst = NA, Median_Fst = NA, SD_Fst = NA, Max_Fst = NA,
    N_outliers = 0, Fst_95th_threshold = NA, 
    #N_samples_low = NA, N_samples_high = NA,
    Focal_mean = NA, Flanking_mean = NA
  )
}

# 4-fold tetraploid
if(!is.null(results_4fold$all_fst_data) && nrow(results_4fold$all_fst_data) > 0) {
  summary_list[[4]] <- create_summary(
    results_4fold$all_fst_data[results_4fold$all_fst_data$ploidy == "tetraploid", ], 
    "4-fold", "Tetraploid"
  )
} else {
  summary_list[[4]] <- data.frame(
    Dataset = "Tetraploid 4-fold",
    N_sites = 0, Mean_Fst = NA, Median_Fst = NA, SD_Fst = NA, Max_Fst = NA,
    N_outliers = 0, Fst_95th_threshold = NA, 
    #N_samples_low = NA, N_samples_high = NA,
    Focal_mean = NA, Flanking_mean = NA
  )
}

# Combine all summaries
summary_stats <- do.call(rbind, summary_list)

# Round numeric columns
summary_stats[, 3:10] <- round(summary_stats[, 3:10], 4)

write.csv(summary_stats, 
          file.path(outdir, paste0(focal_gene_name, "_", focal_id,"_summary_statistics_0fold_vs_4fold.csv")),
          row.names = FALSE)

message("\nSummary Statistics:")
print(summary_stats)

################################################################################
# FINAL SUMMARY
################################################################################

message("\n", strrep("=", 70))
message("ANALYSIS COMPLETE!")
message(strrep("=", 70))

message("\nOutput files created in: ", outdir)
message("\nPlots (PNG with cairo):")
message("  - ", focal_gene_name, "_", focal_id, "_plot_combined_0fold.png")
message("  - ", focal_gene_name, "_", focal_id, "_plot_combined_4fold.png")

message("\nData tables:")
message("  - ", focal_gene_name, "_", focal_id, "_all_fst_statistics_0fold.csv")
message("  - ", focal_gene_name, "_", focal_id, "_all_fst_statistics_4fold.csv")
message("  - ", focal_gene_name, "_", focal_id, "_all_pi_statistics_0fold.csv")
message("  - ", focal_gene_name, "_", focal_id, "_all_pi_statistics_4fold.csv")
message("  - ", focal_gene_name, "_", focal_id, "_all_tajd_statistics_0fold.csv")
message("  - ", focal_gene_name, "_", focal_id, "_all_tajd_statistics_4fold.csv")
message("  - ", focal_gene_name, "_", focal_id, "_summary_statistics_0fold_vs_4fold.csv")

# Check if outlier files were created
if(!is.null(results_0fold$outliers_dip) && sum(results_0fold$outliers_dip, na.rm=TRUE) > 0) {
  message("  - ", focal_gene_name, "_", focal_id, "_diploid_fst_outliers_0fold.csv")
}
if(!is.null(results_0fold$outliers_tet) && sum(results_0fold$outliers_tet, na.rm=TRUE) > 0) {
  message("  - ", focal_gene_name, "_", focal_id, "_tetraploid_fst_outliers_0fold.csv")
}
if(!is.null(results_4fold$outliers_dip) && sum(results_4fold$outliers_dip, na.rm=TRUE) > 0) {
  message("  - ", focal_gene_name, "_", focal_id, "_diploid_fst_outliers_4fold.csv")
}
if(!is.null(results_4fold$outliers_tet) && sum(results_4fold$outliers_tet, na.rm=TRUE) > 0) {
  message("  - ", focal_gene_name, "_", focal_id, "_tetraploid_fst_outliers_4fold.csv")
}

message("\n", strrep("=", 70))
message("All done! Check the output directory for results.")
message(strrep("=", 70))

