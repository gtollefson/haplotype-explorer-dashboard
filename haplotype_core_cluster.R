#!/usr/bin/env Rscript
# =============================================================================
# Haplotype Core Clustering Visualization Tool
# =============================================================================
# Loads VCF.gz files from a directory, filters by coordinate range, and
# visualizes mutant haplotypes clustered around a core allele under positive
# selection. Clustering prioritizes variants immediately adjacent to the core,
# revealing haplotypes most protected from recombination.
# =============================================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(pheatmap)
  library(RColorBrewer)
})

# =============================================================================
# Configuration and CLI
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt_flags <- c("-h", "--help", "--interactive")
  
  if (any(args %in% c("-h", "--help"))) {
    cat("
Haplotype Core Clustering Tool

Usage: Rscript haplotype_core_cluster.R <vcf_dir> [options]

Arguments:
  vcf_dir               Directory containing VCF.gz files (required)

Options (use as Rscript args or set interactively):
  --chr CHR             Chromosome/contig to filter (e.g., 'Pf3D7_13_v3')
  --start START         Start coordinate (bp)
  --end END             End coordinate (bp)
  --core-pos POS        Core allele position (bp) - variant under selection
  --core-allele 0|1     Allele to cluster around: 0=ref, 1=alt (default: 1)
  --pattern PATTERN     Glob pattern for VCF files (default: '*.vcf.gz')
  --metadata FILE       CSV/TSV with sample annotations (sample_id, location, date)
  --gvcf                Use gVCF mode: reports depth at wildtype sites (ref blocks)
  --outdir DIR          Output directory (default: ./haplotype_output)
  --prefix PREFIX       Output file prefix
  --interactive         Launch interactive mode to set parameters

Examples:
  Rscript haplotype_core_cluster.R ./vcf_dir --chr Pf3D7_13_v3 --start 1000000 --end 1500000 --core-pos 1200000
  Rscript haplotype_core_cluster.R ./vcf_dir --interactive
")
    quit(save = "no", status = 0)
  }
  
  # First non-option argument is vcf_dir
  pos_args <- args[!startsWith(args, "-")]
  vcf_dir <- pos_args[1]
  if (is.na(vcf_dir) || vcf_dir == "") {
    message("Error: vcf_dir required. Use --interactive or provide directory path.")
    quit(save = "no", status = 1)
  }
  
  result <- list(
    vcf_dir = vcf_dir,
    chr = "Pf3D7_13_v3",
    start = 1715000,
    end = 1735000,
    core_pos = NULL,
    core_allele = 1L,
    pattern = "*.vcf.gz",
    metadata = NULL,
    gvcf = FALSE,
    outdir = "haplotype_output",
    prefix = "haplotype",
    interactive = "--interactive" %in% args
  )
  
  i <- 2
  while (i <= length(args)) {
    switch(args[i],
      "--chr" = { result$chr <- args[i + 1]; i <- i + 2 },
      "--start" = { result$start <- as.numeric(args[i + 1]); i <- i + 2 },
      "--end" = { result$end <- as.numeric(args[i + 1]); i <- i + 2 },
      "--core-pos" = { result$core_pos <- as.numeric(args[i + 1]); i <- i + 2 },
      "--core-allele" = { result$core_allele <- as.integer(args[i + 1]); i <- i + 2 },
      "--pattern" = { result$pattern <- args[i + 1]; i <- i + 2 },
      "--metadata" = { result$metadata <- args[i + 1]; i <- i + 2 },
      "--gvcf" = { result$gvcf <- TRUE; i <- i + 1 },
      "--outdir" = { result$outdir <- args[i + 1]; i <- i + 2 },
      "--prefix" = { result$prefix <- args[i + 1]; i <- i + 2 },
      "--interactive" = { result$interactive <- TRUE; i <- i + 1 },
      { i <- i + 1 }
    )
  }
  
  result$interactive <- result$interactive || "--interactive" %in% args
  result$gvcf <- result$gvcf || "--gvcf" %in% args
  
  result
}

# =============================================================================
# VCF Loading and Filtering
# =============================================================================

#' Quick summary of VCF: chromosomes and position range (for choosing zoom region)
vcf_summary <- function(vcf_dir, pattern = "*.vcf.gz") {
  vcf <- load_vcf_directory(vcf_dir, pattern)
  fix <- as.data.frame(vcf@fix)
  fix$POS <- as.numeric(fix$POS)
  list(
    chromosomes = unique(fix$CHROM),
    n_variants = nrow(fix),
    position_range = range(fix$POS),
    n_samples = ncol(vcf@gt) - 1,
    sample_names = colnames(vcf@gt)[-1]
  )
}

#' Check if VCF appears to be gVCF (has reference blocks with END in INFO)
is_gvcf <- function(vcf) {
  if (nrow(vcf@fix) == 0) return(FALSE)
  info <- vcf@fix[, "INFO", drop = TRUE]
  any(grepl("END=", info, fixed = TRUE))
}

#' Get genotype and DP at a set of positions from a gVCF (handles reference blocks)
#' vcf: single-sample gVCF; query_pos: numeric positions on same chrom
get_gt_dp_at_positions <- function(vcf, chrom, query_pos) {
  fix <- vcf@fix
  if (nrow(fix) == 0) return(list(gt = character(0), dp = integer(0)))
  pos <- as.numeric(fix[, "POS"])
  info <- fix[, "INFO"]
  end <- rep(NA_integer_, nrow(fix))
  has_end <- grepl("END=", info)
  if (any(has_end)) {
    end[has_end] <- as.integer(sub(".*END=([0-9]+).*", "\\1", info[has_end]))
  }
  end[!has_end] <- pos[!has_end]
  gt_raw <- vcf@gt[, 2]
  dp_raw <- tryCatch(vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE),
                     error = function(e) NULL)
  if (is.null(dp_raw) || all(is.na(dp_raw)))
    dp_raw <- tryCatch(vcfR::extract.gt(vcf, element = "MIN_DP", as.numeric = TRUE),
                       error = function(e) NULL)
  if (is.null(dp_raw)) dp_raw <- rep(NA_integer_, nrow(fix))
  if (is.matrix(dp_raw)) dp_raw <- dp_raw[, 1]
  gt_out <- character(length(query_pos))
  gt_out[] <- "./."
  dp_out <- rep(NA_integer_, length(query_pos))
  for (i in seq_along(query_pos)) {
    p <- query_pos[i]
    in_block <- (fix[, "CHROM"] == chrom) & (pos <= p) & (p <= end)
    if (any(in_block)) {
      idx <- which(in_block)[1]
      gt_out[i] <- gt_raw[idx]
      if (length(dp_raw) >= idx) dp_out[i] <- as.integer(dp_raw[idx])
    }
  }
  list(gt = gt_out, dp = dp_out)
}

#' Load gVCF.gz files from directory and merge
#' Expands reference blocks so wildtype sites have 0/0 + DP (not missing)
load_gvcf_directory <- function(vcf_dir, chr, start, end, pattern = "*.vcf.gz") {
  vcf_files <- list.files(vcf_dir, pattern = glob2rx(pattern), full.names = TRUE)
  vcf_files <- vcf_files[grepl("\\.vcf\\.gz$", vcf_files)]
  if (length(vcf_files) == 0)
    stop("No VCF.gz files found in ", vcf_dir)
  message("Loading ", length(vcf_files), " gVCF file(s)...")
  vcf_list <- lapply(vcf_files, function(f) {
    message("  ", basename(f))
    read.vcfR(f, verbose = FALSE)
  })
  # Variant positions = segregating sites (where ALT is real) + positions in ref blocks
  # that we need for lookup. We only keep actual variant sites to keep matrix small.
  all_var_pos <- integer(0)
  for (v in vcf_list) {
    fix <- v@fix
    if (nrow(fix) == 0) next
    pos <- as.numeric(fix[, "POS"])
    info <- fix[, "INFO"]
    alt <- fix[, "ALT"]
    end_vec <- pos
    has_end <- grepl("END=", info)
    if (any(has_end))
      end_vec[has_end] <- as.integer(sub(".*END=([0-9]+).*", "\\1", info[has_end]))
    in_chr <- fix[, "CHROM"] == chr
    in_range <- (pos <= end) & (end_vec >= start)
    keep <- in_chr & in_range
    if (!any(keep)) next
    pos_keep <- pos[keep]
    alt_keep <- alt[keep]
    is_var <- !grepl("^<NON_REF>$", alt_keep) & alt_keep != "." & nchar(alt_keep) > 0
    if (any(is_var)) all_var_pos <- c(all_var_pos, pos_keep[is_var])
  }
  all_var_pos <- unique(sort(all_var_pos))
  all_var_pos <- all_var_pos[all_var_pos >= start & all_var_pos <= end]
  if (length(all_var_pos) == 0)
    stop("No variant positions in range. Check chr/start/end.")
  n_var <- length(all_var_pos)
  fix_merged <- matrix(".", nrow = n_var, ncol = 8,
                      dimnames = list(NULL, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")))
  fix_merged[, "CHROM"] <- chr
  fix_merged[, "POS"] <- as.character(all_var_pos)
  gt_cols <- character(0)
  dp_list <- list()
  for (v in vcf_list) {
    samp_name <- colnames(v@gt)[2]
    res <- get_gt_dp_at_positions(v, chr, all_var_pos)
    gt_cols <- c(gt_cols, res$gt)
    dp_list[[samp_name]] <- res$dp
  }
  gt_merged <- matrix(c(rep("GT:DP", n_var), gt_cols), nrow = n_var, ncol = 1 + length(vcf_list))
  colnames(gt_merged) <- c("FORMAT", vapply(vcf_list, function(v) colnames(v@gt)[2], character(1)))
  ref_vcf <- vcf_list[[1]]
  ref_vcf@fix <- fix_merged
  ref_vcf@gt <- gt_merged
  attr(ref_vcf, "dp_matrix") <- do.call(rbind, dp_list)
  attr(ref_vcf, "variant_positions") <- all_var_pos
  ref_vcf
}

#' Load VCF.gz files from directory and merge
#' Handles VCFs with different variant sets by taking union of positions
#' Use load_gvcf_directory() for gVCFs to get depth at wildtype sites
#' vcf_files: optional vector of file paths (e.g. from fileInput$datapath); overrides vcf_dir
load_vcf_directory <- function(vcf_dir, pattern = "*.vcf.gz", gvcf = FALSE, chr = NULL, start = NULL, end = NULL, vcf_files = NULL) {
  if (gvcf && !is.null(chr) && !is.null(start) && !is.null(end)) {
    return(load_gvcf_directory(vcf_dir, chr, start, end, pattern))
  }
  if (!is.null(vcf_files) && length(vcf_files) > 0) {
    # Use paths as-is (e.g. from fileInput$datapath; temp paths may lack extension)
  } else {
    vcf_files <- list.files(vcf_dir, pattern = glob2rx(pattern), full.names = TRUE)
    vcf_files <- vcf_files[grepl("\\.vcf\\.gz$", vcf_files)]
  }
  if (length(vcf_files) == 0) {
    stop("No VCF.gz files found.")
  }
  
  message("Loading ", length(vcf_files), " VCF file(s)...")
  vcf_list <- lapply(vcf_files, function(f) {
    message("  ", basename(f))
    read.vcfR(f, verbose = FALSE)
  })
  
  if (length(vcf_list) == 1) {
    return(vcf_list[[1]])
  }
  
  # Union of (CHROM, POS) across all files (samples may have different variant sets)
  all_keys <- unique(do.call(rbind, lapply(vcf_list, function(v) {
    data.frame(CHROM = v@fix[, "CHROM"], POS = as.numeric(v@fix[, "POS"]),
               stringsAsFactors = FALSE)
  })))
  all_keys <- all_keys[order(all_keys$CHROM, all_keys$POS), , drop = FALSE]
  n_var <- nrow(all_keys)
  
  # Build fix: get REF/ALT/etc from first VCF that has each variant
  fix_ref <- vcf_list[[1]]@fix
  fix_merged <- matrix(".", nrow = n_var, ncol = ncol(fix_ref))
  dimnames(fix_merged) <- list(NULL, colnames(fix_ref))
  fix_merged[, "CHROM"] <- all_keys$CHROM
  fix_merged[, "POS"] <- as.character(all_keys$POS)
  u_keys <- paste(all_keys$CHROM, all_keys$POS)
  
  for (v in vcf_list) {
    v_keys <- paste(v@fix[, "CHROM"], v@fix[, "POS"])
    idx <- match(u_keys, v_keys)
    have <- !is.na(idx)
    for (j in colnames(fix_ref)) {
      if (j %in% colnames(v@fix)) {
        empty <- fix_merged[, j] == "." | is.na(fix_merged[, j])
        fill <- have & empty
        if (any(fill)) fix_merged[fill, j] <- v@fix[idx[fill], j]
      }
    }
  }
  
  # Build GT: one column per sample, ././ for missing
  gt_cols <- character(0)
  for (v in vcf_list) {
    v_keys <- paste(v@fix[, "CHROM"], v@fix[, "POS"])
    idx <- match(u_keys, v_keys)
    samp_gt <- rep("./.", n_var)
    have <- !is.na(idx)
    samp_gt[have] <- v@gt[idx[have], 2]
    gt_cols <- c(gt_cols, samp_gt)
  }
  
  format_col <- rep(vcf_list[[1]]@gt[1, 1], n_var)  # FORMAT string same for all
  gt_merged <- matrix(c(format_col, gt_cols), nrow = n_var,
                     ncol = 1 + length(vcf_list))
  samp_names <- vapply(vcf_list, function(v) colnames(v@gt)[2], character(1))
  colnames(gt_merged) <- c("FORMAT", samp_names)
  
  ref_vcf <- vcf_list[[1]]
  ref_vcf@fix <- fix_merged
  ref_vcf@gt <- gt_merged
  ref_vcf
}

#' Filter VCF by chromosome and coordinate range
filter_vcf_by_range <- function(vcf, chr, start, end) {
  fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  fix$POS <- as.numeric(fix$POS)
  
  if (!is.null(chr)) fix <- fix[fix$CHROM == chr, , drop = FALSE]
  if (!is.null(start)) fix <- fix[fix$POS >= start, , drop = FALSE]
  if (!is.null(end)) fix <- fix[fix$POS <= end, , drop = FALSE]
  
  if (nrow(fix) == 0) {
    stop("No variants in specified range. Check chr/start/end.")
  }
  
  idx <- match(paste(fix$CHROM, fix$POS), paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"]))
  vcf@fix <- vcf@fix[idx, , drop = FALSE]
  vcf@gt <- vcf@gt[idx, , drop = FALSE]
  vcf
}

#' Load metadata file (CSV or TSV) for sample annotations.
#' Expects: sample_id (or first column) matching VCF sample names,
#' plus location and date columns (case-insensitive: location, Location,
#' collection_location, etc.; date, Date, collection_date, year, etc.)
load_metadata <- function(metadata_path, sample_ids) {
  if (is.null(metadata_path) || !file.exists(metadata_path)) return(NULL)
  ext <- tolower(tools::file_ext(metadata_path))
  meta <- if (ext == "tsv") {
    read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  nc <- ncol(meta)
  if (nc < 2) return(NULL)
  nms <- tolower(names(meta))
  # Sample ID: prefer seq_sample_id (matches VCF), then sample_id, sample, id
  sid_col <- which(nms %in% c("seq_sample_id", "sample_id", "sample", "id", "samples"))[1]
  if (is.na(sid_col)) sid_col <- 1
  # Location: exact or partial match (origin, location, site, country, region)
  loc_col <- which(nms %in% c("location", "collection_location", "site", "country", "region"))[1]
  if (is.na(loc_col))
    loc_col <- which(grepl("origin|location|site|country|region", nms))[1]
  # Date: prefer date of sample storage, then collection/sequence dates
  date_col <- which(grepl("storage", nms) & grepl("date", nms))[1]
  if (is.na(date_col))
    date_col <- which(nms %in% c("date_of_storage", "date.of.storage", "storage_date",
                                 "date", "collection_date", "year", "collection_year",
                                 "sequence_date"))[1]
  if (is.na(date_col))
    date_col <- which(grepl("date|year", nms))[1]
  meta <- meta[meta[[sid_col]] %in% sample_ids, , drop = FALSE]
  if (nrow(meta) == 0) return(NULL)
  rownames(meta) <- meta[[sid_col]]
  list(
    df = meta,
    sample_col = names(meta)[sid_col],
    location_col = if (!is.na(loc_col)) names(meta)[loc_col] else NULL,
    date_col = if (!is.na(date_col)) names(meta)[date_col] else NULL
  )
}

#' Extract genotype matrix: 0=ref, 1=het, 2=alt, NA=missing
#' Attaches: dp_matrix (from gVCF), af_matrix (AF from FORMAT or GT-derived)
extract_genotype_matrix <- function(vcf) {
  gt <- vcfR::extract.gt(vcf, element = "GT", as.numeric = FALSE)
  
  # Convert GT to numeric (0/0->0, 0/1|1/0->1, 1/1->2)
  gt_num <- apply(gt, c(1, 2), function(x) {
    if (is.na(x) || x == "." || x == "./.") return(NA)
    alleles <- strsplit(x, "[/|]")[[1]]
    sum(as.integer(alleles), na.rm = TRUE)
  })
  
  dimnames(gt_num) <- dimnames(gt)
  out <- t(gt_num)  # samples x variants
  
  # DP: from gVCF attr or extract from FORMAT (GT:DP, etc.)
  dp_mat <- attr(vcf, "dp_matrix")
  if (is.null(dp_mat)) {
    dp_mat <- tryCatch({
      raw <- vcfR::extract.gt(vcf, element = "DP", as.numeric = TRUE)
      if (is.null(raw)) stop("no DP")
      t(raw)
    }, error = function(e) {
      tryCatch({
        raw <- vcfR::extract.gt(vcf, element = "MIN_DP", as.numeric = TRUE)
        if (is.null(raw)) NULL else t(raw)
      }, error = function(e2) NULL)
    })
  }
  if (!is.null(dp_mat))
    attr(out, "dp_matrix") <- dp_mat
  
  # AF: try FORMAT AF, else AD (ref,alt), else derive from GT
  af_mat <- tryCatch({
    raw <- vcfR::extract.gt(vcf, element = "AF", as.numeric = FALSE)
    if (is.null(raw)) stop("no AF")
    m <- apply(raw, c(1, 2), function(x) {
      if (is.na(x) || x == ".") return(NA)
      v <- as.numeric(strsplit(x, ",")[[1]][1])
      if (is.na(v)) NA else v
    })
    t(m)
  }, error = function(e) NULL)
  if (is.null(af_mat)) {
    af_mat <- tryCatch({
      raw <- vcfR::extract.gt(vcf, element = "AD", as.numeric = FALSE)
      if (is.null(raw)) stop("no AD")
      m <- apply(raw, c(1, 2), function(x) {
        if (is.na(x) || x == ".") return(NA)
        parts <- strsplit(x, ",")[[1]]
        ref <- as.numeric(parts[1])
        alt <- if (length(parts) > 1) as.numeric(parts[2]) else 0
        tot <- ref + alt
        if (tot == 0) NA else alt / tot
      })
      t(m)
    }, error = function(e) NULL)
  }
  if (is.null(af_mat))
    af_mat <- out / 2  # 0->0, 1->0.5, 2->1
  attr(out, "af_matrix") <- af_mat
  out
}

# =============================================================================
# Core-Allele Clustering
# =============================================================================

#' Find core variant index by position
find_core_index <- function(vcf, core_pos) {
  pos <- as.numeric(vcf@fix[, "POS"])
  idx <- which.min(abs(pos - core_pos))
  if (length(idx) == 0) stop("Core position ", core_pos, " not in variant set")
  idx
}

#' Proximity-weighted distance: variants closer to core have higher weight.
#' Adjacent variants (immediate neighbors) dominate clustering - haplotypes
#' most protected from recombination share the core + adjacent alleles.
#' Returns distance matrix for samples
proximity_weighted_distance <- function(gt_matrix, variant_positions, core_idx, 
                                        core_allele = 1L) {
  n_var <- ncol(gt_matrix)
  dist_to_core <- abs(variant_positions - variant_positions[core_idx])
  dist_to_core[core_idx] <- 0  # core has max weight
  
  # Weight: inverse distance - adjacent variants (closest) get highest weight
  # Core and immediate neighbors dominate; distant variants contribute less
  nonzero <- dist_to_core[dist_to_core > 0]
  bp_per_unit <- if (length(nonzero) > 0) max(1, min(nonzero) / 10) else 1000
  weights <- 1 / (dist_to_core + bp_per_unit)
  maxw <- max(weights[is.finite(weights)], na.rm = TRUE)
  weights[core_idx] <- maxw * 3  # core gets highest weight
  
  # Carrier status: add large penalty so carrier/non-carrier cluster separately
  core_gt <- gt_matrix[, core_idx]
  is_carrier <- if (core_allele == 1) {
    core_gt >= 1 & !is.na(core_gt)
  } else {
    core_gt <= 1 & !is.na(core_gt)
  }
  
  # Weighted Hamming (vectorized for speed)
  n_samp <- nrow(gt_matrix)
  D <- matrix(0, n_samp, n_samp)
  
  for (j in seq_len(n_var)) {
    g <- gt_matrix[, j]
    g[is.na(g)] <- -9  # treat NA as different from everything for clustering
    diff_mat <- outer(g, g, "!=")
    D <- D + diff_mat * weights[j]
  }
  
  # Penalty: carrier vs non-carrier get large distance (ensures dendrogram groups them)
  carrier_penalty <- max(D) * 10
  carrier_diff <- outer(is_carrier, is_carrier, "!=")
  D <- D + carrier_penalty * carrier_diff
  
  as.dist(D)
}

#' Hierarchical clustering with core-aware ordering.
#' Distance matrix enforces: 1) carriers vs non-carriers split first,
#' 2) within carriers, prioritizes similarity at adjacent variants (near core)
cluster_haplotypes <- function(gt_matrix, variant_positions, core_idx, core_allele = 1L) {
  core_gt <- gt_matrix[, core_idx]
  is_carrier <- if (core_allele == 1) {
    core_gt >= 1 & !is.na(core_gt)
  } else {
    core_gt <= 1 & !is.na(core_gt)
  }
  
  d <- proximity_weighted_distance(gt_matrix, variant_positions, core_idx, core_allele)
  hc <- hclust(d, method = "ward.D2")
  
  # Dendrogram order: carrier penalty ensures carriers cluster together
  # Within carrier/non-carrier groups, order reflects adjacency similarity
  list(
    hclust = hc,
    is_carrier = is_carrier,
    distance = d,
    sample_order = hc$order
  )
}

# =============================================================================
# Visualization
# =============================================================================

#' Main haplotype heatmap with core-allele clustering
#' dp_matrix: optional samples x variants depth matrix (from gVCF); if provided, shows depth as cell labels
plot_haplotype_heatmap <- function(gt_matrix, variant_positions, core_idx,
                                   cluster_result, vcf, chr_label, core_allele = 1L,
                                   metadata = NULL, dp_matrix = NULL, outpath = NULL) {
  
  # Order samples by dendrogram
  ord <- cluster_result$sample_order
  sample_ids <- rownames(gt_matrix)[ord]
  mat_ordered <- gt_matrix[ord, , drop = FALSE]
  
  # Annotation: core carrier status + metadata
  annot_row <- data.frame(
    Core_Carrier = ifelse(cluster_result$is_carrier[ord], "Carrier", "Non-carrier"),
    row.names = rownames(mat_ordered)
  )
  if (!is.null(metadata) && !is.null(metadata$location_col) && metadata$location_col %in% names(metadata$df)) {
    loc_vals <- rep("—", length(sample_ids))
    in_meta <- sample_ids %in% rownames(metadata$df)
    loc_vals[in_meta] <- metadata$df[sample_ids[in_meta], metadata$location_col]
    loc_vals[is.na(loc_vals)] <- "—"
    annot_row$Location <- loc_vals
  }
  if (!is.null(metadata) && !is.null(metadata$date_col) && metadata$date_col %in% names(metadata$df)) {
    date_vals <- rep("—", length(sample_ids))
    in_meta <- sample_ids %in% rownames(metadata$df)
    date_vals[in_meta] <- as.character(metadata$df[sample_ids[in_meta], metadata$date_col])
    date_vals[is.na(date_vals)] <- "—"
    annot_row$Date <- date_vals
  }
  
  # Row labels: sample + location + date for quick identification
  lab_row <- rownames(mat_ordered)
  if (!is.null(metadata)) {
    loc <- rep("", length(sample_ids))
    dat <- rep("", length(sample_ids))
    in_meta <- sample_ids %in% rownames(metadata$df)
    if (!is.null(metadata$location_col) && metadata$location_col %in% names(metadata$df))
      loc[in_meta] <- metadata$df[sample_ids[in_meta], metadata$location_col]
    if (!is.null(metadata$date_col) && metadata$date_col %in% names(metadata$df))
      dat[in_meta] <- as.character(metadata$df[sample_ids[in_meta], metadata$date_col])
    lab_row <- paste0(rownames(mat_ordered), " | ", loc, " | ", dat)
    lab_row[loc == "" & dat == ""] <- rownames(mat_ordered)[loc == "" & dat == ""]
  }
  # Ensure unique row names (pheatmap requires unique rownames for annotation)
  if (anyDuplicated(lab_row)) {
    lab_row <- make.unique(as.character(lab_row), sep = "_")
  }
  rownames(mat_ordered) <- lab_row
  rownames(annot_row) <- lab_row  # keep annotation in sync
  
  # Column labels: full position coordinates
  colnames(mat_ordered) <- format(variant_positions, scientific = FALSE, trim = TRUE)
  
  # Validate dimensions
  stopifnot(nrow(mat_ordered) == nrow(annot_row))
  
  # Gap between carrier and non-carrier groups (if both present)
  carrier_in_order <- cluster_result$is_carrier[ord]
  n_carriers <- sum(carrier_in_order)
  gaps_row <- if (n_carriers > 0 && n_carriers < length(carrier_in_order)) {
    n_carriers  # gap after last carrier
  } else {
    NULL
  }
  
  # Column annotation: mark core allele column in red
  cnames <- colnames(mat_ordered)
  annot_col <- data.frame(
    Core = factor(rep("Other", ncol(mat_ordered)), levels = c("Core", "Other")),
    row.names = cnames
  )
  annot_col[core_idx, "Core"] <- "Core"
  
  # Annotation colors
  annot_colors <- list(
    Core_Carrier = c(Carrier = "#D95F02", "Non-carrier" = "#7570B3"),
    Core = c(Core = "#E63946", Other = "white")  # red for core allele column
  )
  if ("Location" %in% names(annot_row)) {
    loc_vals <- sort(unique(annot_row$Location))
    n_loc <- length(loc_vals)
    if (n_loc > 0 && n_loc <= 12) {
      pals <- brewer.pal(max(3, min(n_loc, 12)), "Set3")
      annot_colors$Location <- setNames(pals[seq_len(n_loc)], loc_vals)
    }
  }
  
  if (!is.null(outpath)) {
    pdf(outpath, width = 16, height = max(8, nrow(mat_ordered) * 0.2))
  }
  
  # Optional: show depth (DP) as cell labels when from gVCF
  pheatmap_args <- list(
    mat = mat_ordered,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = annot_row,
    annotation_col = annot_col,
    annotation_colors = annot_colors,
    color = c("#A0A0A0", "#FFFFFF", "#89CFF0", "#1A5276"),  # missing=grey, ref=white, het=light blue, hom=dark blue
    breaks = c(-0.5, 0.25, 0.75, 1.25, 2.5),
    legend_breaks = c(0, 1, 2),
    legend_labels = c("Ref", "Het", "Hom"),
    na_col = "#A0A0A0",
    main = paste0("Haplotype matrix | Chr ", chr_label, " | Core at ", 
                  format(variant_positions[core_idx], scientific = FALSE), " (selection)"),
    fontsize_row = 6,
    fontsize_col = 8,
    border_color = NA,
    gaps_row = gaps_row,
    show_rownames = TRUE
  )
  if (!is.null(dp_matrix)) {
    dp_ord <- dp_matrix[ord, , drop = FALSE]
    disp_mat <- matrix(as.character(dp_ord), nrow = nrow(dp_ord))
    disp_mat[is.na(dp_ord) | dp_ord == ""] <- ""
    pheatmap_args$display_numbers <- disp_mat
    pheatmap_args$number_format <- "%s"
    pheatmap_args$fontsize_number <- 5
  }
  do.call(pheatmap, pheatmap_args)
  
  if (!is.null(outpath)) {
    dev.off()
    message("Saved: ", outpath)
  }
}

#' Dendrogram of haplotype clustering (carriers cluster together; within groups,
#' similarity at adjacent-to-core variants drives sub-clustering)
plot_haplotype_dendrogram <- function(cluster_result, outpath = NULL) {
  hc <- cluster_result$hclust
  
  if (!is.null(outpath)) {
    pdf(outpath, width = 12, height = max(6, length(cluster_result$is_carrier) * 0.08))
  }
  
  par(mar = c(5, 4, 4, 2))
  plot(hc, hang = -1, cex = 0.5, 
       main = "Haplotype dendrogram (proximity-weighted: adjacent variants prioritized)")
  legend("topright", legend = c("Core carrier", "Non-carrier"), 
         fill = c("#D95F02", "#7570B3"), bty = "n", cex = 0.9)
  
  if (!is.null(outpath)) {
    dev.off()
    message("Saved: ", outpath)
  }
}

#' Summary of core haplotype clusters
summarize_clusters <- function(gt_matrix, cluster_result, variant_positions, core_idx) {
  carriers <- which(cluster_result$is_carrier)
  if (length(carriers) == 0) return(NULL)
  
  carrier_gt <- gt_matrix[carriers, , drop = FALSE]
  # Cluster carriers by adjacent-variant similarity
  d_carrier <- as.dist(as.matrix(cluster_result$distance)[carriers, carriers])
  hc_c <- hclust(d_carrier, "ward.D2")
  
  n_clust <- min(10, max(2, length(carriers) %/% 5))
  clust_assign <- cutree(hc_c, k = n_clust)
  
  tibble(
    Sample = rownames(carrier_gt)[carriers],
    Core_Cluster = clust_assign,
    Core_Genotype = carrier_gt[, core_idx]
  )
}

# =============================================================================
# Interactive Mode
# =============================================================================

run_interactive <- function(vcf_dir) {
  message("\n=== Interactive Haplotype Core Clustering ===\n")
  
  vcf <- load_vcf_directory(vcf_dir)
  
  # Show available chromosomes and position range
  fix <- as.data.frame(vcf@fix)
  fix$POS <- as.numeric(fix$POS)
  chrs <- unique(fix$CHROM)
  pos_range <- range(fix$POS)
  
  message("Available chromosomes: ", paste(head(chrs, 5), collapse = ", "), 
          if (length(chrs) > 5) " ..." else "")
  message("Position range: ", min(pos_range), " - ", max(pos_range))
  
  chr <- readline(prompt = "Enter chromosome: ")
  if (chr == "") chr <- chrs[1]
  
  start <- as.numeric(readline(prompt = "Enter start position: "))
  if (is.na(start)) start <- max(1, pos_range[1] - 50000)
  
  end <- as.numeric(readline(prompt = "Enter end position: "))
  if (is.na(end)) end <- min(pos_range[2] + 50000, Inf)
  
  vcf_filtered <- filter_vcf_by_range(vcf, chr, start, end)
  positions <- as.numeric(vcf_filtered@fix[, "POS"])
  
  message("\nVariants in range: ", length(positions))
  message("Position range: ", min(positions), " - ", max(positions))
  
  core_pos <- as.numeric(readline(prompt = "Enter core allele position (bp): "))
  if (is.na(core_pos)) core_pos <- median(positions)
  
  core_allele <- readline(prompt = "Core allele: 0=ref, 1=alt [1]: ")
  core_allele <- if (core_allele == "0") 0L else 1L
  
  meta_path <- readline(prompt = "Metadata file path (CSV/TSV with sample_id, location, date) [skip]: ")
  meta_path <- if (nzchar(trimws(meta_path))) trimws(meta_path) else NULL
  
  run_analysis(vcf_filtered, positions, core_pos, core_allele, chr, 
               outdir = "haplotype_output", prefix = "haplotype", metadata = meta_path)
}

# =============================================================================
# Main Analysis Pipeline
# =============================================================================

run_analysis <- function(vcf_filtered, variant_positions, core_pos, core_allele,
                        chr_label, outdir, prefix, metadata = NULL) {
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  gt_matrix <- extract_genotype_matrix(vcf_filtered)
  core_idx <- find_core_index(vcf_filtered, core_pos)
  
  message("\nSamples: ", nrow(gt_matrix), " | Variants: ", ncol(gt_matrix))
  message("Core allele at position ", variant_positions[core_idx], 
          " (index ", core_idx, ")")
  
  n_carriers <- sum(gt_matrix[, core_idx] >= 1, na.rm = TRUE)
  message("Core allele carriers: ", n_carriers, "/", nrow(gt_matrix))
  
  cluster_result <- cluster_haplotypes(gt_matrix, variant_positions, core_idx, core_allele)
  
  # Load metadata for annotations if path provided
  meta_obj <- if (is.character(metadata)) {
    load_metadata(metadata, rownames(gt_matrix))
  } else if (is.list(metadata)) metadata else NULL
  
  # Outputs
  heatmap_path <- file.path(outdir, paste0(prefix, "_heatmap.pdf"))
  dendro_path <- file.path(outdir, paste0(prefix, "_dendrogram.pdf"))
  cluster_path <- file.path(outdir, paste0(prefix, "_cluster_summary.csv"))
  
  dp_mat <- attr(gt_matrix, "dp_matrix")
  plot_haplotype_heatmap(gt_matrix, variant_positions, core_idx, cluster_result,
                         vcf_filtered, chr_label, core_allele, meta_obj, dp_mat, heatmap_path)
  
  plot_haplotype_dendrogram(cluster_result, dendro_path)
  
  summary_df <- summarize_clusters(gt_matrix, cluster_result, variant_positions, core_idx)
  if (!is.null(summary_df)) {
    write.csv(summary_df, cluster_path, row.names = FALSE)
    message("Saved: ", cluster_path)
  }
  
  # Return objects for further use
  invisible(list(
    gt_matrix = gt_matrix,
    variant_positions = variant_positions,
    core_idx = core_idx,
    cluster_result = cluster_result,
    summary = summary_df
  ))
}

# =============================================================================
# Entry Point
# =============================================================================

main <- function() {
  args <- parse_args()
  
  if (!dir.exists(args$vcf_dir)) {
    stop("VCF directory not found: ", args$vcf_dir)
  }
  
  if (args$interactive) {
    run_interactive(args$vcf_dir)
    return(invisible())
  }
  
  
  if (args$gvcf) {
    message("gVCF mode: loading with reference block expansion (depth at wildtype sites)")
    vcf_filtered <- load_vcf_directory(args$vcf_dir, args$pattern,
                                       gvcf = TRUE, chr = args$chr,
                                       start = args$start, end = args$end)
    variant_positions <- attr(vcf_filtered, "variant_positions")
    if (is.null(variant_positions)) variant_positions <- as.numeric(vcf_filtered@fix[, "POS"])
  } else {
    vcf <- load_vcf_directory(args$vcf_dir, args$pattern)
    vcf_filtered <- filter_vcf_by_range(vcf, args$chr, args$start, args$end)
    variant_positions <- as.numeric(vcf_filtered@fix[, "POS"])
  }
  
  # Auto core position if not specified
  core_pos <- args$core_pos
  if (is.null(core_pos)) {
    core_pos <- median(variant_positions)
    message("Core position auto-set to median: ", core_pos)
  }
  
  run_analysis(vcf_filtered, variant_positions, core_pos, args$core_allele,
              args$chr, args$outdir, args$prefix, metadata = args$metadata)
}

# Run when executed as script (not when sourced)
if (!interactive()) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  is_direct <- length(file_arg) > 0 && grepl("haplotype_core_cluster", file_arg, ignore.case = TRUE)
  if (is_direct) main()
}
