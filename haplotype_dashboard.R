#!/usr/bin/env Rscript
# =============================================================================
# Interactive Haplotype Core Clustering Dashboard
# =============================================================================
# Launches a Shiny app in the browser for:
#   - Loading VCF + metadata
#   - Interactive heatmap with location/date labels
#   - Drag-and-drop sample reordering
#   - Re-apply unsupervised clustering on demand
#
# Run: Rscript haplotype_dashboard.R
# Or:  shiny::runApp("haplotype_dashboard.R")
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(plotly)
  library(sortable)
  library(vcfR)
  library(dplyr)
})

# Source core functions from haplotype_core_cluster.R (same directory as this script)
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("--file=", args, fixed = TRUE)
  if (length(m)) dirname(normalizePath(sub("--file=", "", args[m]))) else "."
}, error = function(e) ".")
for (d in c(script_dir, ".", getwd(), dirname(getwd()))) {
  p <- file.path(d, "haplotype_core_cluster.R")
  if (file.exists(p)) { source(p); break }
}

# =============================================================================
# UI
# =============================================================================

ui <- fluidPage(
  titlePanel("Haplotype Core Clustering Dashboard"),
  # Top: Data & Region controls
  wellPanel(
    fluidRow(
      column(4,
        h4("Data"),
        fileInput("vcf_files", "VCF files (.vcf.gz)", multiple = TRUE,
                  accept = ""),
        helpText("Select .vcf.gz (and optionally .vcf). .csi and .tbi index files are ignored."),
        fileInput("metadata_file", "Metadata (CSV/TSV)", accept = c(".csv", ".tsv", ".txt"))
      ),
      column(4,
        h4("Region"),
        textInput("chr", "Chromosome", value = "Pf3D7_13_v3"),
        fluidRow(
          column(6, numericInput("start", "Start (bp)", value = 1715000)),
          column(6, numericInput("end", "End (bp)", value = 1735000))
        ),
        fluidRow(
          column(6, numericInput("core_pos", "Core position (bp)", value = NA)),
          column(6, selectInput("core_allele", "Core allele", choices = c("alt (1)" = 1, "ref (0)" = 0)))
        ),
        checkboxInput("filter_variant_sites", "Only variant sites (exclude sites with all ref)",
                     value = TRUE)
      ),
      column(4,
        h4("Load & cluster"),
        actionButton("load_btn", "Load data", class = "btn-primary btn-block"),
        br(),
        actionButton("cluster_btn", "Apply clustering", class = "btn-success btn-block")
      )
    )
  ),
  # Below: Sample order (left) | Heatmap (right)
  fluidRow(
    column(3,
      wellPanel(
        h4("Sample order"),
        p("Drag to reorder. Heatmap updates live."),
        uiOutput("sample_rank_list")
      )
    ),
    column(9,
      plotlyOutput("heatmap", height = "800px"),
      verbatimTextOutput("status")
    )
  )
)

# =============================================================================
# Server
# =============================================================================

server <- function(input, output, session) {
  
  # Reactive values
  vdata <- reactiveValues(
    gt_matrix = NULL,
    variant_positions = NULL,
    core_idx = NULL,
    cluster_result = NULL,
    metadata = NULL,
    sample_order = NULL
  )
  
  # Load data
  observeEvent(input$load_btn, {
    req(input$vcf_files, nrow(input$vcf_files) > 0)
    
    # Filter to VCF files only (exclude .csi, .tbi index files)
    df <- input$vcf_files
    vcf_idx <- grepl("\\.vcf\\.gz$|\\.vcf$", df$name, ignore.case = TRUE)
    if (!any(vcf_idx)) {
      showNotification("No VCF files (.vcf.gz or .vcf) found in selection. Index files (.csi, .tbi) are ignored.", type = "error")
      return()
    }
    vcf_paths <- df$datapath[vcf_idx]
    
    withProgress(message = "Loading VCF...", value = 0, {
      tryCatch({
        vcf <- load_vcf_directory(vcf_dir = NULL, vcf_files = vcf_paths)
        incProgress(0.3)
        
        chr <- if (nzchar(trimws(input$chr))) trimws(input$chr) else NULL
        start <- if (is.na(input$start)) NULL else input$start
        end <- if (is.na(input$end)) NULL else input$end
        
        if (is.null(chr) || is.null(start) || is.null(end)) {
          fix <- as.data.frame(vcf@fix)
          fix$POS <- as.numeric(fix$POS)
          chr <- if (is.null(chr) || !nzchar(chr)) unique(fix$CHROM)[1] else chr
          start <- if (is.null(start) || is.na(start)) max(1, min(fix$POS) - 50000) else start
          end <- if (is.null(end) || is.na(end)) max(fix$POS) + 50000 else end
        }
        
        vcf_sub <- filter_vcf_by_range(vcf, chr, start, end)
        incProgress(0.5)
        
        gt_matrix <- extract_genotype_matrix(vcf_sub)
        positions <- as.numeric(vcf_sub@fix[, "POS"])
        core_pos <- if (is.na(input$core_pos)) median(positions) else input$core_pos
        core_idx <- find_core_index(vcf_sub, core_pos)
        core_allele <- as.integer(input$core_allele)
        
        # Optionally filter to sites variant in at least one sample (always keep core)
        if (isTRUE(input$filter_variant_sites)) {
          is_variant <- colSums(gt_matrix >= 1, na.rm = TRUE) > 0
          keep_core <- (seq_len(ncol(gt_matrix)) == core_idx)
          keep_cols <- is_variant | keep_core
          dp_mat <- attr(gt_matrix, "dp_matrix")
          af_mat <- attr(gt_matrix, "af_matrix")
          gt_matrix <- gt_matrix[, keep_cols, drop = FALSE]
          if (!is.null(dp_mat)) attr(gt_matrix, "dp_matrix") <- dp_mat[, keep_cols, drop = FALSE]
          if (!is.null(af_mat)) attr(gt_matrix, "af_matrix") <- af_mat[, keep_cols, drop = FALSE]
          positions <- positions[keep_cols]
          core_idx <- which.min(abs(positions - core_pos))
        }
        
        cluster_result <- cluster_haplotypes(gt_matrix, positions, core_idx, core_allele)
        incProgress(0.8)
        
        # Metadata
        meta <- NULL
        if (!is.null(input$metadata_file) && nrow(input$metadata_file) > 0) {
          path <- input$metadata_file$datapath[1]
          meta <- load_metadata(path, rownames(gt_matrix))
        }
        incProgress(1)
        
        vdata$gt_matrix <- gt_matrix
        vdata$variant_positions <- positions
        vdata$core_idx <- core_idx
        vdata$cluster_result <- cluster_result
        vdata$metadata <- meta
        vdata$sample_order <- cluster_result$sample_order
        vdata$range_start <- start
        vdata$range_end <- end
        
        updateNumericInput(session, "start", value = start)
        updateNumericInput(session, "end", value = end)
        updateNumericInput(session, "core_pos", value = core_pos)
        updateTextInput(session, "chr", value = chr)
        
        showNotification("Data loaded.", type = "message")
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
    })
  })
  
  # Re-apply clustering
  observeEvent(input$cluster_btn, {
    req(vdata$gt_matrix, vdata$cluster_result)
    vdata$sample_order <- vdata$cluster_result$sample_order
    showNotification("Clustering applied.", type = "message")
  })
  
  # Rank list UI - draggable sample order
  output$sample_rank_list <- renderUI({
    req(vdata$gt_matrix)
    sample_ids <- rownames(vdata$gt_matrix)
    ord <- vdata$sample_order
    if (!is.null(ord)) sample_ids <- sample_ids[ord]
    
    # Labels: sample | location | date
    labels <- sample_ids
    if (!is.null(vdata$metadata)) {
      meta <- vdata$metadata$df
      loc <- rep("", length(sample_ids))
      dat <- rep("", length(sample_ids))
      in_meta <- sample_ids %in% rownames(meta)
      if (!is.null(vdata$metadata$location_col))
        loc[in_meta] <- meta[sample_ids[in_meta], vdata$metadata$location_col]
      if (!is.null(vdata$metadata$date_col))
        dat[in_meta] <- as.character(meta[sample_ids[in_meta], vdata$metadata$date_col])
      labels <- paste0(sample_ids, " | ", loc, " | ", dat)
    }
    
    rank_list(
      text = "Drag to reorder",
      labels = labels,
      input_id = "sample_order",
      options = sortable_options(multiDrag = TRUE)
    )
  })
  
  # Sync rank_list order back when it changes (user dragged)
  observeEvent(input$sample_order, {
    req(vdata$gt_matrix, input$sample_order)
    rl <- input$sample_order
    sample_ids <- rownames(vdata$gt_matrix)
    # Extract sample IDs from "sample | loc | date" format
    extracted <- sub(" \\| .*", "", rl)
    if (all(extracted %in% sample_ids) && length(extracted) == length(sample_ids)) {
      new_ord <- match(extracted, sample_ids)
      vdata$sample_order <- new_ord
    }
  }, ignoreInit = TRUE)
  
  # Interactive heatmap
  output$heatmap <- renderPlotly({
    req(vdata$gt_matrix, vdata$variant_positions, vdata$core_idx)
    
    gt <- vdata$gt_matrix
    pos <- vdata$variant_positions
    core_idx <- vdata$core_idx
    range_start <- if (is.null(vdata$range_start)) min(pos) else vdata$range_start
    range_end <- if (is.null(vdata$range_end)) max(pos) else vdata$range_end
    
    # Order (ensure integer indices)
    ord <- vdata$sample_order
    if (is.null(ord)) ord <- seq_len(nrow(gt))
    ord <- as.integer(ord)
    mat <- as.matrix(gt[ord, , drop = FALSE])
    
    # Row labels with metadata
    rlabels <- rownames(mat)
    if (!is.null(vdata$metadata)) {
      meta <- vdata$metadata$df
      loc <- rep("", nrow(mat))
      dat <- rep("", nrow(mat))
      sids <- rownames(mat)
      in_meta <- sids %in% rownames(meta)
      if (!is.null(vdata$metadata$location_col))
        loc[in_meta] <- meta[sids[in_meta], vdata$metadata$location_col]
      if (!is.null(vdata$metadata$date_col))
        dat[in_meta] <- as.character(meta[sids[in_meta], vdata$metadata$date_col])
      rlabels <- paste0(rownames(mat), " | ", loc, " | ", dat)
    }
    
    mat_plot <- as.matrix(mat)
    mat_plot[is.na(mat_plot)] <- -0.5  # missing -> gray
    n_r <- nrow(mat_plot)
    n_c <- ncol(mat_plot)
    if (length(rlabels) != n_r) rlabels <- if (!is.null(rownames(mat_plot))) rownames(mat_plot) else paste0("Sample_", seq_len(n_r))
    
    # Use integer indices for x so all cells have uniform width/height
    x_idx <- seq_along(pos)
    pos_labels <- format(pos, scientific = FALSE)
    dp_mat <- attr(vdata$gt_matrix, "dp_matrix")
    af_mat <- attr(vdata$gt_matrix, "af_matrix")
    if (!is.null(dp_mat)) dp_mat <- dp_mat[ord, , drop = FALSE]
    if (!is.null(af_mat)) af_mat <- af_mat[ord, , drop = FALSE]
    geno_labels <- c("0/0", "0/1", "1/1", ".")
    geno_str <- matrix(geno_labels[match(mat_plot, c(0, 1, 2, -0.5), nomatch = 4)], nrow = n_r, ncol = n_c)
    text_mat <- matrix("", nrow = n_r, ncol = n_c)
    for (i in seq_len(n_r)) {
      for (j in seq_len(n_c)) {
        parts <- c(paste0("Position: ", pos_labels[j]))
        if (!is.null(dp_mat)) parts <- c(parts, paste0("DP: ", if (is.na(dp_mat[i,j])) "—" else round(dp_mat[i,j])))
        if (!is.null(af_mat)) parts <- c(parts, paste0("AF: ", if (is.na(af_mat[i,j])) "—" else format(round(af_mat[i,j], 3), nsmall = 3)))
        parts <- c(parts, paste0("Genotype: ", geno_str[i,j]))
        text_mat[i,j] <- paste(parts, collapse = "<br>")
      }
    }
    n_ticks <- min(25, length(pos))
    tick_idx <- unique(round(seq(1, length(pos), length.out = n_ticks)))
    tick_idx <- tick_idx[tick_idx >= 1 & tick_idx <= length(pos)]
    p <- plot_ly(
      z = mat_plot,
      x = x_idx,
      y = rlabels,
      text = text_mat,
      type = "heatmap",
      colorscale = list(  # missing=grey, ref=white, het=light blue, hom=dark blue
        c(0, "#A0A0A0"), c(0.2, "#FFFFFF"), c(0.6, "#89CFF0"), c(1, "#1A5276")
      ),
      zmin = -0.5, zmax = 2,
      hovertemplate = "%{y}<br>%{text}<extra></extra>"
    ) %>%
      layout(
        title = paste0("Haplotype matrix | Core at ", format(pos[core_idx], scientific = FALSE), " (red)"),
        xaxis = list(
          title = "Position",
          tickvals = tick_idx,
          ticktext = pos_labels[tick_idx],
          range = c(0.5, length(pos) + 0.5)
        ),
        yaxis = list(title = "", autorange = "reversed", tickfont = list(size = 9)),
        margin = list(l = max(200, max(nchar(rlabels)) * 5)),
        shapes = list(
          list(
            type = "rect",
            x0 = core_idx - 0.5, x1 = core_idx + 0.5,
            y0 = -0.5, y1 = n_r - 0.5,
            xref = "x", yref = "y",
            line = list(color = "#E63946", width = 3),
            fillcolor = "rgba(230, 57, 70, 0.08)"
          )
        )
      ) %>%
      config(toImageButtonOptions = list(format = "png", filename = "haplotype_heatmap"))
    
    p
  })
  
  # Status
  output$status <- renderText({
    req(vdata$gt_matrix)
    n <- nrow(vdata$gt_matrix)
    m <- ncol(vdata$gt_matrix)
    nc <- sum(vdata$cluster_result$is_carrier, na.rm = TRUE)
    paste0("Samples: ", n, " | Variants: ", m, " | Core carriers: ", nc)
  })
}

# =============================================================================
# Run app
# =============================================================================

# Support both runApp() and Rscript
if (interactive()) {
  shinyApp(ui, server)
} else {
  # Rscript mode: launch in browser
  runApp(
    list(ui = ui, server = server),
    launch.browser = TRUE
  )
}
