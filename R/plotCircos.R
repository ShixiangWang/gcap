#' Plot Circos for GCAP object
#'
#' @inheritParams gcap.plotProfile
#' @param highlight_genes gene list to highlight.
#' @param clust_distance a distance as cutoff for different clusters.
#' Default is 1e7, i.e. 10Mb. Note 100 Mb is set to genes on different
#' chromosomes, so please don't set value larger than that.
#' @param col length-2 colors for circular and noncircular.
#' @param genome_build genome version.
#' @param chrs chromosome names.
#' @param ideogram_height ideogram height at default.
#'
#' @return Nothing.
#' @export
gcap.plotCircos <- function(fCNA,
                            highlight_genes = NULL,
                            clust_distance = 1e7,
                            col = c("#FF000080", "#0000FF80"),
                            genome_build = c("hg38", "hg19"),
                            chrs = paste0("chr", 1:22),
                            ideogram_height = 1) {
  .check_install("circlize")
  .check_install("scales")
  genome_build <- match.arg(genome_build)
  target_genes <- readRDS(file.path(
    system.file("extdata", package = "gcap"),
    paste0(genome_build, "_target_genes.rds")
  ))
  if (!startsWith(fCNA$data$gene_id[1], "ENSG")) {
    message("detected you have transformed ENSEMBL ID, also transforming gene annotation data")
    opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
    options(IDConverter.datapath = opts)
    target_genes$gene_id <- IDConverter::convert_hm_genes(target_genes$gene_id, genome_build = genome_build)
    target_genes <- target_genes[!is.na(target_genes$gene_id), ]
  }
  data_bed <- merge(fCNA$data, target_genes, by = "gene_id", all.x = TRUE, sort = FALSE)
  data_bed$amplicon_type <- ifelse(data_bed$amplicon_type %in% c("circular", "possibly_circular"), "circular", "noncircular")
  data_bed <- data_bed[!is.na(data_bed$gene_id) & data_bed$chrom %in% chrs,
    c(
      "chrom", "start", "end", "total_cn", "gene_id",
      "amplicon_type"
    ),
    with = FALSE
  ]
  colnames(data_bed)[1] <- "chr"
  if (nrow(data_bed) < 1) {
    message("no data to plot, please check!")
    return(invisible(NULL))
  }

  cnrange <- range(data_bed$total_cn, na.rm = TRUE)
  # nsamples <- nrow(fCNA$sample_summary)
  bed_list <- split(data_bed, data_bed$amplicon_type)
  bed_list <- lapply(bed_list, function(x) {
    x[, .(freq = .N), by = .(chr, start, end, gene_id)]
  })
  bed_cn <- data.table::dcast(data_bed[, .(total_cn = mean(total_cn, na.rm = T)),
    by = .(chr, start, end, gene_id, amplicon_type)
  ],
  chr + start + end + gene_id ~ amplicon_type,
  value.var = "total_cn"
  )
  bed_cn[, circular := ifelse(is.na(circular), 0, circular)]
  bed_cn[, noncircular := ifelse(is.na(noncircular), 0, -noncircular)]

  gap_after <- c(rep(1, length(chrs) - 1), 12)
  circlize::circos.par(
    "start.degree" = 90,
    "track.height" = 0.15, # % of the circle radius
    "gap.after" = gap_after
  )

  track_height <- 0.1
  circlize::circos.initializeWithIdeogram(
    species = genome_build,
    chromosome.index = chrs,
    track.height = circlize::convert_height(
      if (!is.null(highlight_genes)) 0 else track_height,
      "mm"
    ), ideogram.height = circlize::convert_height(
      if (!is.null(highlight_genes)) 0 else ideogram_height,
      "mm"
    ),
    plotType = if (!is.null(highlight_genes)) NULL else c("ideogram", "axis", "labels")
  )
  on.exit(circlize::circos.clear())

  if (!is.null(highlight_genes)) {
    # circlize::circos.genomicInitialize(
    #   circlize::read.cytoband(species = genome_build, chromosome.index = chrs)$df,
    #   plotType = "labels"
    # )
    if (is.data.frame(highlight_genes)) {
      message("found input a data.frame for highlight genes")
      stopifnot("gene_id" %in% colnames(highlight_genes))
      data <- data.table::as.data.table(highlight_genes)
      highlight_genes <- data$gene_id
    } else {
      data <- data.table::data.table()
    }

    bed <- unique(data_bed[
      data_bed$gene_id %in% highlight_genes,
      c("chr", "start", "end", "gene_id")
    ])
    if (ncol(data) > 0) {
      bed <- merge(bed, data, by = "gene_id", all.x = TRUE, sort = FALSE)
      data.table::setcolorder(bed, c("chr", "start", "end", "gene_id"))
    }

    if (nrow(bed) < 1) {
      message("no data for your selected genes, please check")
      return(invisible(NULL))
    }

    if (is.null(clust_distance) & ncol(data) == 0) {
      bed_col <- as.numeric(factor(bed$gene_id))
    } else if ("cluster" %in% colnames(bed)) {
      message("found 'cluster' column, use it for color mapping")
      bed_col <- as.numeric(factor(bed$cluster))
    } else {
      message("clustering highlight genes by 'hclust' average method with distance data from gene centers")
      message("distance cutoff: ", clust_distance)
      bed_col <- clusterGPosition(bed, clust_distance)
      message("done")
    }

    if ("label" %in% colnames(bed)) {
      message("label detected, add it to gene id")
      bed$gene_id <- paste0(bed$gene_id, " (", bed$label, ")")
    }

    ssize <- max(nchar(highlight_genes)) / 8
    circlize::circos.genomicLabels(bed,
      labels = bed$gene_id, side = "outside",
      cex = 0.5 / (ssize),
      connection_height = circlize::mm_h(5 / ssize),
      col = bed_col,
      line_col = bed_col
    )
    circlize::circos.genomicIdeogram(
      species = genome_build,
      track.height = 0.02 * ideogram_height
    )
    circlize::circos.track(
      track.index = circlize::get.current.track.index(),
      panel.fun = function(x, y) {
        circlize::circos.genomicAxis(h = "top", direction = "outside")
      }
    )
  }

  draw_freq_track(bed_list$circular, col[1])
  draw_track_y_axis(chrs)
  draw_freq_track(bed_list$noncircular, col[2])
  draw_track_y_axis(chrs)
  draw_bi_track(bed_cn, col)
  draw_track_y_axis(chrs, at = scales::breaks_pretty(5)(c(0, cnrange[2])))
}

draw_track_y_axis <- function(chrs, at = NULL) {
  circlize::circos.yaxis(
    at = at,
    labels.cex = 0.3,
    tick.length = circlize::convert_x(
      0.3,
      "mm", circlize::get.cell.meta.data("sector.index"),
      circlize::get.cell.meta.data("track.index")
    ),
    side = "right", sector.index = chrs[length(chrs)]
  )
}

draw_freq_track <- function(data, col) {
  data$gene_id <- NULL # Only keep freq as value
  circlize::circos.genomicTrack(
    data,
    ylim = c(0, max(data$freq, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      # circlize::circos.genomicRect(region, value,
      #   ytop.column = 1,
      #   ybottom = 0,
      #   col = col,
      #   border = NA
      # )
      circlize::circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = col, border = NA)
    }
  )
}

draw_bi_track <- function(data, col) {
  circlize::circos.genomicTrack(
    data,
    # track.height = if (!is.null(highlight_genes)) 0.2 else 0.3, # ylim = c(0, 0.1),
    panel.fun = function(region,
                         value, ...) {
      circlize::circos.genomicPoints(region, value, pch = 16, cex = 0.3, col = col, border = NA)

      # circlize::circos.genomicRect(
      #   region, value,
      #   ytop = value$circular,
      #   ybottom = value$noncircular,
      #   col = ifelse(value$circular > 0, col[1], col[2]), border = NA
      # )
      circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#000040")
    }
  )
}
