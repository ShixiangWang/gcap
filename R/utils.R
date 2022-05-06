.check_install <- function(pkg, bioc = FALSE, ...) {
  lg <- set_logger()
  install_func <- if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) {
    lg$info("installing required package {pkg}")
    install_func(pkg, ...)
  }
}

#' Cluster genomic regions by distance from region center
#' @param dt a genomic region data with first 3 columns for
#' chromosome, start and end.
#' @param distance a distance value for clustering,
#' regions with average distance lower than this value would be clustered
#' into one. Should not greater than 1e8.
#' @param no_annotation if `TRUE`, use builtin data for annotation. This is
#' useful when input a list of genes.
#' @param genome_build genome build version.
#' @param simplify if `FALSE`, return a `data.table` instead of a vector.
#' @return cluster list.
#' @export
#' @examples
#' library(data.table)
#' region <- data.table(
#'   chr = c("chr1", "chr2", "chr1"),
#'   start = c(1032, 992, 10000), end = c(4242, 1e6, 400023)
#' )
#' clusterGPosition(region, distance = 1000)
#' clusterGPosition(region)
#' clusterGPosition(data.table(gene_id = c("KRAS", "EGFR", "MYC", "ERBB2", "GRB7")),
#'   no_annotation = TRUE, simplify = FALSE
#' )
#' clusterGPosition(data.table(gene_id = c("KRAS", "EGFR", "MYC", "ERBB2", "GRB7")),
#'   no_annotation = TRUE, simplify = TRUE
#' )
clusterGPosition <- function(dt, distance = 1e7,
                             no_annotation = FALSE,
                             genome_build = c("hg38", "hg19"),
                             simplify = TRUE) {
  dt <- if (data.table::is.data.table(dt)) data.table::copy(dt) else data.table::as.data.table(dt)
  genome_build <- match.arg(genome_build)
  if (no_annotation) {
    stopifnot("gene_id" %in% colnames(dt))
    target_genes <- readRDS(file.path(
      system.file("extdata", package = "gcap"),
      paste0(genome_build, "_target_genes.rds")
    ))
    if (!startsWith(dt$gene_id[1], "ENSG")) {
      message("detected you have transformed ENSEMBL ID, also transforming gene annotation data")
      opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
      options(IDConverter.datapath = opts)
      target_genes$gene_id <- IDConverter::convert_hm_genes(target_genes$gene_id, genome_build = genome_build)
      target_genes <- target_genes[!is.na(target_genes$gene_id), ]
    }
    dt <- merge(dt, target_genes, by = "gene_id", all.x = TRUE, sort = FALSE)
    data.table::setcolorder(dt, c("chrom", "start", "end"))
  }
  colnames(dt)[1:3] <- c("chr", "start", "end")
  dt[, `:=`(x = as.integer(factor(chr)), y = (start + end) / 2)]
  dst <- stats::as.dist(calc_dist(as.matrix(dt[, .(x, y)])))
  cls <- stats::cutree(stats::hclust(dst, method = "average"), h = distance)
  if (simplify) {
    cls
  } else {
    dt$cluster <- cls
    dt
  }
}

#' Merge a list of data.table
#' @param dt_list a list of `data.table`s.
#' @param by which column used for merging.
#' @param sort should sort the result?
#' @return a `data.table`
#' @export
mergeDTs <- function(dt_list, by = NULL, sort = FALSE) {
  dt_list <- dt_list[lengths(dt_list) != 0]
  Reduce(
    function(...) {
      merge(..., by = by, all = TRUE, sort = sort)
    }, dt_list
  )
}

check_model <- function(m) {
  if (is.character(m)) {
    opts <- c("XGB11", "XGB32", "XGB56", "toy")
    if (length(m) > 1) stop("Multiple inputs for model are not allowed")
    if (!m %in% opts) {
      stop("Internal model can only set to one of ", paste(opts, collapse = ", "))
    }
  }
}

#' Get overlaps of two genomic regions
#' @param x,y a genemic region with data.frame format, the first 3 columns
#' should representing chromosome, start and end position.
#' @return a `data.table`
#' @export
overlaps <- function(x, y) {
  ## Overlaps genome regions from x and y
  if (!is.data.frame(x)) {
    stop("x must be a data.frame")
  }
  if (!is.data.frame(y)) {
    stop("y must be a data.frame")
  }

  x <- data.table::as.data.table(x)
  y <- data.table::as.data.table(y)

  colnames(x)[1:3] <- colnames(y)[1:3] <- c("chr", "start", "end")

  # Determine the intersect size
  data.table::setkey(y, chr, start, end)
  out <- data.table::foverlaps(x, y)[!is.na(start)]

  out[, `:=`(intersect_size, sigminer::get_intersect_size(i.start, i.end, start, end))]
  # Calculate the region cov ratio
  out[, `:=`(intersect_ratio, intersect_size / (abs(end - start) + 1))]
  out
}

#' Set default factor level for fCNA class
#'
#' @param class a vector of fCNA class.
#' @param ref_level the reference level of factor.
#'
#' @return a vector.
#' @export
#'
#' @examples
#' set_default_factor(c("nofocal", "noncircular", "circular", "nofocal", "possibly_circular"))
set_default_factor <- function(class, ref_level = "nofocal") {
  class <- fcase(class %in% c("circular", "possibly_circular"), "circular",
    class == "noncircular", "noncircular",
    default = "nofocal"
  )
  factor(class, levels = c(ref_level, setdiff(c("nofocal", "noncircular", "circular"), ref_level)))
}

# Global variables
utils::globalVariables(
  c(
    ".", "old_sample", "ploidy", "AScore", "cna_burden",
    "age", "chr", "freq_BFB", "freq_Circular", "freq_HR",
    "freq_Linear", "gender", "gene_id", "i.end", "i.start",
    "intersect_ratio", "intersect_size", "minor_cn", "total_cn",
    "len", "center", "ploidy2",
    "eval_aucpr", "prob", "train_aucpr", "background_cn",
    "band", "chrom", "cytoband_cn_median", "status",
    "N_pos", "amplicon_type", "ec_genes", "total_cn_neg", "total_cn_pos",
    "N", "x", "y", "circular", "noncircular"
  )
)
