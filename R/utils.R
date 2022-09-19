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
    opts <- c("XGB11", "XGB32", "XGB56",
              "toy")
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

# Global variables
utils::globalVariables(
  c(
    ".", "old_sample", "ploidy", "AScore", "cna_burden",
    "age", "chr", "freq_BFB", "freq_Circular", "freq_HR",
    "freq_Linear", "gender", "gene_id", "i.end", "i.start",
    "intersect_ratio", "intersect_size", "minor_cn", "total_cn",
    "len", "center", "ploidy2",
    "eval_aucpr", "prob", "train_aucpr", "background_cn",
    "band", "chrom", "cytoband_cn_median", "gene_class", "thlevel"
  )
)
