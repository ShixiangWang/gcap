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

# Copy from sigminer
query_remote_data <- function(x) {
  x_url <- paste0("https://zenodo.org/record/4771552/files/", x)
  x_dest <- file.path(system.file("extdata", package = "sigminer"), x)
  message("Downloading ", x_url, " to ", x_dest)
  tryCatch(
    {
      download.file(
        url = x_url,
        destfile = x_dest
      )
      TRUE
    },
    error = function(e) {
      warning("Failed downloading the data.", immediate. = TRUE)
      FALSE
    }
  )
}

get_ref_data <- function(genome_build = c("hg38", "hg19", "mm10", "mm9")) {
  genome_build <- match.arg(genome_build)
  gene_file <- switch(genome_build,
    mm9 = file.path(
      system.file("extdata", package = "sigminer"),
      "mouse_mm9_gene_info.rds"
    ),
    mm10 = file.path(
      system.file("extdata", package = "sigminer"),
      "mouse_mm10_gene_info.rds"
    ),
    file.path(
      system.file("extdata", package = "sigminer"),
      paste0("human_", genome_build, "_gene_info.rds")
    )
  )
  ok <- TRUE
  if (!file.exists(gene_file)) ok <- query_remote_data(basename(gene_file))
  if (!ok) {
    return(invisible(NULL))
  }
  gene_dt <- readRDS(gene_file)
  gene_dt
}

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
    "band", "chrom", "cytoband_cn_median", "status",
    "N_pos", "amplicon_type", "ec_genes", "total_cn_neg", "total_cn_pos"
  )
)
