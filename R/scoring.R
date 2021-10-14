#' Quantify sample-level ecDNA measueres
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' The column storing prediction result must start with `pred`.
#' @param cutoff a cutoff for converting prob into `0/1` value.
#' @return a `data.table`.
#' - `*_load` for number of gene affected.
#' - `*_burden` for size of genome (per Mb) affected.
#' - `*_clusters` for number of region affected by clustering
#' distance with hard cutoff `1e6` (i.e., 1Mb). For genes from different
#' chromosomes, the distance set to `1e7` (i.e., 10Mb).
#' @export
#'
#' @examples
#' data("ec")
#' ec2 <- ec
#' ec2$pred <- gcap.runPrediction(ec)
#' score <- gcap.runScoring(ec2)
#' score
#' @testexamples
#' expect_equal(nrow(score), 1L)
gcap.runScoring <- function(data, cutoff = 0.9,
                            genome_build = c("hg38", "hg19")) {
  stopifnot(is.data.frame(data))
  genome_build <- match.arg(genome_build)
  lg <- set_logger()

  lg$info("checking input data type")
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  lg$info("checking columns")
  scs <- colnames(data)[startsWith(colnames(data), "pred")]
  if (length(scs) < 1) {
    lg$fatal("no columns start with 'pred' found, please check your input")
    stop()
  }
  if (!"sample" %in% colnames(data)) {
    lg$info("No 'sample' column found, mutating one")
    data$sample <- "sample"
  }

  lg$info("converting prob")
  scs2 <- paste0(scs, "_binary")
  for (i in seq_along(scs)) {
    data[[scs2[i]]] <- ifelse(data[[scs[i]]] > cutoff, 1L, 0L)
  }

  lg$info("calculating load")
  dt_load <- data[, lapply(.SD, sum), .SDcols = scs2, by = "sample"]
  colnames(dt_load) <- sub("binary", "load", colnames(dt_load))

  lg$info("calculating burden")
  scs_burden <- paste0(scs, "_burden")
  dt_burden <- data[, lapply(scs2, function(x) {
    calc_burden(.SD[, c(x, "gene_id"), with = FALSE], genome_build)
  }), by = "sample"]
  colnames(dt_burden)[-1] <- scs_burden

  lg$info("calculating clusters based on distance from gene center with cutoff 1e6")
  scs_clusters <- paste0(scs, "_clusters")
  dt_clusters <- data[, lapply(scs2, function(x) {
    calc_clusters(.SD[, c(x, "gene_id"), with = FALSE], genome_build)
  }), by = "sample"]
  colnames(dt_clusters)[-1] <- scs_clusters

  lg$info("merging and outputing final data")
  mergeDTs(list(dt_load, dt_burden, dt_clusters), by = "sample")
}

calc_burden <- function(dt, genome_build) {
  dt <- dt[dt[[1]] == 1L]
  if (nrow(dt) > 0) {
    ref_file <- system.file(
      "extdata", paste0(genome_build, "_target_genes.rds"),
      package = "gcap", mustWork = TRUE
    )
    y <- readRDS(ref_file)
    colnames(y)[1:3] <- c("chr", "start", "end")
    y <- merge(dt, y, by = "gene_id", all.x = TRUE)
    y[, len := end - start + 1L]
    sum(y$len / 1e6)
  } else {
    0
  }
}

calc_clusters <- function(dt, genome_build) {
  dt <- dt[dt[[1]] == 1L]
  if (nrow(dt) > 0) {
    ref_file <- system.file(
      "extdata", paste0(genome_build, "_target_genes.rds"),
      package = "gcap", mustWork = TRUE
    )
    y <- readRDS(ref_file)
    colnames(y)[1:3] <- c("chr", "start", "end")
    y <- merge(dt, y, by = "gene_id", all.x = TRUE)
    # calculate the distance
    y$chr <- gsub("X", "23", y$chr, ignore.case = TRUE)
    y$chr <- gsub("Y", "24", y$chr, ignore.case = TRUE)
    y$chr <- as.integer(gsub("chr", "", y$chr))
    y$center <- round((y$end - y$start) / 2)
    y <- na.omit(as.matrix(y[, .(chr, center)]))
    if (nrow(y) < 2) {
      nrow(y)
    } else {
      dst <- calc_dist(y)
      length(table(stats::cutree(stats::hclust(stats::as.dist(dst), method = "average"), h = 1e6)))
    }
  } else {
    0L
  }
}
