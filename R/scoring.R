#' Quantify sample-level ecDNA measueres
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' The column storing prediction result must start with `pred`.
#' @return a `data.table`.
#' - `*_load` for number of gene affected.
#' - `*_burden` for size of genome (per Mb) affected.
#' - `*_clusters` for number of region affected by clustering
#' distance with hard cutoff `1e6` (i.e., 1Mb). For genes from different
#' chromosomes, the distance set to `1e8` (i.e., 10Mb).
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
gcap.runScoring <- function(data,
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

  lg$info("using 0.1, 0.5 and 0.9 as cutoffs for low, medium and high-level thresholded amplicon genes")
  scs2 <- paste0(scs, "_levels")
  for (i in seq_along(scs)) {
    data[[scs2[i]]] <- cut(data[[scs[i]]], breaks = c(0.1, 0.5, 0.9, 1), labels = c("low", "medium", "high"))
  }

  lg$info("calculating load with different thresholds")
  th_levels <- c("low", "medium", "high")

  dt_load <- list()
  for (i in th_levels) {
    dt_load[[i]] <- data[, lapply(.SD, function(x) {
      if (i == "low") lvls <- c("low", "medium", "high")
      if (i == "medium") lvls <- c("medium", "high")
      if (i == "high") lvls <- "high"
      sum(x %in% lvls, na.rm = TRUE)
    }),
    .SDcols = scs2, by = "sample"
    ]
    colnames(dt_load[[i]]) <- sub("levels", paste0("load_", i, "_thresholded"), colnames(dt_load[[i]]))
  }
  dt_load <- Reduce(merge, dt_load)

  lg$info("calculating burden with different thresholds")
  dt_burden <- list()
  for (i in th_levels) {
    dt_burden[[i]] <- data[, lapply(scs2, function(x) {
      calc_burden(.SD[, c(x, "gene_id"), with = FALSE], genome_build, th = i)
    }), by = "sample"]
    scs_burden <- paste0(scs, paste0("_burden_", i, "_thresholded"))
    colnames(dt_burden[[i]])[-1] <- scs_burden
  }
  dt_burden <- Reduce(merge, dt_burden)


  lg$info("detecting clusters (<1e6 bp from center) with different thresholds")
  dt_clusters <- list()
  for (i in th_levels) {
    dt_clusters[[i]] <- data[, lapply(scs2, function(x) {
      calc_clusters(.SD[, c(x, "gene_id"), with = FALSE], genome_build, th = i)
    }), by = "sample"]
    scs_clusters <- paste0(scs, paste0("_clusters_", i, "_thresholded"))
    colnames(dt_clusters[[i]])[-1] <- scs_clusters
  }
  dt_clusters <- Reduce(merge, dt_clusters)

  lg$info("merging and outputing final data")
  mergeDTs(list(dt_load, dt_burden, dt_clusters), by = "sample")
}

calc_burden <- function(dt, genome_build, th = "low") {
  if (th == "low") lvls <- c("low", "medium", "high")
  if (th == "medium") lvls <- c("medium", "high")
  if (th == "high") lvls <- "high"
  dt <- dt[dt[[1]] %in% lvls]
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

calc_clusters <- function(dt, genome_build, th = "low") {
  if (th == "low") lvls <- c("low", "medium", "high")
  if (th == "medium") lvls <- c("medium", "high")
  if (th == "high") lvls <- "high"
  dt <- dt[dt[[1]] %in% lvls]
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
