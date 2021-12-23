#' Determine sample-level ecDNA status with aggregated probability
#'
#' Aggregated probability is calculated from gene prob or cluster prob
#' (this first cluster genes and then aggregate the cluster prob).
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' The column storing prediction result must start with `pred`.
#' @param min_ngene,min_ncluster minimal gene or cluster number hosted on
#' amplicon.
#' @return a `data.table`.
#' - `*_gene_based_prob` for gene based probability.
#' - `*_cluster_based_prob` for cluster based probability.
#' We cluster genes based on  distance (hard cutoff `1e6` (i.e., 1Mb)).
#' For genes from different
#' chromosomes, the distance set to `1e8` (i.e., 100Mb).
#' @export
#'
#' @examples
#' data("ec")
#' ec2 <- ec
#' ec2$pred <- gcap.runPrediction(ec)
#' score <- gcap.aggrProb(ec2)
#' score
#' @testexamples
#' expect_equal(nrow(score), 1L)
gcap.aggrProb <- function(data,
                          genome_build = c("hg38", "hg19"),
                          min_ngene = 2L,
                          min_ncluster = 2L) {
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

  lg$info("estimating gene numbers with prob")
  dt_ng <- data[, lapply(scs, function(x) {
    estimate_n(.SD[[x]])
  }), by = "sample"]
  colnames(dt_ng)[-1] <- paste0(scs, "_ngenes")

  lg$info("aggregating prob over genes")
  dt_genes <- data[, lapply(scs, function(x) {
    calc_prob(.SD[[x]], n = min_ngene)
  }), by = "sample"]
  colnames(dt_genes)[-1] <- paste0(scs, "_gene_based_prob")

  lg$info("aggregating prob over clusters (<1e6 bp from center)")

  dt_clusters <- data[, lapply(scs, function(x) {
    calc_clusters2(.SD[, c(x, "gene_id"), with = FALSE], genome_build, n = min_ncluster)
  }), by = "sample"]
  colnames(dt_clusters)[-1] <- paste0(scs, paste0("_cluster_based_prob"))

  lg$info("merging and outputing final data")
  mergeDTs(list(dt_ng, dt_genes, dt_clusters), by = "sample")
}

calc_clusters2 <- function(dt, genome_build, n = 2) {
  dt <- dt[dt[[1]] > 0.01] # set a minimal prob cutoff
  colnames(dt)[1] <- "prob"
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
    y_all <- na.omit(y[, .(chr, center, prob)])
    y <- as.matrix(y[, .(chr, center)])
    if (nrow(y) < n) {
      0
    } else {
      dst <- calc_dist(y)
      cls <- stats::cutree(stats::hclust(stats::as.dist(dst), method = "average"), h = 1e6)
      y_all$cls <- cls
      y_all <- y_all[, .(prob = mean(prob, na.rm = TRUE)), by = "cls"]
      calc_prob(y_all$prob, n)
    }
  } else {
    0
  }
}


# https://stats.stackexchange.com/questions/41247/risk-of-extinction-of-schr%C3%B6dingers-cats
convolve.binomial <- function(p, r = FALSE) {
  # p is a vector of probabilities of Bernoulli distributions.
  # The convolution of these distributions is returned as a vector
  # `z` where z[i] is the probability of i-1, i=1, 2, ..., length(p)+1.
  if (r) {
    n <- length(p) + 1
    z <- c(1, rep(0, n - 1))
    for (pi in p) z <- (1 - pi) * z + pi * c(0, z[-n])
    return(z)
  } else {
    conv_binomial(p)
  }
}

# set.seed(123)
# x = abs(runif(10000))
# r = microbenchmark::microbenchmark(
#   rversion = convolve.binomial(x, r = TRUE),
#   rcpp = convolve.binomial(x), times = 10, check = "equal"
# )
# r
# plot(r)

# n refer to occur >= n
calc_prob <- function(p, n) {
  p <- convolve.binomial(p)
  p <- 1 - sum(p[seq_len(n)])
  if (is.na(p)) 0 else p
}

estimate_n <- function(p, level = 0.01) {
  p <- convolve.binomial(p)
  z <- cumsum(p) < level
  if (any(z)) {
    max(which(z))
  } else 0L
}
