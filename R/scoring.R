#' Summarize prediction result into gene/sample-level
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' @param min_n a minimal cytoband number (default is `1`) to determine
#' sample class. e.g., sample with at least 1 cytoband harboring circular
#' genes would be labelled as "circular".
#' @param tightness If the value is larger, it is more likely a fCNA assigned to `noncircular`
#' instead of `circular`. Range from 1 to Inf.
#' @param gap_cn a gap copy number value, default is `4L` refer to Kim 2020 Nat.Gen.
#' A gene with copy number above `ploidy + gap_cn` would be treated as amplicon.
#' Smaller, more amplicons.
#' @param circ_prob a prob cutoff value to treat a gene as located at circular amplicon.
#' @return a list of `data.table`.
#' @export
#'
#' @examples
#' data("ec")
#' ec2 <- ec
#' ec2$prob <- gcap.runPrediction(ec)
#' score <- gcap.runScoring(ec2)
#' score
#' @testexamples
#' expect_equal(length(score), 3L)
gcap.runScoring <- function(data,
                            genome_build = "hg38",
                            min_n = 1L,
                            tightness = 1L,
                            gap_cn = 4L,
                            circ_prob = 0.95) {
  on.exit(invisible(gc()))
  stopifnot(is.data.frame(data), circ_prob > 0.15, circ_prob <= 1)
  lg <- set_logger()

  lg$info("checking input data type")
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  lg$info("checking columns")
  if (!"prob" %in% colnames(data)) {
    lg$fatal("no columns 'prob' found, please check your input")
    stop()
  }
  if (!"sample" %in% colnames(data)) {
    lg$info("No 'sample' column found, mutating one")
    data$sample <- "sample"
  }

  lg$info("filtering out records without prob result")
  data <- data[!is.na(data$prob)]

  lg$info("joining extra annotation data")
  cytobands <- as.data.table(
    sigminer::get_genome_annotation("cytobands", genome_build = genome_build)
  )[chrom %in% paste0("chr", 1:22)]
  cytobands$stain <- NULL
  cytobands$start <- cytobands$start + 1L

  ref_file <- system.file(
    "extdata", paste0(genome_build, "_target_genes.rds"),
    package = "gcap", mustWork = TRUE
  )
  refgene <- readRDS(ref_file)
  gene_cytobands <- overlaps(refgene, cytobands)
  gene_cytobands[, intersect_ratio := intersect_size / (i.end - i.start + 1)]
  gene_cytobands <- gene_cytobands[, .SD[which.max(intersect_ratio)], by = .(gene_id)][order(chr, i.start)]
  gene_cytobands <- gene_cytobands[, .(gene_id, band = paste(chr, band, sep = ":"))]

  if ("band" %in% colnames(data)) data$band <- NULL
  data <- merge(data, gene_cytobands, by = "gene_id", all.x = TRUE)

  cytoband_cn <- data[
    , .(cytoband_cn_median = median(total_cn, na.rm = TRUE)),
    by = .(sample, band)
  ] # Currently, median is not used
  data <- merge(data, cytoband_cn, by = c("sample", "band"), all.x = TRUE)

  # (background_cn) large threshold (circular)
  # (background_cn2) small threshold (noncircular)
  data$background_cn <- data$ploidy * max(tightness, 1, na.rm = TRUE)
  data$background_cn <- ifelse(is.na(data$background_cn), 2, data$background_cn)
  data$background_cn2 <- data$ploidy
  data$background_cn2 <- ifelse(is.na(data$background_cn2), 2, data$background_cn2)

  flag_amp <- data$total_cn >= data$background_cn + gap_cn * pmax(data$ploidy, 2, na.rm = TRUE) / 2
  flag_amp2 <- data$total_cn >= data$background_cn2 + gap_cn

  flag_circle <- as.integer(cut(data$prob, breaks = c(0, 0.15, circ_prob, 1), include.lowest = TRUE))
  # Classify amplicon
  data$amplicon_type <- data.table::fcase(
    flag_amp & flag_circle == 3, "circular",
    flag_amp & flag_circle == 2, "possibly_circular",
    flag_amp | flag_amp2, "noncircular",
    default = "nofocal"
  )
  data$amplicon_type <- factor(data$amplicon_type, levels = c(
    "nofocal", "noncircular", "possibly_circular", "circular"
  ))

  # Generate input of fCNA class
  sel_cols <- c(
    "sample", "purity", "ploidy", "AScore", "pLOH", "cna_burden",
    paste0("CN", 1:19)
  )
  sel_cols <- sel_cols[sel_cols %in% colnames(data)]
  pdata <- data[match(unique(data$sample), sample), sel_cols, with = FALSE]
  fcna <- data[!amplicon_type %in% c(NA, "nofocal"),
    c(
      "sample", "band", "gene_id", "total_cn",
      "minor_cn", "background_cn", "prob", "amplicon_type"
    ),
    with = FALSE
  ] # Only treat nofocal as (focal) variants

  data.table::setkey(fcna, NULL) # To make sure all equal to rebuild the fCNA object from file
  if (nrow(fcna) == 0) lg$info("No fCNA records detected")
  fCNAobj <- fCNA$new(fcna, pdata, min_n = min_n)
  print(fCNAobj)

  lg$info("done")
  list(
    data = data,
    fCNA = fCNAobj
  )
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
