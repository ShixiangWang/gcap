#' Summarize prediction result into gene/sample-level
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' @param min_prob the minimal aggregated (in cytoband level) probability to determine a circular amplicon.
#' @param tightness a coefficient to times to TCGA somatic CN to set a more strict threshold
#' as a circular amplicon.
#' If the value is larger, it is more likely a fCNA assigned to `noncircular`
#' instead of `circular`. **When it is `NA`, we don't use TCGA somatic CN data as reference**.
#' @param gap_cn a gap copy number value, default is `4L` refer to Kim 2020 Nat.Gen.
#' A gene with copy number above `ploidy + gap_cn` would be treated as focal amplicon.
#' Smaller, more amplicons.
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
                            min_prob = 0.7,
                            tightness = 1L,
                            gap_cn = 4L) {
  on.exit(invisible(gc()))
  stopifnot(is.data.frame(data))
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
  somatic_cn <- readRDS(system.file("extdata", "somatic_gene_cn.rds", package = "gcap"))[, 1:3]
  # colnames(somatic_cn)[2:3] <- c("somatic_cn_mean", "somatic_cn_sd")
  somatic_cn <- somatic_cn[order(somatic_cn$somatic_cn_mean)]

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

  all_ids <- unique(data$gene_id)
  df_ids <- setdiff(all_ids, somatic_cn$gene_id)
  if (length(df_ids) > 0) {
    somatic_cn <- rbind(somatic_cn, data.table::data.table(
      gene_id = df_ids,
      somatic_cn_mean = somatic_cn$somatic_cn_mean[1],
      somatic_cn_sd = somatic_cn$somatic_cn_sd[1]
    ))
  } # Fill with smallest value
  data <- merge(data, somatic_cn, by = "gene_id", all.x = TRUE)

  if ("band" %in% colnames(data)) data$band <- NULL
  data <- merge(data, gene_cytobands, by = "gene_id", all.x = TRUE)

  # cytoband_cn <- data[
  #   , .(cytoband_cn_median = median(total_cn, na.rm = TRUE)),
  #   by = .(sample, band)
  # ] # Currently, median is not used
  # data <- merge(data, cytoband_cn, by = c("sample", "band"), all.x = TRUE)

  # (background_cn) Mean + SD for large threshold (circular), smaller is looser
  # (background_cn2) Ploidy for basic threshold (above noncircular)
  data$background_cn <- pmax((data$somatic_cn_mean + data$somatic_cn_sd) * tightness, data$ploidy, na.rm = TRUE)
  data$background_cn <- ifelse(is.na(data$background_cn), 2, data$background_cn)
  data$background_cn2 <- data$ploidy
  data$background_cn2 <- ifelse(is.na(data$background_cn2), 2, data$background_cn2)

  data$somatic_cn_mean <- NULL
  data$somatic_cn_sd <- NULL

  flag_amp <- data$total_cn >= (data$background_cn + gap_cn) * pmax(data$ploidy, 2, na.rm = TRUE) / 2
  flag_amp2 <- data$total_cn >= data$background_cn2 + gap_cn
  if (is.na(tightness) || tightness == 0) {
    # In such case, remove limit from reference somatic CN
    flag_amp <- flag_amp2
    data$background_cn <- data$background_cn2
  }
  flag_circle <- data$prob > 0.5
  # Classify amplicon
  data$gene_class <- data.table::fcase(
    flag_amp & flag_circle, "circular",
    flag_amp2, "noncircular",
    default = "nofocal"
  )
  #data$gene_class <- factor(data$gene_class, levels = c("nofocal", "noncircular", "circular"))

  # Generate fCNA and output
  sel_cols <- c(
    "sample", "purity", "ploidy", "AScore", "pLOH", "cna_burden",
    paste0("CN", 1:19)
  )
  sel_cols <- sel_cols[sel_cols %in% colnames(data)]
  pdata <- data[match(unique(data$sample), sample), sel_cols, with = FALSE]
  lg$info("only keep genes labeled as amplicons in result fCNA object")
  fcna <- data[!gene_class %in% "nofocal",
    c(
      "sample", "band", "gene_id", "total_cn",
      "minor_cn", "ploidy", "prob", "gene_class"
    ),
    with = FALSE
  ]

  data.table::setkey(fcna, NULL) # To make sure all equal to rebuild the fCNA object from file
  if (nrow(fcna) == 0) lg$info("No fCNA records detected")
  fCNAobj <- fCNA$new(fcna, pdata, min_prob = min_prob)
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
