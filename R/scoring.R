#' Summarize prediction result into gene/sample-level
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' @param min_n a minimal cytoband number (default is `1`) to determine
#' sample class. e.g., sample with at least 1 cytoband harboring circular
#' genes would be labelled as "circular".
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
                            min_n = 1L) {
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
  blood_cn <- readRDS(system.file("extdata", "blood_gene_cn.rds", package = "gcap"))
  colnames(blood_cn)[2] <- c("blood_cn_median")

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

  data <- merge(data, blood_cn, by = "gene_id", all.x = TRUE)
  data$blood_cn_median <- ifelse(is.na(data$blood_cn_median), 2, data$blood_cn_median)
  data <- merge(data, gene_cytobands, by = "gene_id", all.x = TRUE)

  cytoband_cn <- data[, .(cytoband_cn_median = median(total_cn, na.rm = TRUE)), by = .(sample, band)]
  data <- merge(data, cytoband_cn, by = c("sample", "band"), all.x = TRUE)

  data$background_cn <- data$blood_cn_median * data$ploidy / 2
  data$blood_cn_median <- NULL

  # Similar to NG paper
  # As a prerequisite, amplicons must at least 4 copies above background CN to be considered a valid amplicon.
  flag_amp <- data$total_cn >= data$background_cn + 4
  flag_circle <- as.integer(cut(data$prob, breaks = c(0, 0.15, 0.75, 1), include.lowest = TRUE))
  # Classify amplicon
  data$amplicon_type <- data.table::fcase(
    flag_amp & flag_circle == 3, "circular",
    flag_amp & flag_circle == 2, "possibly_circular",
    flag_amp, "noncircular",
    default = "nofocal"
  )
  data$amplicon_type <- factor(data$amplicon_type, levels = c(
    "nofocal", "noncircular", "possibly_circular", "circular"
  ))

  lg$info("summarizing in gene level")
  # Gene level summary
  dt_genes <- data[, .SD[sum(amplicon_type %in% c("possibly_circular", "circular")) > 0], by = .(gene_id)][
    , .(
      N = .N,
      total_cn = mean(total_cn, na.rm = TRUE),
      minor_cn = mean(minor_cn, na.rm = TRUE),
      cytoband_cn = mean(cytoband_cn_median, na.rm = TRUE)
    ),
    by = .(gene_id, amplicon_type)
  ]

  dt_genes <- merge(dt_genes, data[match(unique(dt_genes$gene_id), gene_id), .(gene_id, band)],
    by = "gene_id", all.x = TRUE
  )
  dt_genes <- dt_genes[order(band, -total_cn)]

  data$cytoband_cn_median <- NULL

  lg$info("summarizing in sample level")
  # Sample level summary
  # Number of ec Genes and sample classification
  summarize_sample <- function(data) {
    ec_genes <- sum(data$amplicon_type == "circular", na.rm = TRUE)
    ec_cytobands_detail <- unique(data$band[data$amplicon_type == "circular"])
    ec_cytobands <- length(ec_cytobands_detail)
    ec_cytobands_detail <- paste(sort(ec_cytobands_detail), collapse = ",")
    ec_possibly_genes <- sum(data$amplicon_type == "possibly_circular", na.rm = TRUE)
    # prob_possibly <- data$prob[data$amplicon_type %in% c("possibly_circular", "circular")]
    prob_possibly <- data[data$amplicon_type %in% c("possibly_circular", "circular")]
    if (nrow(prob_possibly) > 0) {
      prob_possibly <- prob_possibly[
        , .(prob = max(prob, na.rm = TRUE)),
        by = .(band)]$prob
    } else {
      prob_possibly <- NULL
    }

    # flags have priority
    # flag_ec <- ec_genes >= min_n
    flag_ec <- ec_cytobands >= min_n
    flag_ec_possibly <- if (length(prob_possibly) > 0) {
      calc_prob(prob_possibly, min_n) > 0.75
    } else {
      FALSE
    }
    # For nonec Amplicons, use cytobands instead of genes to count
    flag_amp <- length(unique(data$band[data$total_cn >= data$background_cn + 4])) >= 1

    class <- if (flag_ec) {
      "circular"
    } else if (flag_ec_possibly) {
      "possibly_circular"
    } else if (flag_amp) {
      "noncircular"
    } else {
      "nofocal"
    }

    data.frame(
      ec_genes = ec_genes,
      ec_cytobands = ec_cytobands,
      ec_cytobands_detail = ec_cytobands_detail,
      ec_possibly_genes = ec_possibly_genes,
      class = class
    )
  }

  sel_cols <- c(
    "sample", "purity", "ploidy", "AScore", "pLOH", "cna_burden",
    paste0("CN", 1:19)
  )
  sel_cols <- sel_cols[sel_cols %in% colnames(data)]
  dt_sample <- data[, summarize_sample(.SD), by = .(sample)]
  dt_sample <- merge(dt_sample, data[match(dt_sample$sample, sample),
    sel_cols,
    with = FALSE
  ], by = "sample", all.x = TRUE)

  lg$info("done")
  list(
    data = data,
    gene = dt_genes,
    sample = dt_sample[order(ec_genes, decreasing = TRUE)]
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
