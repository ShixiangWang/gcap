#' Summarize prediction result into gene/sample-level
#'
#' @inheritParams gcap.collapse2Genes
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' The column storing prediction result must be "prob".
#' @param prob_cutoff a numeric cutoff `[0, 1]` to determine the positive class.
#' @param N_cutoff an integer cutoff to deterimine sample class. e.g.,
#' `2` means a sample with >=2 predicted ecDNA member genes will be assigned
#' to ecDNA class, otherwise >=2 detected focal amplicon genes will be assigned
#' to focal AMP class.
#' @param apply_filter if `TRUE`, apply some filters to ease false positive calls.
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
#' expect_equal(length(score), 2L)
gcap.runScoring <- function(data,
                            genome_build = "hg38",
                            prob_cutoff = 0.5,
                            N_cutoff = 2,
                            apply_filter = TRUE) {
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

  data <- data[!is.na(data$prob)]
  data$status <- ifelse(data$prob > prob_cutoff, 1L, 0L)
  data$prob <- NULL

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
  # gene_cytobands[, gene_id := factor(gene_id, levels = gene_id)]
  gene_cytobands <- gene_cytobands[, .(gene_id, band = paste(chr, band, sep = ":"))]

  data <- merge(data, blood_cn, by = "gene_id", all.x = TRUE)
  data$blood_cn_median <- ifelse(is.na(data$blood_cn_median), 2, data$blood_cn_median)
  data <- merge(data, gene_cytobands, by = "gene_id", all.x = TRUE)

  cytoband_cn <- data[, .(cytoband_cn_median = median(total_cn, na.rm = TRUE)), by = .(sample, band)]
  data <- merge(data, cytoband_cn, by = c("sample", "band"), all.x = TRUE)

  # sample_cn = data[, .(sample_cn_outlier = median(total_cn, na.rm = TRUE) + 1.5*IQR(total_cn, na.rm = TRUE)), by = .(sample)]
  # data = merge(data, sample_cn, by = c("sample"), all.x = TRUE)
  data$background_cn <- round(data$blood_cn_median * data$ploidy / 2)

  if (apply_filter) {
    # !!! 后面再想想，并探索下实际有没有用
    # status should set to 0 if total_cn <= background_cn + 1
    data$status[data$total_cn <= data$background_cn + 1] <- 0L
    # status should set to 0 if total_cn <= cytoband_cn_median - 1
    data$status[data$total_cn <= data$cytoband_cn_median - 1] <- 0L
  }

  # Gene level summary
  dt <- data[, .SD[sum(status) > 0], by = .(gene_id)][
    , .(
      total_cn = mean(total_cn, na.rm = TRUE),
      minor_cn = mean(minor_cn, na.rm = TRUE),
      background_cn = mean(background_cn, na.rm = TRUE),
      cytoband_cn = mean(cytoband_cn_median, na.rm = TRUE),
      N = sum(!is.na(status))
    ),
    by = .(gene_id, status)
  ]
  dt[, status := factor(ifelse(status == 1L, "pos", "neg"),
    levels = c("neg", "pos")
  )]

  dt_genes <- data.table::dcast(dt, gene_id ~ status,
    value.var = c("total_cn", "minor_cn", "background_cn", "cytoband_cn", "N"),
    drop = FALSE
  )

  # Sample level summary
  # Number of ec Genes, ec Status and sample classification
  summarize_sample <- function(data) {
    ec_genes <- sum(data$status)
    ec_status <- as.integer(ec_genes >= N_cutoff)
    N_AMP <- sum(data$total_cn > data$background_cn + 0.5)
    class <- ifelse(ec_genes >= N_cutoff, "ecDNA_AMP",
      ifelse(N_AMP >= N_cutoff, "focal_AMP", "nonfocal")
    )
    data.frame(
      ec_genes = ec_genes,
      ec_status = ec_status,
      class = class
    )
  }

  dt_sample <- data[, summarize_sample(.SD), by = .(sample)]

  list(
    gene = dt_genes,
    sample = dt_sample
  )
}
