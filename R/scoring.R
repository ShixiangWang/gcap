#' Summarize prediction result into gene/sample-level
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' The column storing prediction result must be "prob".
#' @param prob_cutoff a numeric cutoff `[0, 1]` to determine the positive class.
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
#' expect_equal(nrow(score), 1L)
gcap.runScoring <- function(data,
                            genome_build = "hg38",
                            prob_cutoff = 0.5,
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

  data = data[!is.na(data$prob)]
  data$status = ifelse(data$prob > prob_cutoff, 1L, 0L)
  data$prob = NULL

  blood_cn = readRDS(system.file("extdata", "blood_gene_cn.rds", package = "gcap"))
  colnames(blood_cn)[2] = c("blood_cn_median")

  cytobands = as.data.table(
    sigminer::get_genome_annotation("cytobands", genome_build = genome_build)
    )[chrom %in% paste0("chr", 1:22)]
  cytobands$stain = NULL
  cytobands$start = cytobands$start + 1L

  ref_file <- system.file(
    "extdata", paste0(genome_build, "_target_genes.rds"),
    package = "gcap", mustWork = TRUE
  )
  refgene = readRDS(ref_file)
  gene_cytobands = overlaps(refgene, cytobands)
  gene_cytobands[, intersect_ratio := intersect_size / (i.end - i.start + 1) ]
  gene_cytobands = gene_cytobands[, .SD[which.max(intersect_ratio)], by = .(gene_id)][order(chr, i.start)]
  #gene_cytobands[, gene_id := factor(gene_id, levels = gene_id)]
  gene_cytobands = gene_cytobands[, .(gene_id, band = paste(chr, band, sep = ":"))]

  data = merge(data, blood_cn, by = "gene_id", all.x = TRUE)
  data$blood_cn_median = ifelse(is.na(data$blood_cn_median), 2, data$blood_cn_median)
  data = merge(data, gene_cytobands, by = "gene_id", all.x = TRUE)

  cytoband_cn = data[, .(cytoband_cn_median = median(total_cn, na.rm = TRUE)), by = .(sample, band)]
  data = merge(data, cytoband_cn, by = c("sample", "band"), all.x = TRUE)

  # sample_cn = data[, .(sample_cn_outlier = median(total_cn, na.rm = TRUE) + 1.5*IQR(total_cn, na.rm = TRUE)), by = .(sample)]
  # data = merge(data, sample_cn, by = c("sample"), all.x = TRUE)
  data$background_cn = round(data$blood_cn_median * data$ploidy / 2)

  if (apply_filter) {
    # !!! 后面再想想，并探索下实际有没有用
    # status should set to 0 if total_cn <= background_cn + 1
    data$status[data$total_cn <= data$background_cn + 1] = 0L
    # status should set to 0 if total_cn <= cytoband_cn_median - 1
    data$status[data$total_cn <= data$cytoband_cn_median - 1] = 0L
  }

  # Sample level summary

  dt <- gene_result[sample %in% clinfo$Tumor_Sample_Barcode,
                    .(gene_id, sample, status = ifelse(pred_circle > 0.5, 1L, 0L),
                      total_cn)] %>%
    dplyr::group_by(gene_id) %>%
    dplyr::filter(sum(status) != 0) %>%
    dplyr::ungroup()

  dt_genes = dt %>%
    dplyr::group_by(gene_id, sample) %>%
    summarise(status = ifelse(sum(status) > 0, 1L, 0L),
              total_cn = mean(total_cn, na.rm = TRUE), .groups = "drop") %>%
    group_by(gene_id) %>%
    filter(sum(status) >= 1) %>%
    group_by(gene_id, status) %>%
    summarise(total_cn_mean = mean(total_cn, na.rm = TRUE),
              N = sum(status), .groups = "drop") %>%
    dplyr::mutate(status = ifelse(status == 1L, "ecDNA+", "ecDNA-"))

  dt_genes = dt_genes %>%
    select(-N) %>%
    tidyr::pivot_wider(names_from = "status", values_from = "total_cn_mean") %>%
    left_join(dt_genes %>% filter(status == "ecDNA+") %>% select(-total_cn_mean, -status),
              by = "gene_id") %>%
    arrange(desc(N))

  # Sample level summary


}
