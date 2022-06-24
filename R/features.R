#' Extract sample and region level features
#'
#' @param ascat_files a list of file path.
#' Typically the result of [gcap.runASCAT()]
#' @param genome_build genome build version, should be
#' one of 'hg38', 'hg19'.
#' @param ascn_data if `ascat_files` is missing, an alternative
#' `data.frame` can be provided for ASCN data along with
#' purity and ploidy (optional).
#' @import data.table
#' @importFrom sigminer read_copynumber_ascat
#' read_copynumber get_Aneuploidy_score get_cn_ploidy
#' get_pLOH_score sig_tally sig_fit
#'
#' @return a `list`.
#' @export
gcap.extractFeatures <- function(ascat_files,
                                 genome_build = c("hg38", "hg19"),
                                 ascn_data = NULL) {
  genome_build <- match.arg(genome_build)

  lg <- set_logger()

  if (missing(ascat_files)) {
    lg$info("> Extract features from input ascn_data<")
    rvlist <- list()
    ascn_data <- as.data.table(ascn_data)
    rvlist$data <- ascn_data[, c("chromosome", "start", "end", "total_cn", "minor_cn", "sample")]
    if (!"purity" %in% colnames(ascn_data)) ascn_data$purity <- 1
    if (!"ploidy" %in% colnames(ascn_data)) ascn_data$ploidy <- NA_real_
    purity_ploidy <- unique(ascn_data[, c("sample", "purity", "ploidy")])
  } else {
    lg$info("> Extract features from ASCAT results <")
    lg$info()

    lg$info("reading ASCAT file list")
    rvlist <- read_copynumber_ascat(ascat_files)

    lg$info("using unique IDs from file names for avoid the sample name repetition")
    lg$info("back up default sample column to old_sample")
    rvlist$data <- as.data.table(rvlist$data)
    rvlist$data[, old_sample := sample]
    rvlist$data[, sample := sub(".ASCAT.rds", "", source)]
    ids <- unique(rvlist$data$sample)
    names(rvlist$purity) <- names(rvlist$ploidy) <- ids
    rvlist$data$source <- NULL

    lg$info("combining purity and ploidy info as data.frame")
    purity_ploidy <- data.table(
      sample = ids,
      purity = as.numeric(rvlist$purity),
      ploidy = as.numeric(rvlist$ploidy)
    )
  }


  lg$info("generating CopyNumber object in sigminer package")
  cn <- data.table::copy(rvlist$data)
  if ("old_sample" %in% colnames(cn)) cn$old_sample <- NULL
  add_loh <- if (all(is.na(cn$minor_cn))) FALSE else TRUE
  cn <- read_copynumber(
    cn,
    seg_cols = colnames(cn)[1:4],
    max_copynumber = 10000L,
    join_adj_seg = FALSE,
    genome_build = genome_build,
    add_loh = add_loh,
    loh_min_len = 1e3
  )

  lg$info("estimating ploidy from copy number data")
  df_ploidy <- get_cn_ploidy(cn)
  colnames(df_ploidy)[2] <- "ploidy2"
  lg$info("checking if input data contains ploidy and if there are NAs should be overwritten")
  purity_ploidy <- merge(purity_ploidy, df_ploidy, by = "sample", all.x = TRUE)
  purity_ploidy[, ploidy := ifelse(is.na(ploidy), ploidy2, ploidy)]
  purity_ploidy$ploidy2 <- NULL

  lg$info("getting Aneuploidy score")
  # What if ploidy call in ASCAT fail?
  df_aneuploidy <- get_Aneuploidy_score(
    cn,
    ploidy_df = purity_ploidy[, .(sample, ploidy)],
    rm_black_arms = TRUE
  )[, .(sample, AScore)]

  if (add_loh) {
    lg$info("getting pLOH score")
    df_pLOH <- get_pLOH_score(cn)
  } else {
    lg$info("pLOH score is skipped as no LOH data available")
    df_pLOH <- NULL
  }

  lg$info("getting CNA burden")
  df_cna <- cn@summary.per.sample[, .(sample, cna_burden)]

  if (add_loh) {
    lg$info("generating copy number catalog matrix for fitting signature activity")
    cn_tally <- suppressMessages(sig_tally(cn, method = "S"))
    lg$info("fitting copy number signature activity")
    df_act <- suppressMessages(
      sig_fit(
        t(cn_tally$all_matrices$CN_48),
        sig_db = "CNS_TCGA",
        sig_index = "ALL",
        return_class = "data.table"
      )
    )
  } else {
    lg$info("copy number signature activity quantification is skipped as no LOH data available")
    df_act <- NULL
  }

  lg$info("merging data")
  fts_sample <- mergeDTs(list(
    purity_ploidy, df_aneuploidy,
    df_pLOH, df_cna,
    df_act
  ), by = "sample")

  lg$info("feature extraction done")
  lg$info("now you can modify the result and append 'age' and 'gender' columns to the 'fts_sample' element of result list")
  list(
    fts_sample = fts_sample,
    fts_region = rvlist$data
  )
}
