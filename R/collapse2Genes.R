#' Generate unified gene-level feature data
#'
#' @importFrom sigminer get_intersect_size
#' @inheritParams gcap.extractFeatures
#' @param fts (modified) result from [gcap.extractFeatures()]
#' @param extra_info (optional) a `data.frame` with 3 columns 'sample',
#' 'age' and 'gender', for including cancer type, check
#' parameter `include_type`. For gender, should be 'XX' or 'XY',
#' also could be `0` for 'XX' and `1` for 'XY'.
#' @param include_type if `TRUE`, a fourth column named 'type'
#' should be included in `extra_info`, the supported cancer
#' type should be described with [TCGA cancer type abbr.](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations).
#' @param fix_type default is `TRUE`, only cancer types used in pre-trained
#' models are used, others will be convert to `NA`. If `FALSE`, only generating
#' one-hot encoding for cancer types in input data.
#'
#' @return a `data.table`.
#' @export
gcap.collapse2Genes <- function(fts,
                                extra_info = NULL,
                                include_type = FALSE,
                                fix_type = TRUE,
                                genome_build = c("hg38", "hg19")) {
  stopifnot(
    is.list(fts),
    all.equal(names(fts), c("fts_sample", "fts_region"))
  )

  genome_build <- match.arg(genome_build)
  lg <- set_logger()
  lg$info("please make sure the first 3 columns of `fts$fts_region` are for chr, start, end.")

  if (!is.null(extra_info)) {
    stopifnot(
      is.data.frame(extra_info),
      if (include_type) {
        all.equal(colnames(extra_info), c("sample", "age", "gender", "type"))
      } else {
        all.equal(colnames(extra_info), c("sample", "age", "gender"))
      }
    )

    extra_info <- as.data.table(extra_info)
    if (include_type) {
      lg$info("one-hot encoding cancer type")
      if (fix_type) {
        # UPDATE!!!
        types <- c(
          "BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
          "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV",
          "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC", "UVM"
        )
        lg$info("valid cancer types: {paste(types, collapse=',')}")
      } else {
        types <- sort(unique(extra_info$type))
      }
      extra_info$type <- factor(
        as.character(extra_info$type),
        levels = types
      )
      extra_info <- mltools::one_hot(extra_info, cols = "type")
    }

    lg$info("checking data type for age column")
    if (!is.numeric(extra_info$age)) {
      lg$warn("non numeric type found, transforming it")
      extra_info[, age := as.numeric(age)]
    }
    lg$info("checking data type for gender column")
    if (!is.numeric(extra_info$gender)) {
      lg$warn("non numeric type found, transforming XY to 1, otherwise 0")
      extra_info[, gender := data.table::fifelse(
        gender == "XY",
        1L, 0L,
        na = 0L
      )]
    }
  }

  lg$info("collapsing region-level features to gene-level")
  dt <- collapse_to_genes(fts$fts_region, genome_build)

  lg$info("merging gene-level and sample-level data")
  dt <- merge(dt, fts$fts_sample, by = "sample", all.x = TRUE)

  if (!is.null(extra_info)) {
    lg$info("merging data and extra info")
    dt <- merge(dt, extra_info, by = "sample", all.x = TRUE)
  }

  lg$info("merging data and prior amplicon frequency data")
  amp_freq <- readRDS(
    system.file(
      "extdata", "amplicon_freq.rds",
      package = "gcap", mustWork = TRUE
    )
  )
  dt <- merge(dt, amp_freq, by = "gene_id", all.x = TRUE)
  dt[, `:=`(
    freq_Linear = ifelse(is.na(freq_Linear), 0, freq_Linear),
    freq_BFB = ifelse(is.na(freq_BFB), 0, freq_BFB),
    freq_Circular = ifelse(is.na(freq_Circular), 0, freq_Circular),
    freq_HR = ifelse(is.na(freq_HR), 0, freq_HR)
  )]

  lg$info("done")
  dt
}

collapse_to_genes <- function(x, genome_build = "hg38", overlap = 1, drop = TRUE) {
  lg <- set_logger()

  x <- data.table::as.data.table(x)
  lg$info("checking input chromosome names")
  if (!grepl("chr", x[[1]][1])) {
    x[[1]] <- paste0("chr", x[[1]])
  }

  ref_file <- system.file(
    "extdata", paste0(genome_build, "_target_genes.rds"),
    package = "gcap", mustWork = TRUE
  )
  lg$info("reading reference file {ref_file}")
  y <- readRDS(ref_file)
  colnames(x)[1:3] <- colnames(y)[1:3] <- c("chr", "start", "end")

  lg$info("finding overlaps")
  data.table::setkey(y, chr, start, end)
  out <- data.table::foverlaps(x, y)[!is.na(gene_id)]

  lg$info("calculating intersect size")
  out[, intersect_size := get_intersect_size(i.start, i.end, start, end)]
  # Calculate the region cov ratio
  out[, intersect_ratio := intersect_size / (abs(end - start) + 1)]

  if (drop) {
    lg$info("keeping records with >= {round(100 * overlap, 3)}% overlap ratio with a gene")
    out <- out[
      intersect_ratio >= overlap,
      .(
        sample, gene_id,
        total_cn, minor_cn
      )
    ]
  }
  out
}
