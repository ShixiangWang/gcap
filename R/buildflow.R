#' Build data for prediction from ASCAT result files
#'
#' This is is a wrapper of [gcap.extractFeatures()]
#' and [gcap.collapse2Genes()] to combine the feature extraction
#' and predict input generate procedure.
#' If you want to modify the result of [gcap.extractFeatures()],
#' you should always use the two functions instead of this
#' wrapper.
#'
#' @inheritParams gcap.extractFeatures
#' @inheritParams gcap.collapse2Genes
#'
#' @return a `data.table`.
#' @export
gcap.runBuildflow <- function(ascat_files,
                              extra_info,
                              include_type = FALSE,
                              genome_build = c("hg38", "hg19"),
                              overlap = 1) {
  genome_build <- match.arg(genome_build)

  lg <- set_logger()

  lg$info("extracting sample-level and region-level features")
  fts <- gcap.extractFeatures(
    ascat_files = ascat_files,
    genome_build = genome_build
  )

  lg$info("collapsing all data into gene-level prediction input")
  gcap.collapse2Genes(
    fts = fts, extra_info = extra_info,
    include_type = include_type, genome_build = genome_build,
    overlap = overlap
  )
}
