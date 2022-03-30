#' GCAP workflow for gene-level amplicon prediction from ASCN input
#'
#' Unlike [gcap.workflow], this function directly uses the allele-specific
#' copy number data along with some extra sample information to infer
#' ecDNA genes.
#'
#' @inheritParams gcap.runASCAT
#' @inheritParams gcap.extractFeatures
#' @inheritParams gcap.collapse2Genes
#' @inheritParams gcap.runPrediction
#' @inheritParams gcap.workflow
#' @param data a `data.frame` with following columns. The key columns can be obtained
#' from common allele specific CNV calling software, e.g., ASCAT, Sequenza, FACETS.
#' - chromosome: chromosome names starts with 'chr'.
#' - start: start position of the segment.
#' - end: end position of the segment.
#' - total_cn: total integer copy number of the segment.
#' - minor_cn: minor allele integer copy number of the segment. Set it
#' to `NA` if you don't have this data.
#' - sample: sample identifier.
#' - purity: tumor purity of the sample. Set to `1` if you don't know.
#' - ploidy (optinal): ploidy value of the sample tumor genome.
#' - age (optional): age of the case, use along with `gender`.
#' - gender (optional): gender of the case, use along with `age`.
#' - type (optional): cancer type of the case, use along with `age` and `gender`.
#' Please refer to [gcap.collapse2Genes] to see the supported cancer types.
#' This info is only used in 'XGB56' model. If you don't use this model, you
#' don't need to set it.
#' @return a list of invisible `data.table` and corresponding files saved to local machine.
#' @export
#' @examples
#' data("ascn")
#' data <- ascn
#' rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#' data$purity <- 1
#' rv2 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#' data$age <- 60
#' data$gender <- "XY"
#' rv3 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB32")
#' # If you want to use 'XGB56', you should include 'type' column
#' data$type <- "LUAD"
#' rv4 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB56")
#' # If you only have total integer copy number
#' data$minor_cn <- NA
#' rv5 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'
#' # R6 class fCNA --------------------------------
#' print(rv)
#' print(rv$data)
#' print(rv$sample_summary)
#' print(rv$gene_summary)
#' print(rv$cytoband_summary)
#' @testexamples
#' expect_equal(rv, rv2)
#' expect_equal(length(rv3), 2L)
#' expect_equal(length(rv4), 2L)
#' expect_equal(length(rv5), 2L)
#' expect_error(gcap.ASCNworkflow(data, outdir = tempdir()))
gcap.ASCNworkflow <- function(data,
                              genome_build = c("hg38", "hg19"),
                              model = "XGB32",
                              tightness = 1L,
                              outdir = getwd(),
                              result_file_prefix = paste0("gcap_", uuid::UUIDgenerate(TRUE))) {
  genome_build <- match.arg(genome_build)
  check_model(model)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  if (data.table::is.data.table(data)) data.table::setkey(data, NULL)

  lg <- set_logger()
  lg$info("=========================")
  lg$info("   GCAP ASCN WORKFLOW    ")
  lg$info("=========================")
  lg$info()

  lg$info("============================================================")
  lg$info("Step 1: Extract features and collapse features to gene level")
  lg$info("============================================================")

  model_input <- gcap.runASCATBuildflow(
    data,
    genome_build = genome_build
  )

  lg$info("=======================")
  lg$info("Step 2: Run prediction")
  lg$info("=======================")
  model_input$prob <- gcap.runPrediction(
    model_input,
    model = model
  )

  lg$info("====================================")
  lg$info("Step 3: Run scoring and summarizing")
  lg$info("====================================")
  out <- gcap.runScoring(model_input, genome_build, tightness = tightness)

  save_file <- file.path(outdir, paste0(result_file_prefix, "_prediction_result.rds"))
  lg$info("Saving raw prediction result to {save_file}")
  saveRDS(out$data, file = save_file)
  fCNA <- out$fCNA
  rm(out)
  invisible(gc())

  save_file <- file.path(outdir, paste0(result_file_prefix, c("_fCNA_records.csv", "_sample_info.csv")))
  lg$info("Saving fCNA records and sample info to {paste(save_file, collapse = ', ')}")
  fCNA$saveToFiles(outdir, result_file_prefix)

  lg$info("===========================================")
  lg$info(" Done! Thanks for using GCAP ASCN workflow")
  lg$info("===========================================")

  invisible(fCNA)
}

# runASCATBuildflow -----------------------------------

#' Build data for prediction from absolute copy number data
#'
#' This is is a wrapper of [gcap.extractFeatures()]
#' and [gcap.collapse2Genes()] to combine the feature extraction
#' and predict input generate procedure.
#' If you want to modify the result of [gcap.extractFeatures()],
#' you should always use the two functions instead of this
#' wrapper.
#' @inheritParams gcap.ASCNworkflow
#'
#' @seealso [gcap.runBuildflow]
#'
#' @return a `data.table`.
#' @export
gcap.runASCATBuildflow <- function(data,
                                   genome_build = c("hg38", "hg19")) {
  genome_build <- match.arg(genome_build)
  lg <- set_logger()

  lg$info("extracting sample-level and region-level features")
  fts <- gcap.extractFeatures(
    ascn_data = data,
    genome_build = genome_build
  )

  lg$info("checking if input age, gender (and type) columns")
  if (all(c("age", "gender", "type") %in% colnames(data))) {
    extra_info <- unique(data[, c("sample", "age", "gender", "type")])
  } else if (all(c("age", "gender") %in% colnames(data))) {
    extra_info <- unique(data[, c("sample", "age", "gender")])
  } else {
    lg$info("None exist, set to NULL.")
    extra_info <- NULL
  }

  lg$info("collapsing all data into gene-level prediction input")
  gcap.collapse2Genes(
    fts = fts, extra_info = extra_info,
    include_type = if ("type" %in% colnames(extra_info)) TRUE else FALSE,
    genome_build = genome_build
  )
}
