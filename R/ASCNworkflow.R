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
#' @param data a `data.frame` with following rows. The key columns can be obtained
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
#' # If you only have total integer copy number
#' data$minor_cn <- NA
#' rv4 = gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#' @testexamples
#' expect_equal(rv, rv2)
#' expect_equal(ncol(rv3$by_gene), 35L)
#' expect_equal(ncol(rv4$by_gene), 15L)
#' expect_error(gcap.ASCNworkflow(data, outdir = tempdir()))
gcap.ASCNworkflow <- function(data,
                              genome_build = c("hg38", "hg19"),
                              model = "XGB32",
                              target = "circle",
                              outdir = getwd(),
                              result_file_prefix = paste0("gcap_", uuid::UUIDgenerate(TRUE))) {
  genome_build <- match.arg(genome_build)
  target <- match.arg(target, choices = c("circle", "nonLinear"), several.ok = TRUE)
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
  for (t in target) {
    model_input[[paste0("pred_", t)]] <- gcap.runPrediction(
      model_input,
      model = model,
      target = t
    )
  }
  save_file <- file.path(outdir, paste0(result_file_prefix, "_by_gene.csv"))
  lg$info("Saving result to {save_file}")
  data.table::fwrite(model_input, file = save_file)

  lg$info("=======================")
  lg$info("Step 3: Run scoring")
  lg$info("=======================")
  out <- gcap.runScoring(model_input, genome_build)

  save_file <- file.path(outdir, paste0(result_file_prefix, "_by_case.csv"))
  lg$info("Saving result to {save_file}")
  data.table::fwrite(out, file = save_file)

  lg$info("===========================================")
  lg$info(" Done! Thanks for using GCAP ASCN workflow")
  lg$info("===========================================")

  invisible(list(
    by_gene = model_input,
    by_case = out
  ))
}

# runASCATBuildflow -----------------------------------

gcap.runASCATBuildflow <- function(data,
                                   genome_build = c("hg38", "hg19")) {
  genome_build <- match.arg(genome_build)
  lg <- set_logger()

  lg$info("extracting sample-level and region-level features")
  fts <- gcap.extractFeatures(
    ascn_data = data,
    genome_build = genome_build
  )

  lg$info("checking if input age and gender columns")
  if (all(c("age", "gender") %in% colnames(data))) {
    extra_info <- unique(data[, c("sample", "age", "gender")])
  } else {
    lg$info("None exist, set to NULL.")
    extra_info <- NULL
  }

  lg$info("collapsing all data into gene-level prediction input")
  gcap.collapse2Genes(
    fts = fts, extra_info = extra_info,
    genome_build = genome_build
  )
}
