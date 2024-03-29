#' Run gene-level circular prediction
#'
#' @param data data to predict (`data.frame`/`matrix` format), from
#' [gcap.collapse2Genes()] in general.
#' @param model model name ("XGB11", "XGB32", "XGB56") or a custom model
#' from input. 'toy' can be used for test.
#' @return a numeric vector representing prob.
#' @importFrom stats predict
#' @export
#'
#' @examples
#' data("ec")
#' # Use toy model for illustration
#' y_pred <- gcap.runPrediction(ec, "toy")
#' y_pred
#' @testexamples
#' expect_equal(length(y_pred), 2020L)
gcap.runPrediction <- function(data,
                               model = "XGB11") {
  stopifnot(data.table::is.data.table(data))
  if (utils::packageVersion("xgboost") >= "1.6") {
    warning("The gcap model was developed with xgboost <1.6, it is recommended to install the previous version with command:", immediate. = TRUE)
    message('install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz", repos = NULL, type = "source")')
  }

  lg <- set_logger()

  if (is.character(model)) {
    if (model == "toy") {
      model <- readRDS(
        system.file(
          "extdata", "toy_model.rds",
          package = "gcap", mustWork = TRUE
        )
      )
    } else {
      modfile <- switch(model,
        XGB11 = "XGB_NF11.rds",
        XGB32 = "XGB_NF32.rds",
        XGB56 = "XGB_NF56.rds",
        stop("Unsupported model input!")
      )

      lg$info("using model file {modfile}")
      model <- readRDS(
        system.file(
          "extdata", modfile,
          package = "gcap", mustWork = TRUE
        )
      )
    }
  } else {
    lg$info("a custom model is selected from user input")
  }

  lg$info("selecting necessary features from input data")
  if (is.null(data$pLOH)) data$pLOH <- NA
  if (is.null(data$age)) data$age <- NA
  if (is.null(data$gender)) data$gender <- NA
  if (all(is.na(data$minor_cn))) {
    # Set copy number signatures to NA
    z <- paste0("CN", 1:19)
    data[, (z) := rep(list(NA), 19)]
  }
  data <- tryCatch(
    as.matrix(data[, model$feature_names, with = FALSE]),
    error = function(e) {
      lg$fatal(e$message)
      stop("Please try a model with less features if you don't get so many features!")
    }
  )

  lg$info("running prediction")
  if ("best_ntreelimit" %in% names(model)) {
    predict(model, data, ntreelimit = model$best_ntreelimit)
  } else {
    predict(model, data)
  }
}
