#' Run gene-level circular prediction
#'
#' @param data data to predict (`data.frame`/`matrix` format), from
#' [gcap.collapse2Genes()] in general.
#' @param model model name ("XGB11", "XGB32", "XGB54") or a custom model
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
                               model = "XGB32") {
  stopifnot(data.table::is.data.table(data))
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
        XGB11 = "xgb_stepwise_model_NF11.rds",
        XGB32 = "xgb_stepwise_model_NF32.rds",
        XGB54 = "xgb_stepwise_model_NF54.rds",
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
