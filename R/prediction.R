#' Run gene-level amplicon prediction
#'
#' @param data data to predict (`data.frame`/`matrix` format), from
#' [gcap.collapse2Genes()] in general.
#' @param feature one of 'with_type' and 'without_type' to select
#' model whether includes cancer type or not.
#' @param target one of 'circle' and 'nonLinear' to select model
#' if predict circle amplicon or non-linear amplicon.
#' @param use_toy if `TRUE`, use a toy model for showing example,
#' not for practical use.
#' @param custom_model a outer custom xgboost model for prediction.
#'
#' @return a numeric vector representing prob.
#' @importFrom stats predict
#' @export
#'
#' @examples
#' data("ec")
#' y_pred <- gcap.runPrediction(ec)
#' y_pred
#' @testexamples
#' expect_equal(length(y_pred), 2020L)
gcap.runPrediction <- function(data,
                               feature = c("without_type", "with_type"),
                               target = c("circle", "nonLinear"),
                               use_toy = FALSE,
                               custom_model = NULL) {
  feature <- match.arg(feature)
  target <- match.arg(target)

  lg <- set_logger()

  lg$info("reading trained model based on input feature and target setting")
  lg$info("feature: {feature}; target: {target}")
  if (use_toy) {
    model <- readRDS(
      system.file(
        "extdata", "toy_model.rds",
        package = "gcap", mustWork = TRUE
      )
    )
  } else {
    if (!is.null(custom_model)) {
      lg$info("a custom model is selected from user input")
      model <- custom_model
    } else {
      if (feature == "with_type" && target == "circle") {
        model <- readRDS(
          system.file(
            "extdata", "toy_model.rds",
            package = "gcap", mustWork = TRUE
          )
        )
      } else if (feature == "without_type" && target == "circle") {
        model <- readRDS(
          system.file(
            "extdata", "final_xgb_model_circle_without_type.rds",
            package = "gcap", mustWork = TRUE
          )
        )
      } else if (feature == "with_type" && target == "nonLinear") {
        model <- readRDS(
          system.file(
            "extdata", "toy_model.rds",
            package = "gcap", mustWork = TRUE
          )
        )
      } else {
        model <- readRDS(
          system.file(
            "extdata", "final_xgb_model_nonLinear_without_type.rds",
            package = "gcap", mustWork = TRUE
          )
        )
      }
    }
  }

  lg$info("selecting necessary features from input data")
  data <- as.matrix(data[, model$feature_names, with = FALSE])

  lg$info("running prediction")
  if ("best_ntreelimit" %in% names(model)) {
    predict(model, data, ntreelimit = model$best_ntreelimit)
  } else {
    predict(model, data)
  }
}
