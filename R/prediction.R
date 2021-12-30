#' Run gene-level amplicon prediction
#'
#' @param data data to predict (`data.frame`/`matrix` format), from
#' [gcap.collapse2Genes()] in general.
#' @param model model name ("XGB11", "XGB32", "XGB54") or a custom model
#' from input. 'toy' can be used for test.
#' @param target one of 'circle' and 'nonLinear' to select model
#' if predict circle amplicon or non-linear amplicon.
#' @param use_best_ntreelimit if `TRUE`, use default `best_ntreelimit`
#' proposed by XGBOOST, otherwise we determine a iteration number (
#' i.e. tree number) with more careful processing to avoid over-fitting.
#'
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
                               model = "XGB32",
                               target = c("circle", "nonLinear"),
                               use_best_ntreelimit = FALSE) {
  stopifnot(data.table::is.data.table(data))
  target <- match.arg(target)

  lg <- set_logger()

  lg$info("reading trained model based on input target setting: {target}")
  if (is.character(model)) {
    if (model == "toy") {
      model <- readRDS(
        system.file(
          "extdata", "toy_model.rds",
          package = "gcap", mustWork = TRUE
        )
      )
    } else {
      if (target == "circle") {
        modfile <- switch(model,
          XGB11 = "xgb_stepwise_model_NF11.rds",
          XGB32 = "xgb_stepwise_model_NF32.rds",
          XGB54 = "xgb_stepwise_model_NF54.rds",
          stop("Unsupported model input!")
        )
      } else {
        stop("Currently not supported!")
      }

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
    predict(model, data, ntreelimit = if (use_best_ntreelimit) {
      model$best_ntreelimit
    } else {
      determine_iter(model)
    })
  } else {
    predict(model, data)
  }
}

determine_iter <- function(m) {
  log <- m$evaluation_log
  nc <- ncol(log)
  if (nc < 2) {
    return(NULL)
  } else if (nc == 2 || nc > 3) {
    return(order(log[[nc]], decreasing = TRUE)[1])
  } else {
    z <- log[rev(seq_len(.N))][eval_aucpr > 0.5][order(abs(eval_aucpr - train_aucpr), -eval_aucpr, decreasing = FALSE)]$iter[1]
    if (is.na(z)) NULL else z
  }
}
