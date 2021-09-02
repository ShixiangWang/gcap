# https://mlr3learners.mlr-org.com/reference/mlr_learners_classif.xgboost.html

#' Auto tuning for xgboost binary prediction
#'
#' @param data a `data.frame`
#' @param target target column name
#' @param positive positive case representation, like "1", "+", etc.
#' @param seed random seed
#'
#' @return a `list`
#' @export
#' @importFrom mlr3tuning AutoTuner
train <- function(data, target, positive, seed = 2021L) {

  stopifnot(length(positive) == 1L, is.data.frame(data))
  .check_install("mlr3")
  .check_install("mlr3learners")
  .check_install("mlr3tuning")

  # Task
  set.seed(seed)
  lg <- set_logger()

  data <- data.table::as.data.table(data)
  if (!is.factor(data[[target]])) {
    lg$info("transform column {target} to factor, note the out levels")
    data[[target]] <- as.factor(data[[target]])
    lg$info("  levels: ", paste(levels(data[[target]]), collapse = ", "))
  }

  lg$info("create task")
  task = mlr3::as_task_classif(data, target = target, id = "auto-feature-selector", positive = as.character(positive))
  rm(data)
  invisible(gc())

  # Inner CV and hyper-parameter tuning
  learner = mlr3::lrn("classif.xgboost", predict_type = "prob")
  learner$param_set$values$nrounds = 200
  learner$param_set$values$eval_metric = "aucpr"
  learner$param_set$values$early_stopping_rounds = 5

  resampling = mlr3::rsmp("holdout", ratio = 0.9)
  measure = mlr3::msr("classif.prauc")
  search_space = paradox::ps(
     = paradox::p_dbl(lower = 0.001, upper = 0.1))
  terminator = mlr3tuning::trm("evals", n_evals = 5)
  tuner = mlr3tuning::tnr("grid_search", resolution = 10)

  at = AutoTuner$new(learner, resampling, measure, terminator, tuner, search_space)

  # Outer CV
 outer_resampling = rsmp("cv", folds = 10)
 rr = resample(task, at, outer_resampling, store_models = TRUE)

 # Measures
 extract_inner_tuning_results(rr)

 rr$score()
 # Unbiased estimation for measures
 rr$aggregate()

 # Final model
 at$train(task)

}


#' Auto-select features for xgboost binary prediction
#'
#' @param data a `data.frame`
#' @param target target column name
#' @param positive positive case representation, like "1", "+", etc.
#' @param seed random seed
#'
#' @return a `list`
#' @export
#' @importFrom mlr3fselect AutoFSelector
fselect <- function(data, target, positive, seed = 2021L) {
  stopifnot(length(positive) == 1L, is.data.frame(data))
  .check_install("mlr3")
  .check_install("mlr3learners")
  .check_install("mlr3fselect")
  .check_install("mlr3tuning")

  set.seed(seed)
  lg <- set_logger()

  data <- data.table::as.data.table(data)
  if (!is.factor(data[[target]])) {
    lg$info("transform column {target} to factor, note the out levels")
    data[[target]] <- as.factor(data[[target]])
    lg$info("  levels: ", paste(levels(data[[target]]), collapse = ", "))
  }

  lg$info("create task")
  task = mlr3::as_task_classif(data, target = target, id = "auto-feature-selector", positive = as.character(positive))
  rm(data)
  invisible(gc())

  learner = mlr3::lrn("classif.xgboost", predict_type = "prob")
  terminator = mlr3tuning::trm("evals", n_evals = 10)
  fselector = mlr3fselect::fs("random_search")
  resampling = mlr3::rsmp("holdout", ratio = 0.9)
  measure = mlr3::msr("classif.prauc")

  # rewrite parameters
  #learner$param_set$values$nrounds <- paradox::to_tune(c(1, 100))
  learner$param_set$values$verbose <- 1

  at = AutoFSelector$new(
    learner = learner,
    resampling = resampling,
    measure = measure,
    terminator = terminator,
    fselector = fselector
  )

  rr <- at$train(task)
  return(list(
    rr = rr,
    fs = rr$fselect_result$features
  ))

}
