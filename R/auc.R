#' Get AUC value
#'
#' @param y_pred y prediction vector.
#' @param y y true label vector.
#' @param type AUC type, either 'pr' or 'roc'.
#' @param curve if `TRUE`, generate plot data, the result can be
#' plotted by `plot()`.
#'
#' @return A object.
#' @export
#'
#' @examples
#' if (require("PRROC")) {
#'   set.seed(2021)
#'   auc <- get_auc(sample(1:10, 10), c(rep(0, 5), rep(1, 5)))
#'   auc
#' }
#' @testexamples
#' if (require("PRROC")) {
#'   expect_equal(length(auc), 3)
#'   expect_s3_class(auc, "PRROC")
#' }
get_auc <- function(y_pred, y, type = c("pr", "roc"), curve = FALSE) {
  .check_install("PRROC")
  type <- match.arg(type)

  f <- switch (type,
    pr = PRROC::pr.curve,
    roc = PRROC::roc.curve
  )
  f(y_pred, weights.class0 = y, curve = curve)
}
