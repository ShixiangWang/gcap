#' @import lgr
.tf <- tempfile(fileext = ".gcap.log.txt")
.fmt <- "<{logger}> {timestamp} {level_name} [{caller}]: {msg}"

set_logger <- function() {
  lg <- get_logger_glue("gcap")
  lg$set_appenders(list(
    console = AppenderConsole$new(threshold = "info", layout = LayoutGlue$new(.fmt)),
    fileCon = AppenderFile$new(.tf, threshold = "info", layout = LayoutGlue$new(.fmt))
  ))$set_propagate(FALSE)
  return(lg)
}


#' Get log file
#'
#' @return a file path
#' @export
#'
#' @examples
#' f <- get_log_file()
#' f
#' @testexamples
#' expect_equal(length(f), 1)
get_log_file <- function() {
  .tf
}
