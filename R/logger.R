.app <- rappdirs::app_dir("gcap", "ShixiangWang")
.log <- file.path(.app$log(), "gcap.log")

#' @import lgr
#'
.fmt <- "<{logger}> {timestamp} {level_name} [{caller}]: {msg}"

set_logger <- function() {
  if (!dir.exists(.app$log())) dir.create(.app$log(), recursive = TRUE)
  lg <- get_logger_glue("gcap")
  lg$set_appenders(list(
    console = AppenderConsole$new(threshold = "info", layout = LayoutGlue$new(.fmt)),
    fileCon = AppenderFile$new(.log, threshold = "info", layout = LayoutGlue$new(.fmt))
  ))$set_propagate(FALSE)
  return(lg)
}

get_log_file <- function() {
  .log
}

cat_log_file <- function() {
  if (!file.exists(.log)) {
    message("No log available.")
    invisible()
  } else {
    cat(readLines(.log), sep = "\n")
  }
}
