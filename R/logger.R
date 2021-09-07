#' @import lgr
#'
.fmt <- "<{logger}> {timestamp} {level_name} [{caller}]: {msg}"

set_logger <- function() {
  lg <- get_logger_glue("gcap")
  lg$set_appenders(list(
    console = AppenderConsole$new(threshold = "info", layout = LayoutGlue$new(.fmt))
    #fileCon = AppenderFile$new(get_log_file(), threshold = "info", layout = LayoutGlue$new(.fmt))
  ))$set_propagate(FALSE)
  return(lg)
}

