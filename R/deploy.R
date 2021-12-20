#' Deploy Command Line Interface to System Local Path
#'
#' Only should be used in Unix-like system.
#' For details of the arguments passing to CLI, please check [gcap.workflow()]
#' and [gcap.ASCNworkflow()].
#'
#' @return Nothing.
#' @export
deploy <- function() {
  if (.Platform$OS.type != "unix") stop("This is designed for Unix-like system.")

  dir = system.file(package = "gcap")
  cmd1 = qq("ln -sf @{dir}/gcap-wes.R  /usr/local/bin/gcap-wes.R")
  cmd2 = qq("ln -sf @{dir}/gcap-ascn.R  /usr/local/bin/gcap-ascn.R")

  message("Linking gcap-wes.R command")
  system(cmd1)
  message("Linking gcap-ascn.R command")
  system(cmd2)
  message("Done")

  message("Now you shall run gcap-wes.R and gcap-ascn.R from anywhere.")
}

