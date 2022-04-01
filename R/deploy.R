#' Deploy Command Line Interface to System Local Path
#'
#' Only should be used in Unix-like system.
#' For details of the arguments passing to CLI, please check [gcap.workflow()]
#' and [gcap.ASCNworkflow()].
#' @importFrom GetoptLong qq
#' @importFrom utils packageVersion
#'
#' @return Nothing.
#' @export
deploy <- function() {
  if (.Platform$OS.type != "unix") stop("This is designed for Unix-like system.")

  dir <- system.file(package = "gcap", mustWork = TRUE)
  if (packageVersion("GetoptLong") >= "1.1.0") {
    cmd <- qq("ln -sf @{dir}/gcap.R  /usr/local/bin/gcap")
  } else {
    cmd <- NULL
  }
  cmd1 <- qq("ln -sf @{dir}/gcap-bam.R  /usr/local/bin/gcap-bam.R")
  cmd2 <- qq("ln -sf @{dir}/gcap-ascn.R  /usr/local/bin/gcap-ascn.R")

  message("Linking gcap-bam.R command")
  system(cmd1)
  message("Linking gcap-ascn.R command")
  system(cmd2)

  if (!is.null(cmd)) {
    message("Linking main command gcap")
    system(cmd)
  }

  message("Done")
  message("Now you shall run gcap-bam.R and gcap-ascn.R from anywhere.")
  if (!is.null(cmd)) {
    message("Restart your terminal and type 'gcap' command to see how it works")
  }
}
