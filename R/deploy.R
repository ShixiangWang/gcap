#' Deploy Command Line Interface to System Local Path
#'
#' Only should be used in Unix-like system.
#' For details of the arguments passing to CLI, please check [gcap.workflow()]
#' and [gcap.ASCNworkflow()].
#' @importFrom GetoptLong qq
#'
#' @return Nothing.
#' @export
deploy <- function() {
  if (.Platform$OS.type != "unix") stop("This is designed for Unix-like system.")

  dir <- system.file(package = "gcap")
  cmd1 <- qq("ln -sf @{dir}/gcap-bam.R  /usr/local/bin/gcap-bam.R")
  cmd2 <- qq("ln -sf @{dir}/gcap-ascn.R  /usr/local/bin/gcap-ascn.R")

  message("Linking gcap-bam.R command")
  system(cmd1)
  message("Linking gcap-ascn.R command")
  system(cmd2)
  message("Done")

  message("Now you shall run gcap-bam.R and gcap-ascn.R from anywhere.")
}
