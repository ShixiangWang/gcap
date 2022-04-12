#' Deploy Command Line Interface to System Local Path
#'
#' Only should be used in Unix-like system.
#' For details of the arguments passing to CLI, please check [gcap.workflow()]
#' and [gcap.ASCNworkflow()].
#' @importFrom GetoptLong qq
#' @importFrom utils packageVersion
#' @param target the target path to deploy the CLI.
#' @return Nothing.
#' @export
deploy <- function(target = "/usr/local/bin") {
  stopifnot(length(target) == 1L, dir.exists(target))
  if (file.access(target, 2) != 0) {
    stop("Cannot link to a target path without write permission")
  }

  if (.Platform$OS.type != "unix") stop("This is designed for Unix-like system.")

  dir <- system.file(package = "gcap", mustWork = TRUE)
  if (packageVersion("GetoptLong") >= "1.1.0") {
    cmd <- qq("ln -sf @{dir}/gcap.R @{target}/gcap")
  } else {
    cmd <- NULL
  }
  cmd1 <- qq("ln -sf @{dir}/gcap-bam.R  @{target}/gcap-bam.R")
  cmd2 <- qq("ln -sf @{dir}/gcap-ascn.R  @{target}/gcap-ascn.R")

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
