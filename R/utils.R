.check_install <- function(pkg, bioc = FALSE, ...) {
  lg <- set_logger()
  install_func <- if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) install_func(pkg, ...)
  lg$info("Required package ", pkg, " has been installed.")
}


mergeDTs <- function(dt_list, by = NULL, sort = FALSE) {
  Reduce(
    function(...) {
      merge(..., by = by, all = TRUE, sort = sort)
    }, dt_list)
}

# Global variables
utils::globalVariables(
  c(".", "old_sample", "ploidy", "AScore", "cna_burden")
)
