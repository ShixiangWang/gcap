.check_install <- function(pkg, bioc = FALSE, ...) {
  lg <- set_logger()
  install_func <- if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) install_func(pkg, ...)
  lg$info("Required package ", pkg, " has been installed.")
}
