.check_install <- function(pkg, bioc = FALSE, ...) {
  lg <- set_logger()
  install_func <- if (bioc) BiocManager::install else utils::install.packages
  if (bioc) {
    .check_install("BiocManager")
  }
  if (!requireNamespace(pkg)) {
    lg$info("installing required package {pkg}")
    install_func(pkg, ...)
  }
}


mergeDTs <- function(dt_list, by = NULL, sort = FALSE) {
  Reduce(
    function(...) {
      merge(..., by = by, all = TRUE, sort = sort)
    }, dt_list
  )
}

# Global variables
utils::globalVariables(
  c(
    ".", "old_sample", "ploidy", "AScore", "cna_burden",
    "age", "chr", "freq_BFB", "freq_Circular", "freq_HR",
    "freq_Linear", "gender", "gene_id", "i.end", "i.start",
    "intersect_ratio", "intersect_size", "minor_cn", "total_cn"
  )
)
