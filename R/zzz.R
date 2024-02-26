.onAttach <- function(libname, pkgname) {
  version <- packageDescription(pkgname, fields = "Version")

  msg <- paste0(
    pkgname, " version ", version,
    "\n- Project URL at https://github.com/ShixiangWang/gcap",
    "\n\nCitation:\n",
    "\tWang, S., Wu, CY., He, MM. et al. Machine learning-based extrachromosomal DNA identification in \n\tlarge-scale cohorts reveals its clinical implications in cancer. Nat Commun 15, 1515 (2024). https://doi.org/10.1038/s41467-024-45479-6"
  )
  packageStartupMessage(msg)
}


utils::globalVariables(
  c("oncogenes")
)
