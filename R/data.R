## how to document datasets: you need to specify @docType and @name; do not
## forget NULL in the end

#' Example ecDNA training data
#' @docType data
#' @name ec
#' @format A `data.table`
#' @source Generate from `data-raw/`
#' @examples
#' data("ec")
NULL

#' Example allele specific copy number (ASCN) data
#' @docType data
#' @name ascn
#' @format A `data.frame`
#' @source Generate from `data-raw/`, raw source from our study by
#' calling ASCAT v3.0 alpha on corresponding WES sequencing data.
#' @examples
#' data("ascn")
NULL
