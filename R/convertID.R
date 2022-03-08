#' Convert Gene IDs between Ensembl and Hugo Symbol System
#'
#' @param IDs a character vector to convert.
#' @param type type of input `IDs`, could be 'ensembl' or 'symbol'.
#' @param genome_build reference genome build.
#' @param multiple 	if `TRUE`, return a `data.table` instead of a string vector,
#' so multiple identifier mappings can be kept.
#'
#' @return a vector or a `data.table`.
#' @export
#'
#' @examples
#' \donttest{
#' convertID("ENSG00000243485")
#' convertID("ENSG00000243485", multiple = TRUE)
#' convertID(c("TP53", "KRAS", "EGFR", "MYC"), type = "symbol")
#' }
convertID <- function(IDs, type = c("ensembl", "symbol"),
                      genome_build = c("hg38", "hg19", "mm10", "mm9"),
                      multiple = FALSE) {
  type <- match.arg(type)
  genome_build <- match.arg(genome_build)
  ref_data <- get_ref_data(genome_build)
  if (genome_build %in% c("hg38", "hg19")) {
    ref_data$gene_id <- substr(ref_data$gene_id, 1, 15)
  } else {
    ref_data$gene_id <- substr(ref_data$gene_id, 1, 18)
  }

  if (type == "symbol") {
    from <- "gene_name"
    to <- "gene_id"
  } else {
    from <- "gene_id"
    to <- "gene_name"
  }

  IDConverter::convert_custom(IDs, from, to, ref_data, multiple)
}
