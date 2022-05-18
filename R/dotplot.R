#' Plot CN Mean vs. Frequency Dotplot
#'
#' @param fCNA a [fCNA] object.
#' @param by the level of focal amplicons.
#' @param filter a filter expression based on `cn` (CN mean) and `N` (frequency).
#' @param include amplicon type to include for plotting.
#' @param ... other parameters passing to `ggrepel::geom_label_repel()`.
#'
#' @return a ggplot.
#' @export
gcap.dotplot <- function(fCNA, by = c("gene_id", "band", "chr"), filter = cn > 50 | (N > 1 & cn > 20), include = c("circular", "possibly_circular"), ...) {
  .check_install("ggrepel")
  .check_install("cowplot")
  data <- data.table::copy(fCNA$data)
  by <- match.arg(by)
  if (by == "chr") {
    data$chr <- gsub("(.*):(.*)", "\\1", data$band)
  }
  genes_summary <- data[
    , .(
      cn = mean(total_cn[amplicon_type %in% include], na.rm = TRUE),
      N = sum(amplicon_type %in% include, na.rm = TRUE)
    ),
    by = by
  ][N > 0]
  genes_summary <- genes_summary[!is.na(genes_summary[[by]])]

  e <- substitute(filter)
  p <- ggplot2::ggplot(data = genes_summary, ggplot2::aes(x = N, y = cn)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.2, col = "black") +
    ggrepel::geom_label_repel(
      ggplot2::aes(label = .data[[by]]),
      data = subset(genes_summary, eval(e)),
      ...
    ) +
    ggplot2::labs(x = "Frequency", y = "Copy number mean") +
    cowplot::theme_cowplot()
  p
}
