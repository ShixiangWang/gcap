#' Draw Gene Amplicon Distribution By Gene or User Specified Class
#'
#' See [gcap.plotProfile()] for examples.
#'
#' @inheritParams gcap.plotProfile
#' @param x a list of gene or a column name in `fCNA$sample_summary`.
#' @param x_size font size of x axis text.
#' @param by "gene" (default) or "sample" mode. When "sample" is set, it uses the column specified
#' by `x` to grouping samples.
#' @param fill if `TRUE`, show the percentage instead of count.
#' @param palette color palette.
#' @param ... other parameters passing to `ggplot2::geom_bar`.
#'
#' @return a ggplot object.
#' @export
#' @seealso [gcap.plotProfile] for plot landscape of fCNA, [fCNA] for building object.
gcap.plotDistribution <- function(fCNA,
                                  x = NULL,
                                  x_size = 8,
                                  by = c("gene", "sample"),
                                  merge_circular = TRUE,
                                  fill = TRUE,
                                  palette = c("#CCCCCC", "#0066CC", "#FFCCCC", "#CC0033"),
                                  ...) {
  stopifnot(inherits(fCNA, "fCNA"))
  .check_install("ggplot2")

  by <- match.arg(by)
  nsamples <- nrow(fCNA$sample_summary)

  if (by == "gene") {
    data <- fCNA$data[, c("sample", "amplicon_type", "gene_id")]
    colnames(data)[2] <- "class"
    if (!is.null(x)) {
      data <- data[gene_id %in% x]
      data$gene_id <- factor(data$gene_id, x)
    }
  } else {
    if (length(x) != 1 && !is.character(x)) {
      warning("When plot by 'sample', x should be a column name of `fCNA$sample_summary`", immediate. = TRUE)
      return(NULL)
    }
    data <- fCNA$sample_summary[, c("sample", "class", x), with = FALSE]
  }
  colnames(data)[3] <- "by"

  if (merge_circular) {
    data[, class := data.table::fcase(
      class %in% c("circular", "possibly_circular"), "circular",
      class == "noncircular", "noncircular",
      default = "nofocal"
    )]
    class_lvls <- c("nofocal", "noncircular", "circular")
  } else {
    class_lvls <- c("nofocal", "noncircular", "possibly_circular", "circular")
  }
  data[, class := factor(class, class_lvls)]

  dt_n <- data[, .N, by = list(by, class)]

  if (by == "gene") {
    dt_nofocal <- dt_n[, list(N = sum(N)), by = list(by)]
    dt_nofocal$class <- "nofocal"
    dt_nofocal$N <- nsamples - dt_nofocal$N
    if (any(dt_nofocal$N < 0)) {
      message("A negative value is obtained for gene nonfocal counts, if you are not runing an example, please check your data or options")
      dt_nofocal$N <- ifelse(dt_nofocal$N < 0, 0L, dt_nofocal$N)
    }
    dt_nofocal <- dt_nofocal[, c("by", "class", "N")]
    dt_n <- rbind(dt_n, dt_nofocal)
  }

  if (fill) {
    dt_n <- dt_n[, list(class, N = N / sum(N)), by = list(by)]
  }


  if (length(class_lvls) == 3 && identical(palette, c("#CCCCCC", "#0066CC", "#FFCCCC", "#CC0033"))) {
    palette <- palette[-3]
  }

  p <- ggplot2::ggplot(dt_n, ggplot2::aes(by, N, fill = class)) +
    ggplot2::geom_bar(stat = "identity", ...) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion()) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::scale_x_discrete(position = "top") + # guide = ggplot2::guide_axis(n.dodge = 2)
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        size = x_size,
        vjust = 0,
        hjust = 0,
        color = "black",
        face = "bold"
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL)

  p
}
