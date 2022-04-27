#' Draw Oncoprint for fCNA Profile Visualization
#'
#' See [oncoprint](https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html)
#' and [heatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#size-of-the-heatmap)
#' to modify the result plot in depth.
#'
#' @param fCNA a [fCNA] object.
#' @param genes a list of genes. Must be subsets of `fCNA$data$gene_id`.
#' @param top_n top N genes to show.
#' @param top_n_by how to order the genes in `gene_summary` to get top N genes.
#' @param samples a vector to subset samples.
#' @param only_circular only show circular amplicon.
#' @param merge_circular if `TRUE`, merge 'possibly_circular' into 'circular' type.
#' @param show_column_names see `?ComplexHeatmap::oncoPrint`.
#' @param remove_empty_columns see `?ComplexHeatmap::oncoPrint`.
#' @param remove_empty_rows see `?ComplexHeatmap::oncoPrint`.
#' @param heatmap_legend_param see `?ComplexHeatmap::oncoPrint`.
#' @param col see `?ComplexHeatmap::oncoPrint`.
#' @param ... Other parameters passing to `ComplexHeatmap::oncoPrint`.
#'
#' @return a `Heatmap` object, you can go further modify it.
#' @export
#' @seealso [fCNA] for building object.
#'
#' @examples
#' \donttest{
#' if (require("ComplexHeatmap") && require("IDConverter")) {
#'   data("ascn")
#'   data <- ascn
#'
#'   # Create fake data
#'   set.seed(1234)
#'   data$sample <- sample(LETTERS[1:10], nrow(data), replace = TRUE)
#'   rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
#'
#'   gcap.plotProfile(rv, samples = c("B", "A", "C"))
#'
#'   rv$convertGeneID()
#'   ht <- gcap.plotProfile(rv,
#'     samples = c("B", "A", "C", "D", "F"),
#'     genes = unique(rv$data$gene_id)[1:10],
#'     top_annotation = ComplexHeatmap::HeatmapAnnotation(
#'       cbar = ComplexHeatmap::anno_oncoprint_barplot(),
#'       foo1 = c("g1", "g1", "g2", "g2", "g3"),
#'       bar1 = ComplexHeatmap::anno_points(1:5)
#'     )
#'   )
#'   ht
#'   ComplexHeatmap::draw(ht, merge_legends = TRUE)
#'
#'   #----------------
#'   # Plot KM curve
#'   #----------------
#'   surv_data <- data.frame(
#'     sample = rv$sample_summary$sample,
#'     time = 3000 * abs(rnorm(nrow(rv$sample_summary))),
#'     status = sample(c(0, 1), nrow(rv$sample_summary), replace = TRUE)
#'   )
#'   surv_data
#'   p <- gcap.plotKMcurve(rv, surv_data)
#'   p
#'
#'   p2 <- gcap.plotKMcurve(rv, surv_data, genes = "MYC")
#'   p2
#'
#'   # ---------------
#'   # Plot forest
#'   # ---------------
#'   p3 <- gcap.plotForest(rv, surv_data)
#'
#'   # Fake some data for gene analysis
#'   rv$data[amplicon_type == "circular"]$gene_id[1:10] <- "MYC"
#'   rv$data[amplicon_type == "circular"]$gene_id[11:20] <- "EGFR"
#'   rv$data[amplicon_type == "circular"]$gene_id[31:35] <- "GSX2"
#'   p4 <- gcap.plotForest(rv, surv_data,
#'     x = c("TP53", "MYC", "ABC", "EGFR", "GSX2"), x_is_gene = TRUE,
#'     ref_line = 1, xlim = c(0, 10)
#'   )
#'
#'   # Plot on gene clusters
#'   gcap.plotForest(rv, surv_data,
#'     x = data.frame(
#'       cluster = paste0("cluster", c(1, 3, 2, 1, 3)),
#'       gene_id = c("TP53", "MYC", "ABC", "EGFR", "GSX2")
#'     ), x_is_gene = TRUE,
#'     ref_line = 1, xlim = c(0, 10), optimize_model = TRUE
#'   ) -> zz
#'   # zz$plot
#'
#'   # -----------------
#'   # Plot distribution
#'   # -----------------
#'   gcap.plotDistribution(rv)
#'   p5 <- gcap.plotDistribution(rv,
#'     x = c("MYC", "EGFR", "CDK4", "AKT3"),
#'     width = 0.5, x_size = 5, fill = FALSE
#'   )
#'   p5
#'   rv$sample_summary[, ploidy_class := ifelse(ploidy > 2.5, "2+", "2")]
#'   p6 <- gcap.plotDistribution(rv, x = "ploidy_class", by = "sample")
#'   p6
#' }
#' }
#' @testexamples
#' expect_is(ht, "Heatmap")
#' expect_is(p, "ggsurvplot")
#' if (!is.null(p2)) {
#'   expect_is(p2, "ggsurvplot")
#' }
#' expect_is(p3, "list")
#' expect_is(p4, "list")
#' expect_is(p5, "ggplot")
#' expect_is(p6, "ggplot")
#' expect_is(zz, "list")
gcap.plotProfile <- function(fCNA,
                             genes = NULL,
                             samples = NULL,
                             top_n = NULL,
                             top_n_by = c("circular", "Total", "noncircular"),
                             only_circular = FALSE,
                             merge_circular = FALSE,
                             show_column_names = TRUE,
                             remove_empty_columns = FALSE,
                             remove_empty_rows = TRUE,
                             heatmap_legend_param = list(title = "fCNA"),
                             col = c(
                               "noncircular" = "#0066CC",
                               "circular" = "#CC0033",
                               "possibly_circular" = "#FFCCCC"
                             ),
                             ...) {
  stopifnot(inherits(fCNA, "fCNA"))
  .check_install("ComplexHeatmap")
  if (nrow(fCNA$data) == 0) {
    warning("No data to plot", immediate. = TRUE)
    return(NULL)
  }

  data <- data.table::copy(fCNA$data)[!is.na(gene_id)]
  if (only_circular) {
    data[, amplicon_type := data.table::fcase(
      amplicon_type %in% c("circular", "possibly_circular"), amplicon_type,
      default = NA
    )]
  }
  if (merge_circular) {
    data[, amplicon_type := data.table::fcase(
      amplicon_type %in% c("circular", "possibly_circular"), "circular",
      amplicon_type %in% "noncircular", "noncircular",
      default = NA
    )]
    data <- data[!is.na(amplicon_type)]
    data[, amplicon_type := factor(amplicon_type, c("noncircular", "circular"))]
  }
  data <- data.table::dcast(
    data,
    gene_id ~ sample,
    value.var = "amplicon_type", fill = NA,
    drop = FALSE, fun.aggregate = function(x) x[1]
  )

  if (!is.null(top_n)) {
    gene_info <- fCNA$gene_summary[!is.na(gene_id)]
    orders <- do.call("order", args = c(lapply(top_n_by, function(x) gene_info[[x]]), decreasing = TRUE))
    genes <- gene_info[orders]$gene_id[seq_len(top_n)]
  }

  data <- data.frame(data[, -1], row.names = data[[1]])
  if (!is.null(genes)) {
    data <- data[genes, ]
    if (nrow(data) == 0) {
      warning("No gene left to plot", immediate. = TRUE)
      return(NULL)
    }
  }

  if (!is.null(samples)) {
    data <- data[, samples, drop = FALSE]
  }

  alter_fun <- list(
    background = ComplexHeatmap::alter_graphic("rect", fill = "#CCCCCC"),
    circular = ComplexHeatmap::alter_graphic("rect", fill = col["circular"]),
    noncircular = ComplexHeatmap::alter_graphic("rect", fill = col["noncircular"]),
    possibly_circular = ComplexHeatmap::alter_graphic("rect", fill = col["possibly_circular"])
  )

  ht <- ComplexHeatmap::oncoPrint(data,
    alter_fun = alter_fun, col = col,
    heatmap_legend_param = heatmap_legend_param,
    remove_empty_columns = remove_empty_columns,
    remove_empty_rows = remove_empty_rows,
    show_column_names = show_column_names,
    ...
  )
  return(ht)
}
