#' Quantify sample-level ecDNA measueres
#'
#' @param data a `data.table` containing result from [gcap.runPrediction].
#' The column storing prediction result must start with `pred`.
#' @param cutoff a cutoff for converting prob into `0/1` value.
#'
#' @return a `data.frame`
#' @export
#'
#' @examples
#' data("ec")
#' ec2 <- ec
#' ec2$pred <- gcap.runPrediction(ec)
#' @testexamples
#'
gcap.runScoring <- function(data, cutoff = 0.9,
                            measures = c("load", "burden", "agglomeration"),
                            genome_build = c("hg38", "hg19")) {
  stopifnot(is.data.frame(data))
  measures <- match.arg(measures, several.ok = TRUE)
  genome_build <- match.arg(genome_build)
  lg <- set_logger()

  lg$info("checking input data type")
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }

  lg$info("checking columns")
  scs <- colnames(data)[startsWith(colnames(data), "pred")]
  if (length(scs) < 1) {
    lg$fatal("no columns start with 'pred' found, please check your input")
    stop()
  }
  if (!"sample" %in% colnames(data)) {
    lg$info("No 'sample' column found, mutating one")
    data$sample <- "sample"
  }

  lg$info("converting prob")
  scs2 <- paste0(scs, "_binary")
  for (i in seq_along(scs)) {
    data[[scs2[i]]] <- ifelse(data[[scs[i]]] > cutoff, 1L, 0L)
  }

  # Calculate load
  dt_load <- data[,  lapply(.SD, sum), .SDcols = scs2, by = "sample"]
  colnames(dt_load) <- sub("binary", "load", colnames(dt_load))

  # Calculate burden (per Mb)
  scs_burden <- paste0(scs, "_burden")
  dt_burden <- data[, lapply(scs2, function(x) {
    print(x)
    calc_burden(.SD[, c(x, "gene_id"), with = FALSE])
    sum(rnorm(nrow(.SD)))
  }), by = "sample"]

  # Calculate agglomeration index
  scs_agglomeration <- paste0(scs, "_agglomeration")
}

calc_burden <- function(dt) {

}
