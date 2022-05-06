#' Draw K-M Curve for fCNA Survival Comparison by Sample or Gene
#'
#' See [gcap.plotProfile()] for examples.
#'
#' @inheritParams gcap.plotProfile
#' @param surv_data survival data, eithor a 3-column `data.frame` to store
#' sample, time and status, or a length-2 string to specify the colnames
#' representing time and status in `fCNA$sample_summary`.
#' - sample must be identical to sample ID in `fCNA`.
#' - time must be numeric.
#' - status must be 0 or 1.
#' @param gene_focus focal amplication type you focus on.
#' Can be 'fCNA' or 'circular'. If 'fCNA' selected,
#' noncircular and circular genes are included to classify samples.
#' @param palette plot color palette.
#' @param class_col column name in `sample_summary` field for classification.
#' If you set to other column (you want to run survival analysis with custom column),
#' parameters like `merge_circular`, `genes`, `gene_focus`
#' etc. will be omitted.
#' @param ending_time survival analysis ending time. If a numeric ending
#' is typed, all survival data longer than the ending time will be rewritten.
#' @param ... other parameters passing to `survminer::ggsurvplot`.
#'
#' @return a plot.
#' @export
#' @seealso [gcap.plotProfile] for plot landscape of fCNA, [fCNA] for building object.
gcap.plotKMcurve <- function(fCNA,
                             surv_data,
                             merge_circular = TRUE,
                             genes = NULL,
                             gene_focus = c("fCNA", "circular"),
                             palette = c("grey", "#0066CC", "#CC0033"),
                             class_col = "class",
                             ending_time = NULL,
                             ...) {
  stopifnot(inherits(fCNA, "fCNA"))
  .check_install("survminer")
  gene_focus <- match.arg(gene_focus)

  if (is.character(surv_data)) {
    surv_data <- fCNA$sample_summary[, c("sample", surv_data), with = FALSE]
  }
  colnames(surv_data)[2:3] <- c("time", "status")

  if (is.null(genes)) {
    data <- fCNA$sample_summary[, c("sample", class_col), with = FALSE]
    colnames(data)[2] = "class"
    if (merge_circular & class_col == "class") {
      data[, class := data.table::fcase(
        class %in% c("circular", "possibly_circular"), "circular",
        class == "noncircular", "noncircular",
        default = "nofocal"
      )]
      data[, class := factor(class, c("nofocal", "noncircular", "circular"))]
    } else if (class_col == "class") {
      data[, class := factor(class, c("nofocal", "noncircular", "possibly_circular", "circular"))]
    }
  } else {
    # Extract class based on gene
    all_samples <- fCNA$sample_summary$sample
    # AMP samples
    if (gene_focus == "fCNA") {
      types <- c("noncircular", "possibly_circular", "circular")
      labels <- c("fCNA-", "fCNA+")
    } else if (merge_circular) {
      types <- c("possibly_circular", "circular")
      labels <- c("circular-", "circular+")
    } else {
      types <- "circular"
      labels <- c("circular-", "circular+")
    }
    amp_samples <- unique(fCNA$data[gene_id %in% genes & amplicon_type %in% types]$sample)
    if (length(amp_samples) == 0L) {
      warning("No sample amplified based on input genes and options", immediate. = TRUE)
      return(NULL)
    }
    data <- data.table::data.table(
      sample = c(amp_samples, setdiff(all_samples, amp_samples)),
      class = c(
        rep(labels[2], length(amp_samples)),
        rep(labels[1], length(setdiff(all_samples, amp_samples)))
      )
    )
    data$class <- factor(data$class, levels = labels)
  }

  data <- merge(data, data.table::as.data.table(surv_data),
    by = "sample",
    all.x = TRUE
  )

  if (!is.null(ending_time)) {
    data[, status := ifelse(time >= ending_time, 0, status)]
    data[, time := ifelse(time >= ending_time, ending_time, time)]
  }

  fit <- survminer::surv_fit(survival::Surv(time, status) ~ class, data = data)
  print(survminer::surv_pvalue(fit = fit))

  cls_lvls <- sub("class=", "", names(fit$strata))
  if (identical(palette, c("grey", "#0066CC", "#CC0033"))) {
    if (length(cls_lvls) == 2) palette <- palette[2:3]
  }

  p <- survminer::ggsurvplot(fit,
    pval = TRUE, data = data,
    palette = palette,
    risk.table = TRUE,
    legend.labs = cls_lvls,
    ...
  )
  p
}
