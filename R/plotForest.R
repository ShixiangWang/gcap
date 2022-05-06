#' Draw Forest Plot for fCNA+/- Comparison by Sample or Gene
#'
#' See [gcap.plotProfile()] for examples.
#'
#' @inheritParams gcap.plotKMcurve
#' @param response_data a `data.frame` or specified column names in `fCNA$sample_summary` prepared
#' regression modeling, default is CoxPH model, you can use same data as [gcap.plotKMcurve].
#' @param f a length-1 string specifying modeling function or family of [glm()], default is 'coxph'. Other options are members of GLM family, see [stats::family()].
#' 'binomial' is logistic, and 'gaussian' is linear.
#' @param x focal variables (terms), the column names of `fCNA$sample_summary` or a list of genes (vector).
#' @param covars covariables in `fCNA$sample_summary`.
#' @param y predicted variables or expression in string format,
#' e.g., `"sex"` or `"factor(sex)"` is acceptable.
#' @param x_is_gene if `TRUE`, treat `x` as a list of genes.
#' @param gene_focus focal amplication type you focus on.
#' Can be 'fCNA' or 'circular' and 'vs'. If 'fCNA' selected,
#' noncircular and circular genes are included to classify samples.
#' 'vs' for comparing circular and noncircular. (at least 2 samples in each group)
#' @param optimize_model if `TRUE`, optimize the variable list with `stepAIC` in
#' 'both' direction.
#' @param parallel if `TRUE`, use N-1 cores to run the task.
#' @param exp logical, indicating whether or not to exponentiate the the coefficients.
#' @param ref_line reference line, default is 1 for HR.
#' @param xlim limits of x axis.
#' @param ... other plot options passing to `?forestploter::forest()`.
#' Also check <https://github.com/adayim/forestploter> to see more complex adjustment of the result plot.
#'
#' @return a list containing model and forest plot.
#' @export
gcap.plotForest <- function(fCNA,
                            response_data,
                            f = c(
                              "coxph", "binomial", "gaussian", "Gamma", "inverse.gaussian", "poisson",
                              "quasi", "quasibinomial", "quasipoisson"
                            ),
                            x = "class",
                            covars = NULL,
                            y = if (is.data.frame(response_data)) {
                              colnames(response_data)[-1]
                            } else {
                              response_data
                            },
                            merge_circular = TRUE,
                            x_is_gene = FALSE,
                            gene_focus = c("fCNA", "circular", "vs"),
                            optimize_model = FALSE,
                            ending_time = NULL,
                            parallel = FALSE,
                            exp = NULL,
                            ref_line = NULL,
                            xlim = NULL,
                            ...) {
  stopifnot(inherits(fCNA, "fCNA"))
  .check_install("regport")
  force(y)
  f <- f[1]
  gene_focus <- match.arg(gene_focus)

  if (is.character(response_data)) {
    response_data <- fCNA$sample_summary[, c("sample", response_data), with = FALSE]
  }

  if (!x_is_gene) {
    data <- fCNA$sample_summary[, c("sample", c(x, covars)), with = FALSE]
    if (merge_circular & "class" %in% c(x, covars)) {
      data[, class := data.table::fcase(
        class %in% c("circular", "possibly_circular"), "circular",
        class == "noncircular", "noncircular",
        default = "nofocal"
      )]
      data[, class := factor(class, c("nofocal", "noncircular", "circular"))]
    } else {
      if ("class" %in% c(x, covars)) {
        data[, class := factor(class, c("nofocal", "noncircular", "possibly_circular", "circular"))]
      }
    }
  } else {
    # Extract class based on gene
    all_samples <- fCNA$sample_summary$sample
    # AMP samples
    if (gene_focus == "fCNA") {
      types <- c("noncircular", "possibly_circular", "circular")
    } else if (merge_circular) {
      types <- c("possibly_circular", "circular")
    } else {
      types <- "circular"
    }
    labels <- c("-", "+")

    data <- fCNA$data
    # If input a gene cluster instead of a gene list
    if (is.data.frame(x)) {
      message("found input a data.frame when x is gene, treat as gene clusters")
      colnames(x)[1:2] <- c("cluster", "gene_id")
      cMap <- x$cluster
      names(cMap) <- x$gene_id
      data$gene_id <- cMap[data$gene_id]
      data$gene_id <- ifelse(is.na(data$gene_id), ".others", data$gene_id)
      x <- unique(x$cluster)
    }

    dt_list <- list()
    xc <- x
    for (gene in xc) {
      message("Classifying based on gene/cluster ", gene)
      amp_samples <- unique(data[gene_id %in% gene & amplicon_type %in% types]$sample)
      if (gene_focus == "vs") {
        all_samples <- unique(
          data[
            gene_id %in% gene & amplicon_type %in% c("noncircular", types)
          ]$sample
        )
      }
      dt <- create_label_dt(amp_samples, all_samples, labels)
      if (is.null(dt) || (gene_focus == "vs" && min(table(dt[[2]])) < 2) || min(table(dt[[2]])) < 1) {
        message("No proper data for gene ", gene, " skipping.")
        x <- setdiff(x, gene)
        next()
      }
      colnames(dt)[2] <- gene
      dt_list[[gene]] <- dt
    }

    if (length(dt_list) == 0) {
      message("no data left to analyze, for type 'vs', at least 2 for each group, other at least 1")
      return(invisible(NULL))
    } else if (length(x) == 1) {
      data <- dt_list[[1]]
    } else {
      data <- mergeDTs(dt_list, by = "sample")
    }
  }

  response_data <- data.table::as.data.table(response_data)
  if (!is.null(ending_time)) {
    # Assume the 2nd is time and 3rd is status
    response_data[[2]] <- ifelse(response_data[[2]] >= ending_time, ending_time, response_data[[2]])
    response_data[[3]] <- ifelse(response_data[[2]] >= ending_time, 0, response_data[[3]])
  }

  data <- merge(data, response_data,
    by = "sample",
    all.x = TRUE
  )

  # 这里聚合的数据可以做任意合适的分析

  if (optimize_model) {
    # 按照影响的样本大小排序，然后使用forward
    message("try constructing an optimized model")
    x2 <- x[order(sapply(data[, x, with = FALSE], function(x) sum(table(x)[-1], na.rm = TRUE)), decreasing = TRUE)] # 有2级以上的因子的数目
    optmodel <- tryCatch(
      {
        mo <- regport::REGModel$new(data, recipe = list(
          x = x, y = y
        ))
        stats::step(mo$model, steps = 1e4)
      },
      error = function(e) {
        message("Failed in building the optimized model, error:")
        message(e$mess)
        NULL
      }
    )
  } else {
    optmodel <- NULL
  }

  ml <- regport::REGModelList$new(data, y = y, x = x, covars = covars)
  ml$build(f, exp = exp, parallel = parallel)
  ml$print()

  p <- tryCatch(
    ml$plot_forest(
      ref_line = ref_line,
      xlim = xlim,
      ...
    ),
    error = function(e) {
      message("\nPlotting failed. Check error message and model result for details")
      message(e$message)
      NULL
    }
  )
  invisible(list(data = data, plot = p, model = ml, optmodel = optmodel))
}

create_label_dt <- function(amp_samples, all_samples, labels) {
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
  data
}
