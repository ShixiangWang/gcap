#' R6 class representing focal copy number amplification list predicted from a cohort
#'
#' @description
#' Contains fields storing data and methods to get, process and visualize
#' fCNA information. Examples please see [gcap.ASCNworkflow()].
#'
#' @export
fCNA <- R6::R6Class(
  "fCNA",
  inherit = NULL,
  cloneable = FALSE,
  lock_objects = TRUE,
  lock_class = TRUE,
  public = list(
    #' @field data a `data.table` storing fCNA list, which typically contains following columns:
    #' - `sample` sample or case ID.
    #' - `band` chromosome cytoband.
    #' - `gene_id` gene ID, typically Ensembl ID. You can convert the ID with R package `IDConverter`.
    #' - `total_cn` total copy number value.
    #' - `minor_cn` copy number value for minor allele.
    #' - `background_cn` background copy number value (for circular amplicon), calculated as `pmax(blood_cn_top5 + tightness * blood_cn_top5_sd, 2) * ploidy / 2`. (**Optional**)
    #' At default, `blood_cn_top5(_sd)` is obtained from SNP array data of ~2000 TCGA diploidy blood samples,
    #' they are mean and sd of top 5% copy number values for each gene (in ~2000 samples);
    #' `ploidy` is tumor ploidy.
    #' - `prob` the probability the gene located in circular DNA.
    #' - `amplicon_type` the type of DNA amplicon.
    #' @field sample_summary a `data.table` storing sample summary data, which typically contains
    #' at least the following columns:
    #' - `sample` sample or case ID. **Should only include cases have been called with GCAP workflow,
    #' otherwise the extra cases would be automatically classified as 'nofocal' (i.e. `NA` in `sample_summary` field) class**.
    #' - `purity`, `ploidy` for tumor purity or ploidy.
    #' - `AScore` aneuploidy score.
    #' - `pLOH` genome percentage harboring LOH events.
    #' - `CN1 ... CN19` activity of copy number signatures.
    #' - **`class`** the sample class based on amplicon type.
    #' - `ec_genes` number of genes predicted as located on circular DNA.
    #' - `ec_possibly_genes` same with `ec_genes` but with less confidence.
    #' - `ec_cytobands` number of cytobands predicted as located on circular DNA.
    #' (the regions of `ec_possibly_genes` are not included in computation)
    #' @field gene_summary a `data.table` storing gene summary data.
    #' @field cytoband_summary a `data.table` storing cytoband summary data.
    data = NULL,
    sample_summary = NULL,
    gene_summary = NULL,
    cytoband_summary = NULL,
    #' @description Create a `fCNA` object.
    #' Typically, you can obtain this object from [gcap.workflow()] or [gcap.ASCNworkflow()].
    #' @param fcna a `data.frame` storing focal copy number amplicon list.
    #' @param pdata a `data.frame` storing phenotype or sample-level related data. (Optional)
    #' @param min_n a minimal cytoband number (default is `1`) to determine
    #' sample class. e.g., sample with at least 1 cytoband harboring circular
    #' genes would be labelled as "circular". **This only affect `sample_summary` field**.
    initialize = function(fcna, pdata = fcna[, "sample", drop = FALSE], min_n = 1L) {
      stopifnot(
        is.data.frame(fcna), is.data.frame(pdata),
        all(c(
          "sample", "band", "gene_id", "total_cn",
          "prob", "amplicon_type"
        ) %in% colnames(fcna)),
        all(c("sample") %in% colnames(pdata)),
        !is.null(min_n) && min_n >= 1
      )

      private$n <- min_n

      message("set fixed levels for 'amplicon_type'")
      amplvls <- c("noncircular", "possibly_circular", "circular")
      fcna <- fcna[fcna$amplicon_type %in% amplvls, ]
      fcna$amplicon_type <- factor(fcna$amplicon_type, levels = amplvls)
      self$data <- data.table::as.data.table(fcna)
      message("summarizing sample...")
      ss <- self$getSampleSummary()
      pdata <- data.table::as.data.table(pdata)
      overlap_cols <- setdiff(intersect(colnames(pdata), colnames(ss)), "sample")
      if (length(overlap_cols) > 0) {
        message("  overlapping columns detected, remove them from input 'pdata'")
        pdata <- pdata[, !colnames(pdata) %in% overlap_cols, with = FALSE]
      }
      self$sample_summary <- merge(pdata, ss,
        by = "sample", all = TRUE
      )
      # Set default label 'nofocal' to samples with NA
      self$sample_summary[, `:=`(class = data.table::fifelse(is.na(class), "nofocal", class))]

      message("summarizing gene and cytoband...")
      self$gene_summary <- self$getGeneSummary()
      self$cytoband_summary <- self$getCytobandSummary()
      message("done")
    },
    #' @description Return a subset `fCNA` object
    #' @param ... subset expressions on `fCNA$data` or `fCNA$sample_summary`.
    #' @param on if it is "data", subset operations are on data field of `fCNA` object,
    #' same for "sample_summary".
    #' @return a `fCNA`
    subset = function(..., on = c("data", "sample_summary")) {
      on <- match.arg(on)
      data <- if (on == "data") self$data else self$sample_summary
      subs <- unique(subset(data, ...)$sample)
      fcna <- if (on == "data") subset(data, ...) else self$data[sample %in% subs]
      pdata <- self$sample_summary[sample %in% subs]
      fCNA$new(fcna, pdata, min_n = private$n)
    },
    #' @description Get sample summary of fCNA
    #' @return a `data.table`
    getSampleSummary = function() {
      stopifnot(!is.null(private$n) && private$n >= 1L)
      message("  classifying samples with min_n=", private$n)
      self$data[
        , summarize_sample(.SD,
          min_n = private$n
        ),
        by = .(sample)
      ]
    },
    #' @description Get gene level summary of fCNA type
    #' @return a `data.table`
    getGeneSummary = function() {
      if (nrow(self$data) > 0) {
        rv <- data.table::dcast(self$data[, .N, by = .(gene_id, amplicon_type)],
          gene_id ~ amplicon_type,
          value.var = "N", fill = 0L,
          drop = FALSE
        )
        rv$Total <- rowSums(rv[, -1])
      } else {
        rv <- data.table::data.table()
      }
      rv
    },
    #' @description Get cytoband level summary of fCNA type
    #' @return a `data.table`
    getCytobandSummary = function() {
      if (nrow(self$data) > 0) {
        rv <- data.table::dcast(self$data[, .N, by = .(band, amplicon_type)],
          band ~ amplicon_type,
          value.var = "N", fill = 0L,
          drop = FALSE
        )
        rv$Total <- rowSums(rv[, -1])
      } else {
        rv <- data.table::data.table()
      }
      rv
    },
    #' @description Save the key data to local files
    #' @param dirpath directory path storing output files.
    #' @param fileprefix file prefix. Two result files shall be generated.
    saveToFiles = function(dirpath, fileprefix = "fCNA") {
      fl_records <- file.path(dirpath, paste0(fileprefix, "_fCNA_records.csv"))
      fl_sample <- file.path(dirpath, paste0(fileprefix, "_sample_info.csv"))

      data.table::fwrite(self$data, file = fl_records)
      data.table::fwrite(self$sample_summary, file = fl_sample)
    },
    #' @description Convert Gene IDs between Ensembl and Hugo Symbol System
    #' @param type type of input IDs, could be 'ensembl' or 'symbol'.
    #' @param genome_build reference genome build.
    convertGeneID = function(type = c("ensembl", "symbol"),
                             genome_build = c("hg38", "hg19", "mm10", "mm9")) {
      if (!require("IDConverter")) {
        message("package 'IDConverter' is required to convert IDs")
        return(NULL)
      }
      opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
      options(IDConverter.datapath = opts)

      type <- match.arg(type)
      genome_build <- match.arg(genome_build)
      message("converting gene IDs, will update 'data' and 'gene_summary' fields")
      self$data$gene_id <- IDConverter::convert_hm_genes(
        self$data$gene_id,
        type = type,
        genome_build = genome_build
      )
      self$gene_summary$gene_id <- IDConverter::convert_hm_genes(
        self$gene_summary$gene_id,
        type = type,
        genome_build = genome_build
      )
    },
    #' @description print the fCNA object
    #' @param ... unused.
    print = function(...) {
      tbl <- table(self$data$amplicon_type)
      ss <- self$sample_summary
      cat("======================\nA <")
      cat(cli::col_br_cyan("fCNA"))
      cat("> object\n")
      cat(sprintf("%8s: %s\n", "case", cli::col_cyan(nrow(ss))))
      cat(sprintf("%8s: %s\n", "gene", cli::col_green(nrow(self$gene_summary))))
      cat(sprintf("%8s: %s\n", "fCNA", cli::col_green(nrow(self$data))))
      cat("     |__ ", cli::col_green(tbl[1]),
        " (", cli::col_cyan(sum(ss$class == "noncircular", na.rm = TRUE)), ") noncircular\n",
        sep = ""
      )
      cat("     |__ ", cli::col_green(tbl[2]),
        " (", cli::col_cyan(sum(ss$class == "possibly_circular", na.rm = TRUE)), ") possibly_circular\n",
        sep = ""
      )
      cat("     |__ ", cli::col_green(tbl[3]),
        " (", cli::col_cyan(sum(ss$class == "circular", na.rm = TRUE)), ") circular\n",
        sep = ""
      )
      cat("======================\n")
    }
  ),
  private = list(
    n = NULL
  ),
  active = list(
    #' @field min_n check `$new()` method for details. If you updated this value,
    #' a function will be called to update the sample summary.
    min_n = function(x) {
      if (missing(x)) {
        private$n
      } else {
        message("A new 'min_n' is specified, try updating sample summary")
        private$n <- x
        # Update sample summary
        ss <- self$getSampleSummary()
        pdata <- self$sample_summary
        overlap_cols <- setdiff(intersect(colnames(pdata), colnames(ss)), "sample")
        if (length(overlap_cols) > 0) {
          pdata <- pdata[, !colnames(pdata) %in% overlap_cols, with = FALSE]
        }
        self$sample_summary <- merge(pdata, ss, by = "sample", all = TRUE)
        message("done")
      }
    }
  )
)

# summarize sample --------------------------------------------------------

summarize_sample <- function(data, min_n) {
  ec_genes <- sum(data$amplicon_type == "circular", na.rm = TRUE)
  ec_cytobands_detail <- na.omit(unique(data$band[data$amplicon_type == "circular"]))
  ec_cytobands <- length(ec_cytobands_detail)
  ec_cytobands_detail <- paste(sort(ec_cytobands_detail), collapse = ",")
  ec_possibly_genes <- sum(data$amplicon_type == "possibly_circular", na.rm = TRUE)
  prob_possibly <- data[data$amplicon_type %in% c("possibly_circular", "circular")]

  if (nrow(prob_possibly) > 0) {
    # Collapse probs by cytobands
    prob_possibly <- prob_possibly[
      , .(prob = max(prob, na.rm = TRUE)),
      by = .(band)
    ]$prob
  } else {
    prob_possibly <- NULL
  }

  # flags have priority
  # flag_ec <- ec_genes >= min_n
  flag_ec <- ec_cytobands >= min_n
  flag_ec_possibly <- if (!flag_ec && length(prob_possibly) > 0) {
    calc_prob(prob_possibly, min_n) > 0.5
  } else {
    FALSE
  }

  class <- if (flag_ec) {
    "circular"
  } else if (flag_ec_possibly) {
    "possibly_circular"
  } else {
    # Use cytobands instead of genes to count
    flag_amp <- length(na.omit(unique(
      data$band[data$amplicon_type %in% c("noncircular", "possibly_circular", "circular")]
    ))) >= 1
    if (flag_amp) {
      "noncircular"
    } else {
      "nofocal"
    }
  }

  data.frame(
    class = class,
    ec_genes = ec_genes,
    ec_possibly_genes = ec_possibly_genes,
    ec_cytobands = ec_cytobands,
    ec_cytobands_detail = ec_cytobands_detail
  )
}


# Archive -----------------
# Checker comes from https://github.com/r-lib/R6/issues/48#issuecomment-339872876
# #' R6 class representing a Asserter
# #'
# #' @description
# #' This class includes many publc functions for validating
# #' the modification of a R6 object.
# #' @export
# Asserter <- R6::R6Class(
#   "Asserter",
#   inherit = NULL,
#   cloneable = FALSE,
#   lock_objects = TRUE,
#   lock_class = TRUE,
#   public = list(
#     #' @description
#     #' Assert input is a `data.table`
#     #' @param x modified/updated R6 field/active-binding
#     #' @param varnm corresponding variable name of input `x`, e.g.,
#     #' for command `obj$year <- 10`, `10` means `x` and `year` means `varnm`.
#     assertDT = function(x, varnm = deparse(substitute(x))) {
#       `if`(!all(data.table::is.data.table(x)), stop(sprintf("'%s' must be data.table", varnm), call. = FALSE), invisible(x))
#     }
#   )
# )
