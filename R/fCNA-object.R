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
  lock_objects = FALSE,
  lock_class = TRUE,
  public = list(
    #' @field data a `data.table` storing fCNA list, which typically contains following columns:
    #' - `sample` sample or case ID.
    #' - `band` chromosome cytoband.
    #' - `gene_id` gene ID, typically Ensembl ID. You can convert the ID with R package `IDConverter`.
    #' - `total_cn` total copy number value.
    #' - `minor_cn` copy number value for minor allele.
    #' - `prob` the probability the gene located in circular DNA.
    #' - `gene_class` gene level amplicon classification.
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
    data = NULL,
    sample_summary = NULL,
    #' @description Create a `fCNA` object.
    #' Typically, you can obtain this object from [gcap.workflow()] or [gcap.ASCNworkflow()].
    #' @param fcna a `data.frame` storing focal copy number amplicon list.
    #' @param pdata a `data.frame` storing phenotype or sample-level related data. (Optional)
    #' @param min_prob the minimal aggregated (in cytoband level) probability to determine a circular amplicon.
    #' @param only_oncogenes only_oncogenes if `TRUE`, only known oncogenes are kept for circular prediction.
    #' @param genome_build genome version
    initialize = function(fcna, pdata = fcna[, "sample", drop = FALSE],
                          min_prob = 0.6,
                          only_oncogenes = FALSE,
                          genome_build = c("hg38", "hg19", "mm10")) {
      stopifnot(
        is.data.frame(fcna), is.data.frame(pdata),
        all(c(
          "sample", "band", "gene_id", "total_cn",
          "ploidy", "prob", "gene_class"
        ) %in% colnames(fcna)),
        all(c("sample") %in% colnames(pdata)),
        !is.null(min_prob) && min_prob >= 0.5 && min_prob < 1
      )
      genome_build = match.arg(genome_build)

      private$prob <- min_prob
      amplvls <- c("noncircular", "circular")
      fcna <- fcna[fcna$gene_class %in% amplvls, ]
      fcna$gene_class <- factor(fcna$gene_class, levels = amplvls)
      self$data <- data.table::as.data.table(fcna)
      message("summarizing sample...")
      ss <- self$getSampleSummary(only_oncogenes = only_oncogenes,
                                  genome_build = genome_build)
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
      fCNA$new(fcna, pdata, min_prob = private$prob)
    },
    #' @description Get sample summary of fCNA
    #' @param only_oncogenes only_oncogenes if `TRUE`, only known oncogenes are kept for circular prediction.
    #' @param genome_build genome version.
    #' @return a `data.table`
    getSampleSummary = function(only_oncogenes = FALSE, genome_build = c("hg38", "hg19", "mm10")) {
      genome_build = match.arg(genome_build)
      message("  classifying samples with min_prob=", private$prob)
      if (only_oncogenes) message("  classifying 'circular' only based on oncogenes")
      self$data[
        , summarize_sample(.SD,
          min_prob = private$prob,
          only_oncogenes = only_oncogenes,
          genome_build = genome_build
        ),
        by = .(sample)
      ]
    },
    #' @description Get gene level summary of fCNA type
    #' @param return_mat if `TRUE`, return a cytoband by sample matrix instead of a summary.
    #' @return a `data.table` or a `matrix`.
    getGeneSummary = function(return_mat = FALSE) {
      data <- copy(self$data)
      allsamps <- self$sample_summary$sample
      getGeneSummary(data, allsamps, return_mat)
    },

    #' @description Get cytoband level summary of fCNA type
    #' @param unique if `TRUE`, count sample frequency instead of gene frequency.
    #' @param return_mat if `TRUE`, return a cytoband by sample matrix instead of a summary.
    #' @return a `data.table`
    getCytobandSummary = function(unique = FALSE,
                                  return_mat = FALSE) {
      data <- copy(self$data)
      allsamps <- self$sample_summary$sample
      getCytobandSummary(data, allsamps, unique, return_mat)
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
                             genome_build = c("hg38", "hg19", "mm10")) {
      if (!require("IDConverter")) {
        message("package 'IDConverter' is required to convert IDs")
        return(NULL)
      }
      opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
      options(IDConverter.datapath = opts)

      type <- match.arg(type)
      genome_build <- match.arg(genome_build)
      message("converting gene IDs, will update 'data' fields")
      self$data$gene_id <- IDConverter::convert_hm_genes(
        self$data$gene_id,
        type = type,
        genome_build = genome_build
      )
    },
    #' @description print the fCNA object
    #' @param ... unused.
    print = function(...) {
      ss <- self$sample_summary
      cat("======================\nA <")
      cat(cli::col_br_cyan("fCNA"))
      cat("> object\n")
      cat(sprintf("%8s: %s\n", "record", cli::col_green(nrow(self$data))))
      cat(sprintf("%8s: %s\n", "case", cli::col_cyan(nrow(ss))))
      cat("     |__ ",
        "(", cli::col_cyan(sum(ss$class == "noncircular", na.rm = TRUE)), ") ",
        cli::col_green(nrow(self$data[gene_class %in% "noncircular"])),
        " noncircular\n",
        sep = ""
      )
      cat("     |__ ",
        "(", cli::col_cyan(sum(ss$class == "circular", na.rm = TRUE)), ") ",
        cli::col_green(nrow(self$data[gene_class %in% "circular"])),
        " circular\n",
        sep = ""
      )
      cat("======================\n")
    }
  ),
  private = list(
    min_probs = NULL
  ),
  active = list(
    #' @field min_prob check `$new()` method for details. If you updated this value,
    #' a function will be called to update the sample summary.
    min_prob = function(x) {
      if (missing(x)) {
        private$prob
      } else {
        message("A new 'min_prob' is specified, try updating sample summary")
        private$prob <- x
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

summarize_sample <- function(data, min_prob,
                             only_oncogenes = FALSE,
                             genome_build = "hg38") {
  if (nrow(data) == 0) {
    return(data.frame(class = "nofocal"))
  }

  if (only_oncogenes) {
    if (startsWith(genome_build, "mm")) {
      ids <- readRDS(
        system.file(
          "extdata", "oncogenes_mouse.rds",
          package = "gcap", mustWork = TRUE
        )
      )
    } else {
      ids <- na.omit(unique(oncogenes$gene_id))
    }

    prob_possibly <- data[data$gene_id %in% ids & data$gene_class %in% "circular" & !is.na(data$prob)]
  } else {
    prob_possibly <- data[data$gene_class %in% "circular" & !is.na(data$prob)]
  }

  if (nrow(prob_possibly) == 0) {
    return(data.frame(class = "noncircular"))
  }

  # Use cytobands instead of genes to calculate
  prob_possibly <- prob_possibly[
    , .(prob = max(prob, na.rm = TRUE)),
    by = .(band)
  ]$prob

  if (max(prob_possibly, na.rm = TRUE) >= min_prob) {
    flag_ec <- TRUE
  } else {
    flag_ec <- calc_prob(prob_possibly, 2) >= min_prob
  }
  # flag_ec <- calc_prob(prob_possibly, 1) >= min_prob

  class <- if (flag_ec) {
    "circular"
  } else {
    "noncircular"
  }

  data.frame(class = class)
}


getGeneSummary <- function(data, allsamps, return_mat) {
  if (nrow(data) > 0) {
    rv <- data.table::dcast(data[!is.na(gene_id), .N, by = .(gene_id, gene_class)],
      gene_id ~ gene_class,
      value.var = "N", fill = 0L,
      drop = FALSE
    )
    rv$Total <- rowSums(rv[, -1])
    rv <- rv[order(rv$circular, rv$Total, rv$gene_id, decreasing = TRUE), ]
  } else {
    rv <- data.table::data.table()
  }
  if (!return_mat) {
    return(rv)
  }
  if (nrow(rv) == 0) {
    message("No data")
    return(NULL)
  }

  data$sample <- factor(data$sample, allsamps)
  data <- data.table::dcast(
    data[!is.na(gene_id)],
    gene_id ~ sample,
    value.var = "gene_class", fill = NA,
    drop = FALSE, fun.aggregate = function(x) x[1]
  )
  data <- data.frame(data[, -1], row.names = data[[1]])
  # Reorder data
  data <- data[rv$gene_id, , drop = FALSE]
  return(data)
}

getCytobandSummary <- function(data, allsamps, unique, return_mat) {
  if (nrow(data) > 0) {
    if (unique) {
      data2 <- data[!is.na(band), .(N = length(unique(sample))), by = .(band, gene_class)]
    } else {
      data2 <- data[!is.na(band), .N, by = .(band, gene_class)]
    }

    rv <- data.table::dcast(data2,
      band ~ gene_class,
      value.var = "N", fill = 0L,
      drop = FALSE
    )
    rv$Total <- rowSums(rv[, -1])
    rv <- rv[order(rv$circular, rv$Total, rv$band, decreasing = TRUE), ]
  } else {
    rv <- data.table::data.table()
  }
  if (!return_mat) {
    return(rv)
  }
  if (nrow(rv) == 0) {
    message("No data")
    return(NULL)
  }

  data$sample <- factor(data$sample, allsamps)
  data <- data.table::dcast(
    data[!is.na(band)],
    band ~ sample,
    value.var = "gene_class", fill = NA,
    drop = FALSE, fun.aggregate = function(x) {
      if ("circular" %in% x) {
        return("circular")
      }
      if ("noncircular" %in% x) {
        return("noncircular")
      }
      return(NA)
    }
  )
  data <- data.frame(data[, -1], row.names = data[[1]])
  # Reorder data
  data <- data[rv$band, , drop = FALSE]
  return(data)
}
