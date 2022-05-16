# Reference: https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html#msigdb-analysis
# https://igordot.github.io/msigdbr/articles/msigdbr-intro.html
# extra visualization: https://github.com/noriakis/CBNplot
# plot GSEA: https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01],
#                                       examplePathways, exampleRanks)
# mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
#                          order(-NES), pathway]
# plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes,
#               gseaParam = 0.5)

#' Enrich analysis
#'
#' @param data a [fCNA] object, or a vector (list) of genes.
#' If the vector is named, the names must be genes and values must numeric.
#' @param target `NA` for general, and "circular" or "noncircular" for specific fCNA type.
#' @param class_by if a [fCNA] object provided, this specify a column name in `fCNA$sample_summary` to
#' classify genes into different sample groups for enrichment. NOTE: any group with <10 genes is ignored.
#' @param analysis_func analysis function to run.
#' @param gene_encode gene ID type from input.
#' @param species species, not to change.
#' @param category gene set main category name.
#' @param subcategory gene set sub-category name.
#' @param genome_build genome version.
#' @param ... other parameters passing to `clusterProfiler::enricher()`
#' or `fgsea::fgsea`.
#'
#' @return depends on `data` and `analysis_func`.
#' @export
gcap.enrich <- function(data,
                        target = c(NA, "circular", "noncircular"),
                        class_by = NULL,
                        analysis_func = c("enricher", "fgsea"),
                        gene_encode = c("ensembl", "symbol"),
                        species = "Homo sapiens", category = "H", subcategory = "",
                        genome_build = c("hg38", "hg19"),
                        ...) {
  .check_install("fgsea", bioc = TRUE)
  .check_install("clusterProfiler", bioc = TRUE)
  .check_install("msigdbr")
  message("Check `msigdbr::msigdbr_collections()` for category list")

  analysis_func <- match.arg(analysis_func)
  gene_encode <- match.arg(gene_encode)
  genome_build <- match.arg(genome_build)

  msigdbr_df <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory) %>%
    data.table::as.data.table()

  if (gene_encode == "ensembl") {
    cols <- c("gs_name", "ensembl_gene")
  } else {
    cols <- c("gs_name", "gene_symbol")
  }

  msigdbr_df <- msigdbr_df[, cols, with = FALSE]
  colnames(msigdbr_df) <- c("term", "gene")

  target <- match.arg(target)
  if (inherits(data, "fCNA")) {
    if (analysis_func == "gsea") stop("gsea is not supported for fCNA object, please input a ranked list")
    dts <- data$sample_summary
    data <- data.table::copy(data$data)
    data[, amplicon_type := set_default_factor(amplicon_type)]

    if (!is.na(target)) {
      data <- data[data$amplicon_type %in% target]
    }

    if (!is.null(class_by)) {
      dts <- dts[, c("sample", class_by), with = FALSE]
      colnames(dts)[2] <- "class"
      data <- merge(data, dts, by = "sample", all.x = TRUE)
      geneList <- lapply(split(data, f = data$class, drop = TRUE), function(dt) {
        na.omit(unique(dt$gene_id))
      })
      # Only keep gene list with >=10 genes
      geneList <- geneList[sapply(geneList, length) >= 10]
      if (length(geneList) < 1) {
        warning("No gene list left to run enrichment analysis", immediate. = TRUE)
        return(invisible(NULL))
      }
    } else {
      geneList <- lapply(split(data, f = data$amplicon_type, drop = TRUE), function(dt) {
        na.omit(unique(dt$gene_id))
      })
    }
    if (is.list(geneList) && length(geneList) == 1) geneList <- geneList[[1]]
  } else {
    # A list/vector of (ranked) gene list
    geneList <- data
  }

  if (is.list(geneList)) {
    geneList <- lapply(geneList, genList)
  } else {
    geneList <- genList(geneList)
  }

  if (analysis_func == "enricher") {
    if (is.list(geneList)) {
      clusterProfiler::compareCluster(
        geneCluster = lapply(geneList, names),
        fun = clusterProfiler::enricher,
        TERM2GENE = msigdbr_df,
        # universe = y,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        ...
      )
    } else {
      clusterProfiler::enricher(
        gene = names(geneList), TERM2GENE = msigdbr_df,
        # universe = y,
        pvalueCutoff = 1, qvalueCutoff = 1, ...
      )
    }
  } else {
    # 需要全部基因有个rank值
    msigdbr_list <- split(msigdbr_df$gene, f = msigdbr_df$term)
    message("NOTE: check https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html for visualization")
    if (is.list(geneList)) {
      lapply(geneList, function(x) {
        fgsea::fgsea(msigdbr_list, x, ...)
      })
    } else {
      fgsea::fgsea(msigdbr_list, geneList, ...)
    }
  }
}

genList <- function(x) {
  if (is.null(names(x))) {
    y <- seq_along(x)
    names(y) <- x
  } else {
    y <- x
  }
  y
}

# gcap.plotEnrichment <- function() {
#
# }