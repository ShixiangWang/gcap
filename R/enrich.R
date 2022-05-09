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

gcap.enrich <- function(data,
                        target = c(NA, "circular", "noncircular"),
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

  # ref_file <- system.file(
  #   "extdata", paste0(genome_build, "_target_genes.rds"),
  #   package = "gcap", mustWork = TRUE
  # )
  # y <- readRDS(ref_file)$gene_id

  target <- match.arg(target)
  if (inherits(data, "fCNA")) {
    if (analysis_func == "gsea") message("")
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # The gene rank is useless currently
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    data <- data.table::copy(data$data)
    data[, amplicon_type := set_default_factor(amplicon_type)]
    dt_summary <- data[, .(
      cn = mean(total_cn, na.rm = TRUE),
      prob = mean(prob, na.rm = TRUE),
      N = .N
    ), by = .(gene_id, amplicon_type)][!is.na(gene_id)]

    geneList <- lapply(split(dt_summary, f = dt_summary$amplicon_type, drop = TRUE), function(dt) {
      if (dt$amplicon_type[1] == "circular") {
        rv <- dt[order(prob, N, cn, decreasing = TRUE)]
      } else {
        rv <- dt[order(cn, N)]
      }
      r <- seq_len(nrow(rv))
      names(r) <- rv$gene_id
      r
    })
    if (!is.na(target)) {
      geneList <- geneList[[target]]
    }
  } else {
    # A list/vector of (ranked) gene list
    geneList <- data
  }

  if (analysis_func == "enricher") {
    if (is.list(geneList)) {
      clusterProfiler::compareCluster(
        geneCluster = lapply(geneList, names),
        fun = clusterProfiler::enricher,
        TERM2GENE = msigdbr_df,
        #universe = y,
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        ...
      )
    } else {
      clusterProfiler::enricher(
        gene = names(geneList), TERM2GENE = msigdbr_df, 
        #universe = y,
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

# gcap.plotEnrichment <- function() {
#
# }