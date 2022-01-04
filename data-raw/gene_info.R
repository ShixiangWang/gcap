hg38 <- readRDS("data-raw/hg38_gene_info.rds")
hg19 <- readRDS("data-raw/hg19_gene_info.rds")


process_geneinfo <- function(ref_genes) {
  ref_genes <- ref_genes[chrom %in% paste0("chr", 1:22) & gene_type == "protein_coding"]
  ref_genes$strand <- NULL
  ref_genes[, gene_id := gsub("(\\..+)", "", gene_id)]
  ref_genes$gene_type <- NULL
  ref_genes$gene_name <- NULL

  ref_genes
}

hg38 <- process_geneinfo(hg38)
hg19 <- process_geneinfo(hg19)

saveRDS(hg38, file = "inst/extdata/hg38_target_genes.rds")
saveRDS(hg19, file = "inst/extdata/hg19_target_genes.rds")

freq_df <- readRDS("~/Downloads/model_gene_amplicon_freq.rds")


library(tidyverse)

fts_freq <- freq_df %>%
  tidyr::pivot_wider(names_from = "type", values_from = "freq", values_fill = 0)
# %>%
#   dplyr::full_join(df_genes %>% dplyr::select(gene_id), by = "gene_id") %>%
#   dplyr::mutate_if(is.numeric, ~ifelse(is.na(.), 0, .))
colnames(fts_freq)[-1] <- paste0("freq_", colnames(fts_freq)[-1])
colnames(fts_freq)[5] <- "freq_HR"

setDT(fts_freq)

fts_freq
saveRDS(fts_freq, file = "inst/extdata/amplicon_freq.rds")
