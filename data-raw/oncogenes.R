## code to prepare `oncogenes` dataset goes here
oncogenes <- data.table::fread("http://ongene.bioinfo-minzhao.org/ongene_human.txt")
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

oncogenes[, gene_id := IDConverter::convert_hm_genes(oncogenes$OncogeneName, type = "symbol")]

oncogenes[!is.na(gene_id)]


usethis::use_data(oncogenes, overwrite = TRUE)

# mouse
map <- readRDS("data-raw/human-mouse-id-map.rds")
library(IDConverter)
oncogenes_mouse <- convert_custom(unique(oncogenes$gene_id), from = "gene_id.x", to = "gene_id.y", dt = map)
oncogenes_mouse <- unique(na.omit(oncogenes_mouse))
nrow(oncogenes)

saveRDS(oncogenes_mouse, file = "inst/extdata/oncogenes_mouse.rds")
