## code to prepare `oncogenes` dataset goes here
oncogenes = data.table::fread("http://ongene.bioinfo-minzhao.org/ongene_human.txt")
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

oncogenes[, gene_id := IDConverter::convert_hm_genes(oncogenes$OncogeneName, type = "symbol")]

oncogenes[!is.na(gene_id)]


usethis::use_data(oncogenes, overwrite = TRUE)
