# https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt
hom = data.table::fread("data-raw/HOM_MouseHumanSequence.rpt.txt")
hom

table(hom$`NCBI Taxon ID`)

hom = hom[, c(1, 3, 4)]
hom
hom_n = hom[, .N, by = .(`DB Class Key`)]
hom_n

hom = hom[`DB Class Key` %in% hom_n[N == 2]$`DB Class Key`]
hom

mouse = readRDS("data-raw/mouse_mm10_gene_info.rds")
nchar(mouse$gene_id[1]) - 2

mouse[, gene_id := substr(gene_id, 1, 18)]
mouse

human = readRDS("data-raw/human_hg38_gene_info.rds")
human[, gene_id := substr(gene_id, 1, 15)]
human

head(mouse)
head(human)

mouse = mouse[, c(1, 2, 3, 5, 6)]
mouse
human = human[, c(5, 6)]
human

hom
hom2 = data.table::dcast(hom,  `DB Class Key` ~ `NCBI Taxon ID`, value.var = "Symbol")
colnames(hom2) = c("key", "human", "mouse")
hom2

hom2_m1 = merge(hom2, human, by.x = "human", by.y = "gene_name", all.x = FALSE, all.y = FALSE)
hom2_m2 = merge(hom2_m1, mouse, by.x = "mouse", by.y = "gene_name", all.x = FALSE, all.y = FALSE)
hom2_m2

amplicon_freq_mm10 = data.table::copy(amplicon_freq)
amplicon_freq_mm10[, gene_id := IDConverter::convert_custom(gene_id, from = "gene_id.x", to = "gene_id.y", dt = hom2_m2)]
amplicon_freq_mm10 = amplicon_freq_mm10[!is.na(gene_id)]
amplicon_freq_mm10

library(data.table)
saveRDS(amplicon_freq_mm10, file = "inst/extdata/amplicon_freq_mm10.rds")

amplicon_freq
hg38_target_genes

mouse2 = readRDS("data-raw/mouse_mm10_gene_info.rds")
mouse2
mm10_target_genes = mouse2[gene_type == "protein_coding" & chrom %in% paste0("chr", 1:19), .(chrom, start, end, gene_id)]
mm10_target_genes[, gene_id := substr(gene_id, 1, 18)]
mm10_target_genes
hg38_target_genes

saveRDS(mm10_target_genes, file = "inst/extdata/mm10_target_genes.rds")

somatic_gene_cn_mm10 = data.table::copy(somatic_gene_cn)
somatic_gene_cn_mm10[, gene_id := IDConverter::convert_custom(gene_id, from = "gene_id.x", to = "gene_id.y", dt = hom2_m2)]
somatic_gene_cn_mm10 = somatic_gene_cn_mm10[!is.na(gene_id)]
somatic_gene_cn_mm10

saveRDS(somatic_gene_cn_mm10, file = "inst/extdata/somatic_gene_cn_mm10.rds")

saveRDS(hom2_m2, file = "data-raw/human-mouse-id-map.rds")
