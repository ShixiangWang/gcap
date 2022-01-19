rmsk <- data.table::fread("data-raw/rmsk_hg38.bed", header = FALSE)
rmsk[, V2 := V2 + 1L]

rmsk_gene <- collapse_to_genes(rmsk, drop = FALSE)
rmsk_gene

rmsk_gene[, intersect_ratio := intersect_size / (i.end - i.start + 1)]
summary(rmsk_gene$intersect_ratio)
table(rmsk_gene$V5)

rmsk_gene <- rmsk_gene[intersect_ratio >= 1]
rmsk_gene_sum <- rmsk_gene[, .N, by = .(gene_id, V4)]
saveRDS(rmsk_gene, file = "/Volumes/Extra/ec-materials/rmsk_gene_test.rds")

rmsk_gene_sum$gene_id <- factor(rmsk_gene_sum$gene_id, levels = unique(substr(hg38_gene_info$gene_id, 1, 15)))
data.table::dcast(rmsk_gene_sum, gene_id ~ V4, fill = 0L, value.var = "N", drop = FALSE)
