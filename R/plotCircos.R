gcap.plotCircos  <- function(fCNA,
                            genome_build = c("hg38", "hg19"),
                            chrs = paste0("chr", 1:22),
                            track_height = 0.2,
                            ideogram_height = 1,
                            ...) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Please install 'circlize' package firstly.")
  }
  genome_build = match.arg(genome_build)
  target_genes = readRDS(file.path(system.file("extdata", package = "gcap"),
                                   paste0(genome_build, "_target_genes.rds")))

  data_bed = merge(fCNA$data, target_genes, by = "gene_id", all.x = TRUE)
  data_bed = data_bed[!is.na(data_bed$gene_id), c("chrom", "start", "end", "total_cn", "gene_id", "amplicon_type"), with = FALSE]
  colnames(data_bed)[1] = "chr"

  bed_list = split(data_bed, data_bed$amplicon_type)
  bed_cn = data.table::dcast(data_bed[, .(total_cn = mean(total_cn, na.rm = T)),
                                      by = .(chr, start, end, gene_id, amplicon_type)],
                             chr + start + end + gene_id ~ amplicon_type, value.var = "total_cn")
  bed_cn[, circular := ifelse(is.na(circular), 0, circular)]
  bed_cn[, noncircular := ifelse(is.na(noncircular), 0, -noncircular)]

  gap_after <- c(rep(2, length(chrs) - 1), 10)
  circlize::circos.par("track.height" = 0.1,
                       "gap.after" = gap_after)  # 10% of the circle radius
  circlize::circos.initializeWithIdeogram(
    species = species,
    chromosome.index = chrs,
    track.height = circlize::convert_height(
      track_height,
      "mm"
    ), ideogram.height = circlize::convert_height(
      ideogram_height,
      "mm"
    )
  )
  on.exit(circlize::circos.clear())

  #circlize::circos.genomicLabels(zz[circular > 0], labels = zz[circular > 0]$gene_id, side = "outside", niceFacing = TRUE)
  circos.genomicDensity(bed_list[[1]], col = c("#FF000080"), track.height = 0.2)
  draw_y_axis(chrs)
  circos.genomicDensity(bed_list[[3]], col = c("#0000FF80"), track.height = 0.2)
  draw_y_axis(chrs)

  circlize::circos.genomicTrack(zz, track.height = 0.4, #ylim = c(0, 0.1),
                                panel.fun = function(region,
                                                     value, ...) {
    circlize::circos.genomicRect(region, value, ytop = value$circular,
                                 ybottom = value$noncircular,
                                 col = ifelse(value$circular > 0, "red", "blue"), border = NA)
    circlize::circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#000040")

  })
  draw_y_axis(chrs)
}

draw_y_axis = function(chrs) {
  circlize::circos.yaxis(
    labels.cex = 0.3, tick.length = circlize::convert_x(
      0.3,
      "mm", circlize::get.cell.meta.data("sector.index"),
      circlize::get.cell.meta.data("track.index")
    ),
    side = "right", sector.index = chrs[length(chrs)]
  )
}
