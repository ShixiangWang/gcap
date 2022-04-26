gcap.plotCircos  <- function(fCNA,
highlight_genes = NULL,
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
  if (!startsWith(fCNA$data$gene_id[1], "ENSG")) {
    message("detected you have transformed ENSEMBL ID, also transforming gene annotation data")
    opts <- getOption("IDConverter.datapath", default = system.file("extdata", package = "IDConverter"))
    options(IDConverter.datapath = opts)
    target_genes$gene_id = IDConverter::convert_hm_genes(target_genes$gene_id, genome_build = genome_build)
    target_genes = target_genes[!is.na(target_genes$gene_id), ]
  }
  data_bed = merge(fCNA$data, target_genes, by = "gene_id", all.x = TRUE)
  data_bed$amplicon_type = ifelse(data_bed$amplicon_type %in% c("circular", "possibly_circular"), "circular", "noncircular")
  data_bed = data_bed[!is.na(data_bed$gene_id) & data_bed$chrom %in% chrs, c("chrom", "start", "end", "total_cn", "gene_id", "amplicon_type"), with = FALSE]
  colnames(data_bed)[1] = "chr"

  bed_list = split(data_bed, data_bed$amplicon_type)
  bed_cn = data.table::dcast(data_bed[, .(total_cn = mean(total_cn, na.rm = T)),
                                      by = .(chr, start, end, gene_id, amplicon_type)],
                             chr + start + end + gene_id ~ amplicon_type, value.var = "total_cn")
  bed_cn[, circular := ifelse(is.na(circular), 0, circular)]
  bed_cn[, noncircular := ifelse(is.na(noncircular), 0, -noncircular)]

  gap_after <- c(rep(2, length(chrs) - 1), 10)
  circlize::circos.par(
    "start.degree" = 90,
    "track.height" = 0.1,
    "gap.after" = gap_after)  # 10% of the circle radius

  circlize::circos.initializeWithIdeogram(
    species = genome_build,
    chromosome.index = chrs,
    track.height = circlize::convert_height(
      track_height,
      "mm"
    ), ideogram.height = circlize::convert_height(
      ideogram_height,
      "mm"
    ),
    plotType = if (!is.null(highlight_genes)) NULL else c("ideogram", "axis", "labels")
  )
  on.exit(circlize::circos.clear())
  
  if (!is.null(highlight_genes)) {
    bed = unique(data_bed[data_bed$gene_id %in% highlight_genes, c("chr", "start", "end", "gene_id")])
    circlize::circos.genomicLabels(bed, labels = bed$gene_id, side = "outside", cex = 0.5,
      col = as.numeric(factor(bed$gene_id)), line_col = as.numeric(factor(bed$gene_id))
      )
    circlize::circos.genomicIdeogram(species = genome_build,
    track.height = circlize::convert_height(
      track_height,
      "mm"
    ))
  }

  #circlize::circos.genomicLabels(zz[circular > 0], labels = zz[circular > 0]$gene_id, side = "outside", niceFacing = TRUE)
  circlize::circos.genomicDensity(bed_list[["circular"]], col = c("#FF000080"), track.height = 0.15)
  draw_y_axis(chrs)
  # draw_track_label(1, "circular freq")
  circlize::circos.genomicDensity(bed_list[["noncircular"]], col = c("#0000FF80"), track.height = 0.15)
  draw_y_axis(chrs)
  # draw_track_label(2, "noncircular freq")

  circlize::circos.genomicTrack(bed_cn, track.height = if (!is.null(highlight_genes)) 0.2 else 0.3, #ylim = c(0, 0.1),
                                panel.fun = function(region,
                                                     value, ...) {
    circlize::circos.genomicRect(region, value, ytop = value$circular,
                                 ybottom = value$noncircular,
                                 col = ifelse(value$circular > 0, "red", "blue"), border = NA)
    circlize::circos.lines(circlize::CELL_META$cell.xlim, c(0, 0), lty = 2, col = "#000040")
  })
  draw_y_axis(chrs)
  # draw_track_label(3, "average CN")
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

# draw_track_label = function(track_index, track_label, sector.index="chr1") {
#   circlize::circos.text(
#     circlize::get.cell.meta.data("cell.xlim")-mean(circlize::get.cell.meta.data("cell.xlim"))/2,
#     circlize::get.cell.meta.data("cell.ylim")-max(circlize::get.cell.meta.data("cell.ylim"))/2,
#     sector.index="chr1", track.index = track_index,
#     labels = track_label,facing = "clockwise", 
#     niceFacing = TRUE, adj = c(0,0), cex = 0.5)
# }
