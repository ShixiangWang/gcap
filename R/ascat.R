# ref:
#   #1: https://github.com/tlesluyes/ascat/blob/v3.0/ASCAT/R/ascat.prepareHTS.R
#   #2: https://github.com/VanLoo-lab/ascat/blob/master/ExampleData/ASCAT_examplePipeline.R
# remotes::install_github("ShixiangWang/ascat@v3.0", subdir = "ASCAT")
# data <- na.omit(readr::read_csv("~/proj/ecDNA/data/tcga_wes_pair_info.csv"))
# raw_data_dir <- "~/data/gdc/wes"

#' Run ASCAT on tumor-normal pair WES data files
#'
#' A wrapper calling ASCAT on WES data on one or more tumor-normal pairs.
#' Note, for multiple tumor-normal pairs, the first 5 arguments should
#' be a vector with same length.
#'
#' @inheritParams ASCAT::ascat.prepareHTS
#' @inheritParams ASCAT::ascat.aspcf
#' @inheritParams ASCAT::ascat.GCcorrect
#' @param jobname job name, typically an unique name for a tumor-normal pair.
#' @param outdir result output path.
#' @param gender a vector of gender for each cases ("XX" or "XY").
#' Default = all female ("XX"). Ignore this if you don't include sex
#' chromosomes.
#'
#' @importFrom ASCAT ascat.prepareHTS ascat.loadData ascat.GCcorrect
#' ascat.plotRawData ascat.aspcf ascat.plotSegmentedData
#' ascat.runAscat
#'
#' @return Nothing. Check the `outdir` for results.
#' @export
gcap.runASCAT <- function(tumourseqfile, normalseqfile,
                          tumourname, normalname, jobname = tumourname,
                          outdir = getwd(),
                          allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter",
                          g1000allelesprefix = file.path(
                            "~/data/snp/1000G_loci_hg38",
                            "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
                          ),
                          g1000lociprefix = file.path(
                            "~/data/snp/1000G_loci_hg38",
                            "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
                          ),
                          GCcontentfile = "~/data/snp/GC_correction_hg38.txt",
                          replictimingfile = "~/data/snp/RT_correction_hg38.txt",
                          nthreads = 22,
                          minCounts = 10,
                          BED_file = NA,
                          probloci_file = NA,
                          chrom_names = 1:22,
                          gender = "XX",
                          min_base_qual = 20,
                          min_map_qual = 35,
                          penalty = 70) {
  cwd <- getwd()
  setwd(outdir)
  on.exit(setwd(cwd))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  lg <- set_logger()
  lg$info("> Run ASCAT for WES data <")
  lg$info()
  lg$info("Configs:")
  lg$info("  result path set to {outdir}")
  lg$info("  allelecounter_exe set to {allelecounter_exe}")
  lg$info("  g1000allelesprefix set to {g1000allelesprefix}")
  lg$info("  g1000lociprefix set to {g1000lociprefix}")
  lg$info("  GCcontentfile set to {GCcontentfile}")
  lg$info("  replictimingfile set to {replictimingfile}")
  lg$info("  nthreads set to {nthreads}")
  lg$info("  minCounts set to {minCounts}")
  lg$info("  BED_file set to {BED_file}")
  lg$info("  probloci_file set to {probloci_file}")
  lg$info("  chrom_names set to <{paste(chrom_names, collapse = ' ')}>")
  lg$info("  min_base_qual set to {min_base_qual}")
  lg$info("  penalty set to {penalty}")

  stopifnot(
    file.exists(allelecounter_exe), file.exists(GCcontentfile),
    file.exists(replictimingfile)
  )

  lg$info("{length(jobname)} jobs detected")

  run_one <- function(i) {
    tfile <- tumourseqfile[i]
    nfile <- normalseqfile[i]
    tn <- tumourname[i]
    nn <- normalname[i]
    id <- jobname[i]

    lg$info("start submitting job {id}")
    lg$info("     tumor data file: {tfile}")
    lg$info("    normal data file: {nfile}")
    lg$info("   tumor sample name: {tn}")
    lg$info("  normal sample name: {nn}")

    ascat.prepareHTS(
      tumourseqfile = tfile,
      normalseqfile = nfile,
      tumourname = tn,
      normalname = nn,
      allelecounter_exe = allelecounter_exe,
      g1000allelesprefix = g1000allelesprefix,
      g1000lociprefix = g1000lociprefix,
      nthreads = nthreads,
      minCounts = minCounts,
      BED_file = BED_file,
      probloci_file = probloci_file,
      chrom_names = chrom_names,
      min_base_qual = min_base_qual,
      min_map_qual = min_map_qual,
      skip_allele_counting_tumour = FALSE,
      skip_allele_counting_normal = FALSE
    )

    ascat.bc <- ascat.loadData(
      paste0(tn, "_tumourLogR.txt"), paste0(tn, "_tumourBAF.txt"),
      paste0(tn, "_normalLogR.txt"), paste0(tn, "_normalBAF.txt"),
      chrs = chrom_names,
      gender = gender
    )
    ascat.bc <- ascat.GCcorrect(ascat.bc, GCcontentfile, replictimingfile)
    ascat.plotRawData(ascat.bc)
    ascat.bc <- ascat.aspcf(ascat.bc, penalty = penalty)
    ascat.plotSegmentedData(ascat.bc)
    ascat.output <- ascat.runAscat(ascat.bc, gamma = 1L, pdfPlot = TRUE)
    saveRDS(ascat.output, file = paste0(id, ".ASCAT.rds"))

    lg$info("job {id} done")
  }

  lapply(seq_along(jobname), run_one)
  lg$info("analysis done, check {outdir} for results")
}
