# ref:
#   #1: https://github.com/tlesluyes/ascat/blob/v3.0/ASCAT/R/ascat.prepareHTS.R
#   #2: https://github.com/VanLoo-lab/ascat/blob/master/ExampleData/ASCAT_examplePipeline.R
# remotes::install_github("ShixiangWang/ascat@v3.0", subdir = "ASCAT")
# data <- na.omit(readr::read_csv("~/proj/ecDNA/data/tcga_wes_pair_info.csv"))
# raw_data_dir <- "~/data/gdc/wes"

#' Run ASCAT on tumor-normal pair WES data files
#'
#' A wrapper calling ASCAT on WES data on one or more tumor(-normal paired) bam data.
#' Note, for multiple tumor-normal pairs, the first 5 arguments should
#' be a vector with same length.
#'
#' @inheritParams ASCAT::ascat.prepareHTS
#' @inheritParams ASCAT::ascat.aspcf
#' @inheritParams ASCAT::ascat.GCcorrect
#' @param tumourseqfile Full path to the tumour BAM file.
#' @param normalseqfile Full path to the normal BAM file.
#' @param jobname job name, typically an unique name for a tumor-normal pair.
#' @param outdir result output path.
#' @param gender a vector of gender for each cases ("XX" or "XY").
#' Default = all female ("XX"). Ignore this if you don't include sex
#' chromosomes.
#' @param chrom_names A vector containing the names of chromosomes to be
#' considered (optional, default=1:22).
#' @param skip_finished_ASCAT if `TRUE`, skipped finished ASCAT calls
#' to save time.
#'
#' @importFrom ASCAT ascat.prepareHTS ascat.loadData ascat.GCcorrect
#' ascat.plotRawData ascat.aspcf ascat.plotSegmentedData
#' ascat.runAscat
#'
#' @return Nothing. Check the `outdir` for results.
#' @export
gcap.runASCAT <- function(tumourseqfile,
                          normalseqfile = NA_character_,
                          tumourname,
                          normalname = NA_character_, jobname = tumourname,
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
                          penalty = 70,
                          skip_finished_ASCAT = FALSE) {
  stopifnot(all(gender %in% c("XX", "XY")))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  cwd <- getwd()
  setwd(outdir)
  on.exit(setwd(cwd))
  if (length(gender) == 1 && length(tumourseqfile) > 1) {
    gender <- rep(gender, length(tumourseqfile))
  }
  if (length(normalseqfile) == 1 && length(tumourseqfile) > 1) {
    normalseqfile <- rep(normalseqfile, length(tumourseqfile))
    normalname <- rep(normalname, length(tumourseqfile))s
  }

  lg <- set_logger()
  lg$info("> Run ASCAT on WES data <")
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
  lg$info("  gender set to <{paste(gender, collapse = ' ')}>")
  lg$info("  min_base_qual set to {min_base_qual}")
  lg$info("  min_map_qual set to {min_map_qual}")
  lg$info("  penalty set to {penalty}")
  lg$info("  skip_finished_ASCAT set to {skip_finished_ASCAT}")

  stopifnot(
    file.exists(allelecounter_exe), file.exists(GCcontentfile),
    file.exists(replictimingfile),
    length(gender) == length(tumourseqfile)
  )

  lg$info("{length(jobname)} jobs detected")

  run_one <- function(i) {
    tfile <- tumourseqfile[i]
    nfile <- normalseqfile[i]
    tn <- tumourname[i]
    nn <- normalname[i]
    id <- jobname[i]
    gd <- gender[i]

    lg$info("start submitting job {id}")
    lg$info("     tumor data file: {tfile}")
    lg$info("    normal data file: {nfile}")
    lg$info("   tumor sample name: {tn}")
    lg$info("  normal sample name: {nn}")

    if (is.na(normalseqfile) || is.na(normalname)) {
      lg$info("run with tumor only model, this is only supported in specified ASCAT version")
      skip_norm = TRUE
    } else skip_norm = FALSEs

    # In some special cases, ASCAT failed after alleleCounter.
    # Maybe we should handle the ASCAT source code to fix corresponding issues
    tryCatch(
      {
        if (!all(file.exists(tfile, nfile))) {
          lg$fatal("Not all bam files exist")
        }
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
          gender = gd,
          chrom_names = chrom_names,
          min_base_qual = min_base_qual,
          min_map_qual = min_map_qual,
          skip_allele_counting_tumour = FALSE,
          skip_allele_counting_normal = skip_norm
        )

        ascat.bc <- ascat.loadData(
          paste0(tn, "_tumourLogR.txt"), paste0(tn, "_tumourBAF.txt"),
          if (skip_norm) NULL else paste0(tn, "_normalLogR.txt"),
          if (skip_norm) NULL else paste0(tn, "_normalBAF.txt"),
          chrs = chrom_names,
          gender = gender
        ) # New parameter genomeVersion in the latest version of ASCAT
        ascat.bc <- ascat.GCcorrect(ascat.bc, GCcontentfile, replictimingfile)
        ascat.plotRawData(ascat.bc)
        if (skip_norm) {
          # https://github.com/VanLoo-lab/ascat/issues/73
          # gg = ascat.predictGermlineGenotypes(ascat.bc, platform = "AffySNP6")
          gg = ascat.predictGermlineGenotypes_custom(
                ascat.bc, 
                maxHomozygous = 0.05, 
                proportionHetero = 0.59, 
                proportionHomo = 0.38, 
                proportionOpen = 0.02, 
                segmentLength = 100)
        } else gg = NULL
        ascat.bc <- ascat.aspcf(ascat.bc, ascat.gg = gg, penalty = penalty)
        ascat.plotSegmentedData(ascat.bc)
        ascat.output <- ascat.runAscat(ascat.bc, gamma = 1L, pdfPlot = TRUE)
        saveRDS(ascat.output, file = paste0(id, ".ASCAT.rds"))

        lg$info("job {id} done")
      },
      error = function(e) {
        lg$fatal("job {id} failed in ASCAT due to following error")
        lg$info(e$message)
        lg$info("=====")
        lg$info("Please check your input bam files (if missing bam index? if its alignment quality is lower?)")
        lg$info("=====")
      }
    )
  }

  if (skip_finished_ASCAT) {
    drop_idx <- sapply(paste0(jobname, ".ASCAT.rds"), file.exists)
    if (sum(drop_idx) > 0) {
      # ALL input parameter vectors need to be update!!
      tumourseqfile <- tumourseqfile[!drop_idx]
      normalseqfile <- normalseqfile[!drop_idx]
      jobname <- jobname[!drop_idx]
      tumourname <- tumourname[!drop_idx]
      normalname <- normalname[!drop_idx]
      gender <- gender[!drop_idx]

      lg$info("{sum(drop_idx)} ASCAT job(s) skipped, {sum(!drop_idx)} to run.")
    } else {
      lg$info("No ASCAT job to skip.")
    }
  }

  if (length(jobname) > 0) lapply(seq_along(jobname), run_one)
  lg$info("ASCAT analysis done, check {outdir} for results")
}
