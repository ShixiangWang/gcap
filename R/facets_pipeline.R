# A file designed for mouse cancer genome analysis workflow
#
# cd /data3/wsx/R/x86_64-pc-linux-gnu-library/4.2/facets/extcode/
# g++ -std=c++11 -I/data3/wsx/miniconda3/envs/circlemap/include snp-pileup.cpp -L/data3/wsx/miniconda3/envs/circlemap/lib -lhts -Wl,-rpath=/data3/wsx/miniconda3/envs/circlemap/lib -o snp-pileup
gcap.workflow.facets <- function(tumourseqfile, normalseqfile,
                                 jobname,
                                 extra_info = NULL,
                                 include_type = FALSE,
                                 genome_build = c("mm10", "hg38", "hg19"),
                                 model = "XGB11",
                                 tightness = 1L,
                                 gap_cn = 3L,
                                 overlap = 1,
                                 pro_cval = 100,
                                 only_oncogenes = FALSE,
                                 snp_file = "path/to/genome_build_responding.vcf.gz",
                                 outdir = getwd(),
                                 result_file_prefix = paste0("gcap_", uuid::UUIDgenerate(TRUE)),
                                 util_exe = system.file("extcode", "snp-pileup", package = "facets"),
                                 nthreads = 1,
                                 skip_finished_facets = TRUE,
                                 skip_facets_call = FALSE) {
  stopifnot(sum(duplicated(jobname), na.rm = TRUE) == 0L, sum(is.na(jobname)) == 0, .Platform$OS.type == "unix")
  genome_build <- match.arg(genome_build)
  check_model(model)

  if (!requireNamespace("facets")) {
    stop("Package 'facets' cannot be found, please install it firstly.")
  }

  lg <- set_logger()
  lg$info("=========================================")
  lg$info("   GCAP WORKFLOW backend by FACETS")
  lg$info("=========================================")
  lg$info("")

  lg$info("=====================")
  lg$info("Step 1: Run FACETS")
  lg$info("=====================")

  if (is.character(extra_info) && file.exists(extra_info)) {
    extra_info <- data.table::fread(extra_info, header = TRUE)
  }
  if (!is.null(extra_info)) {
    lg$info("Check extra_info values")
    if (!all(jobname %in% extra_info$sample)) {
      lg$warn("Cannot find all 'jobname' identifiers in 'sample' column of extra_info")
      lg$info("Try subsetting")
      extra_info <- extra_info[extra_info$sample %in% jobname, ]
      if (nrow(extra_info) == 0L) {
        lg$fatal("'sample' column is consistent with 'jobname', if no extra data, set it to NULL!")
        stop("Bad input")
      } else {
        if (nrow(extra_info) != length(jobname)) {
          lg$warn("Non-consistent sample records detected, try filling with default values")
          extra_info <- merge(
            data.table::data.table(sample = jobname),
            extra_info,
            by = "sample", all.x = TRUE
          )
          if (is.character(extra_info$gender)) {
            extra_info$gender <- data.table::fifelse(is.na(extra_info$gender), "XX", extra_info$gender, na = "XX")
          } else if (is.numeric(extra_info$gender)) {
            extra_info$gender <- data.table::fifelse(is.na(extra_info$gender), 0, extra_info$gender, na = 0)
          }
        }
      }
    } else {
      # Make sure the useful rows kept
      extra_info <- extra_info[extra_info$sample %in% jobname, ]
    }

    lg$info("Check sample order consistence in jobname and extra_info$sample")
    if (!all(jobname == extra_info$sample)) {
      lg$info("Sample order not consistent, try reordering")
      extra_info$sample <- factor(extra_info$sample, levels = jobname)
      extra_info <- extra_info[order(extra_info$sample), ]
      extra_info$sample <- as.character(extra_info$sample)
    }
  }

  totalT <- parallel::detectCores()
  if (nthreads > totalT) {
    lg$info("{nthreads} threads are set but only {totalT} cores are available, reset it to {totalT}")
    nthreads <- totalT
  }

  if (!skip_facets_call) {
    if (!file.exists(util_exe)) {
      stop("snp-pileup not found, please compile it by following https://github.com/mskcc/facets/blob/master/inst/extcode/README.txt")
    }

    if (!file.exists(snp_file)) {
      stop("snp file for genome not found, please download one based your genome_build from NCBI, e.g., https://ftp.ncbi.nlm.nih.gov/snp/organisms/archive/mouse_10090/VCF/")
    }

    lg$info("creating pileup files...")
    facets_dir = file.path(outdir, "facets")
    if (!dir.exists(facets_dir)) dir.create(facets_dir, recursive = TRUE)

    run_one <- function(i) {
      tfile <- tumourseqfile[i]
      nfile <- normalseqfile[i]
      id <- jobname[i]

      tryCatch(
        {
          if (!all(file.exists(tfile, nfile))) {
            lg$fatal("Not all bam files exist")
          }
          facets_file = file.path(facets_dir, paste0(id, ".facets.gz"))


          lg$info("Genrating snp-pileup facets.gz file from linux shell...")
          cmd1 = sprintf("%s -g -q15 -Q20 -P100 -r20,0 %s %s %s %s",
                         util_exe, snp_file, facets_file, nfile, tfile)
          system(cmd1)

          if (!(file.exists(facets_file) && file.size(facets_file) > 200)) {
            lg$fatal("file {facets_file} not be properly generated")
          }

          lg$info("Running FACETS standard analysis...")
          #if (is.null(Sys.getenv("VROOM_CONNECTION_SIZE"))) Sys.setenv("VROOM_CONNECTION_SIZE" = 28311552 * 6)
          #library("pctGCdata")
          #library("facets")
          eval(parse(text = "library('pctGCdata')"))
          set.seed(1234)
          rcmat = facets::readSnpMatrix(facets_file)
          xx = facets::preProcSample(rcmat, gbuild = genome_build)
          oo = facets::procSample(xx, cval = pro_cval)
          fit = facets::emcncf(oo)

          #plot
          pdf(file.path(outdir, paste0(id, "_facets.pdf")))
          facets::plotSample(x=oo, emfit=fit)
          facets::logRlogORspider(oo$out, oo$dipLogR)
          while (!is.null(dev.list()))  dev.off()

          # output purity and ploidy -----
          purity = fit$purity
          purity = round(purity, 2)
          ploidy = fit$ploidy
          ploidy = round(ploidy, 1)

          data = fit$cncf[, c("chrom", "start", "end", "tcn.em", "lcn.em")]
          colnames(data) = c("chromosome", "start", "end", "total_cn", "minor_cn")
          data$sample = id
          data = subset(data, subset = data$chromosome %in% if (startsWith(genome_build, "mm")) 1:19 else 1:22)
          data$chromosome = paste0("chr", data$chromosome)
          # NOTE: I don't know why there exists data with end < start
          dataList = list(data = data.table::as.data.table(data)[end > start], pp = data.table::data.table(
            sample = id,
            purity = purity,
            ploidy = ploidy
          ))
          saveRDS(dataList, file = file.path(outdir, paste0(id, "_facets.rds")))

          lg$info("job {id} done")
        },
        error = function(e) {
          lg$fatal("job {id} failed in FACETS due to following error")
          lg$info(e$message)
          lg$info("=====")
          lg$info("Please check your input bam files (if missing bam index? if its alignment quality is lower?)")
          lg$info("=====")
        }
      )
    }

    jobname2 = jobname
    if (skip_finished_facets) {
      drop_idx <- sapply(file.path(outdir, paste0(jobname2, "_facets.rds")), file.exists)
      if (sum(drop_idx) > 0) {
        # ALL input parameter vectors need to be update!!
        tumourseqfile <- tumourseqfile[!drop_idx]
        normalseqfile <- normalseqfile[!drop_idx]
        jobname2 <- jobname2[!drop_idx]

        lg$info("{sum(drop_idx)} job(s) skipped, {sum(!drop_idx)} to run.")
      } else {
        lg$info("No job to skip.")
      }
    }

    if (length(jobname2) > 0) parallel::mclapply(seq_along(jobname2), run_one, mc.cores = nthreads)
    lg$info("FACETS analysis done, check {outdir} for results")

  }

  lg$info("checking FACETS result files")
  FACETS_files <- file.path(outdir, paste0(jobname, "_facets.rds"))
  keep_idx <- sapply(FACETS_files, function(x) {
    keep <- if (file.exists(x)) {
      file.info(x)$size > 200
    } else {
      lg$warn("result file {x} does not exist, the corresponding FACETS calling has error occurred ")
      FALSE
    }
    if (!keep) {
      lg$warn("{x} contains a failed FACETS job, will discard it before next step")
    }
    keep
  })
  FACETS_files <- FACETS_files[keep_idx]
  if (length(FACETS_files) < 1) {
    lg$fatal("no sucessful FACETS result file to proceed!")
    lg$fatal("check your FACETS setting before make sure this case could not be used!")
    stop()
  }

  lg$info("============================================================")
  lg$info("Step 2: Extract features and collapse features to gene level")
  lg$info("============================================================")

  model_input <- gcap.runBuildflow(
    FACETS_files,
    extra_info,
    include_type = include_type,
    genome_build = genome_build,
    overlap = overlap
  )

  lg$info("=======================")
  lg$info("Step 3: Run prediction")
  lg$info("=======================")
  model_input$prob <- gcap.runPrediction(
    model_input,
    model = model
  )

  lg$info("====================================")
  lg$info("Step 4: Run scoring and summarizing")
  lg$info("====================================")
  out <- gcap.runScoring(model_input, genome_build,
                         tightness = tightness, gap_cn = gap_cn,
                         only_oncogenes = only_oncogenes)

  save_file <- file.path(outdir, paste0(result_file_prefix, "_prediction_result.rds"))
  lg$info("Saving raw prediction result to {save_file}")
  saveRDS(out$data, file = save_file)
  fCNA <- out$fCNA
  rm(out)
  invisible(gc())

  save_file <- file.path(outdir, paste0(result_file_prefix, c("_fCNA_records.csv", "_sample_info.csv")))
  lg$info("Saving fCNA records and sample info to {paste(save_file, collapse = ', ')}")
  fCNA$saveToFiles(outdir, result_file_prefix)

  lg$info("=======================================")
  lg$info(" Done! Thanks for using GCAP workflow")
  lg$info("=======================================")

  invisible(fCNA)
}

read_copynumber_facets = function(x) {
  stopifnot(length(x) >= 1)
  rv = list()
  rv$data = purrr::map_df(x, ~readRDS(.)$data)
  rv$purity_ploidy = purrr::map_df(x, ~readRDS(.)$pp)

  return(rv)
}
