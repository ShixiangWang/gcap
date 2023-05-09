#' GCAP sequenza workflow for gene-level amplicon prediction
#'
#' @inheritParams gcap.workflow
#' @param genome_build genome build version, should be one of 'hg38', 'hg19' and 'mm10'.
#' @param ref_file a reference genome file, should be consistent with `genome_build` option.
#' @param data_tmp_dir a directory path for storing temp data for reuse in handling multiple samples.
#' @param util_exe the path to `sequenza-utils`.
#' @param samtools_exe the path to `samtools_exe`.
#' @param tabix_exe the path to `tabix`.
#' @param skip_finished_sequenza if `TRUE`, skip finished sequenza runs.
#' @param skip_sequenza_call if `TRUE`, skip calling sequenza.
#' This is useful when you have done this step and just want
#' to run next steps.
#' @return a list of invisible `data.table` and corresponding files saved to local machine.
#' @export
gcap.workflow.seqz <- function(tumourseqfile, normalseqfile,
                               jobname,
                               extra_info = NULL,
                               include_type = FALSE,
                               genome_build = c("mm10", "hg38", "hg19"),
                               model = "XGB11",
                               tightness = 1L,
                               gap_cn = 3L,
                               overlap = 1,
                               only_oncogenes = FALSE,
                               ref_file = "path/to/reference.fa",
                               data_tmp_dir = "~/gcap_data",
                               outdir = getwd(),
                               result_file_prefix = paste0("gcap_", uuid::UUIDgenerate(TRUE)),
                               util_exe = "~/miniconda3/bin/sequenza-utils",
                               samtools_exe = "~/miniconda3/bin/samtools",
                               tabix_exe = "~/miniconda3/bin/tabix",
                               nthreads = 1,
                               skip_finished_sequenza = TRUE,
                               skip_sequenza_call = FALSE) {
  stopifnot(sum(duplicated(jobname), na.rm = TRUE) == 0L, sum(is.na(jobname)) == 0, .Platform$OS.type == "unix")
  genome_build <- match.arg(genome_build)
  check_model(model)

  if (!requireNamespace("sequenza")) {
    stop("Package 'sequenza' cannot be found, please install it firstly.")
  }

  lg <- set_logger()
  lg$info("=========================================")
  lg$info("   GCAP WORKFLOW backend by Sequenza")
  lg$info("=========================================")
  lg$info("")

  lg$info("=====================")
  lg$info("Step 1: Run Sequenza")
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

  if (!skip_sequenza_call) {
    gc_file = file.path(data_tmp_dir, paste0(genome_build, ".gc50Base.wig.gz"))
    if (!dir.exists(data_tmp_dir)) {
      dir.create(data_tmp_dir, recursive = TRUE)
      lg$info("created data directory for storing reference relavant annotation")
    }
    if (!file.exists(gc_file)) {
      lg$info("creating gc50 data...")
      if (!file.exists(util_exe)) {
        stop("Sequenza utils not found, please install it by `pip install sequenza-utils`")
      }
      cmd = sprintf("%s gc_wiggle -w 50 --fasta %s -o %s", util_exe, ref_file, gc_file)
      system(cmd)
      lg$info("done")
    } else {
      lg$info("gc50 data already exists, we will directly use it, please make sure it is complete and correct")
    }

    lg$info("creating sequenza files...")
    seqz_dir = file.path(outdir, "seqz")
    sseqz_dir = file.path(outdir, "small-seqz")
    if (!dir.exists(seqz_dir)) dir.create(seqz_dir, recursive = TRUE)
    if (!dir.exists(sseqz_dir)) dir.create(sseqz_dir, recursive = TRUE)

    run_one <- function(i) {
      tfile <- tumourseqfile[i]
      nfile <- normalseqfile[i]
      id <- jobname[i]

      tryCatch(
        {
          if (!all(file.exists(tfile, nfile))) {
            lg$fatal("Not all bam files exist")
            return(NULL)
          }
          seqz_file = file.path(seqz_dir, paste0(id, ".seqz.gz"))
          sseqz_file = file.path(sseqz_dir, paste0(id, ".small.seqz.gz"))
          sseqz_file2 = file.path(sseqz_dir, paste0(id, ".small_filtered.seqz.gz"))


          lg$info("Genrating proper small seqz.gz file from linux shell...")
          cmd1 = sprintf("%s bam2seqz -n %s -t %s --fasta %s -gc %s -S %s -T %s -o %s",
                         util_exe, nfile, tfile, ref_file, gc_file, samtools_exe, tabix_exe, seqz_file)
          lg$info(cmd1)
          system(cmd1)

          cmd2 = sprintf("%s seqz_binning --seqz %s -w 50 -T %s -o %s",
                         util_exe, seqz_file, tabix_exe, sseqz_file)
          lg$info(cmd2)
          system(cmd2)

          cmd3 = sprintf('zcat %s | awk "/^chromosome|chr[12]?[0-9XY]\t/{print}" | gzip > %s',
                         sseqz_file, sseqz_file2)
          lg$info(cmd3)
          system(cmd3)

          if (!(file.exists(sseqz_file2) && file.size(sseqz_file2) > 200)) {
            lg$fatal("file {sseqz_file2} not be properly generated")
            return(NULL)
          }

          lg$info("Running sequenza standard analysis...")
          if (is.null(Sys.getenv("VROOM_CONNECTION_SIZE"))) Sys.setenv("VROOM_CONNECTION_SIZE" = 28311552 * 6)
          s1 = sequenza::sequenza.extract(sseqz_file2, assembly = genome_build)
          s2 = sequenza::sequenza.fit(s1, chromosome.list = if (startsWith(genome_build, "mm")) 1:19 else 1:22)
          sequenza::sequenza.results(
            sequenza.extract = s1,
            cp.table = s2,
            sample.id = id,
            out.dir = outdir,
            CNt.max = 2000,
            chromosome.list = if (startsWith(genome_build, "mm")) 1:19 else 1:22
          )

          lg$info("job {id} done")
        },
        error = function(e) {
          lg$fatal("job {id} failed in sequenza due to following error")
          lg$info(e$message)
          lg$info("=====")
          lg$info("Please check your input bam files (if missing bam index? if its alignment quality is lower?)")
          lg$info("=====")
        }
      )
    }

    jobname2 = jobname
    if (skip_finished_sequenza) {
      drop_idx <- sapply(file.path(outdir, paste0(jobname2, "_segments.txt")), file.exists)
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
    lg$info("sequenza analysis done, check {outdir} for results")

  }

  lg$info("checking sequenza result files")
  sequenza_files <- file.path(outdir, paste0(jobname, "_segments.txt"))
  keep_idx <- sapply(sequenza_files, function(x) {
    keep <- if (file.exists(x)) {
      file.info(x)$size > 200
    } else {
      lg$warn("result file {x} does not exist, the corresponding sequenza calling has error occurred ")
      FALSE
    }
    if (!keep) {
      lg$warn("{x} contains a failed sequenza job, will discard it before next step")
    }
    keep
  })
  sequenza_files <- sequenza_files[keep_idx]
  if (length(sequenza_files) < 1) {
    lg$fatal("no sucessful sequenza result file to proceed!")
    lg$fatal("check your sequenza setting before make sure this case could not be used!")
    stop()
  }

  lg$info("============================================================")
  lg$info("Step 2: Extract features and collapse features to gene level")
  lg$info("============================================================")

  model_input <- gcap.runBuildflow(
    sequenza_files,
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

read_copynumber_seqz_cn = function(x) {
  stopifnot(length(x) >= 1)
  SAMPLE = sub("_segments.txt", "", basename(x))
  res <- purrr::map2_df(x, SAMPLE, function(x, y) {
    df <- data.table::fread(x)
    df <- df[, c("chromosome", "start.pos",
                 "end.pos", "CNt", "B"), with = FALSE]
    df$sample = y
    colnames(df) <- c("Chromosome", "Start.bp", "End.bp",
                      "modal_cn", "minor_cn", "sample")
    df
  })

  return(res)
}

read_copynumber_seqz_pp = function(x) {
  stopifnot(length(x) >= 1)
  x2 = sub("_segments.txt", "_alternative_solutions.txt", x)
  SAMPLE = sub("_segments.txt", "", basename(x))
  res <- purrr::map2_df(x2, SAMPLE, function(x, y) {
    df <- data.table::fread(x)[1, 1:2]
    df$sample = y
    colnames(df) <- c("purity", "ploidy", "sample")
    df[, c(3, 1, 2)]
  })

  return(res)
}
