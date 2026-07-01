#' Validate BAM files, index files, and chromosome naming before running ASCAT
#'
#' This helper checks common causes of ASCAT/alleleCounter failures:
#' - BAM file existence and readability
#' - BAI index file existence and freshness (stale index warning)
#' - Chromosome naming format in BAM headers
#' - Basic BAM format validity (via samtools if available)
#'
#' @param bam_files character vector of BAM file paths to validate.
#' @param chrom_names expected chromosome names (e.g., 1:22).
#' @param genome_build genome build string ("hg38" or "hg19").
#'
#' @return invisibly returns a list of diagnostic info. Warnings/errors are
#' raised via the logger.
#' @keywords internal
validate_bam_inputs <- function(bam_files, chrom_names = 1:22, genome_build = "hg38") {
  lg <- set_logger()
  results <- list()

  for (bf in stats::na.omit(bam_files)) {
    lg$info("Validating BAM file: {bf}")

    # 1. BAM existence and readability
    if (!file.exists(bf)) {
      lg$fatal("BAM file does not exist: {bf}")
      stop("BAM file does not exist: ", bf)
    }
    if (file.access(bf, mode = 4) != 0) {
      lg$fatal("BAM file is not readable: {bf}")
      stop("BAM file is not readable: ", bf)
    }
    bam_size <- file.info(bf)$size
    if (is.na(bam_size) || bam_size == 0) {
      lg$fatal("BAM file is empty: {bf}")
      stop("BAM file is empty: ", bf)
    }
    lg$info("  BAM file exists, readable, size={bam_size} bytes")

    # 2. BAI index file check
    bai_file <- paste0(bf, ".bai")
    if (!file.exists(bai_file)) {
      lg$warn("BAI index file NOT FOUND: {bai_file}")
      lg$warn("  This will cause alleleCounter to fail with 'Null iterator' or 'Error scanning through bam'")
      lg$warn("  Generate the index with: samtools index {bf}")
      results[[bf]] <- c(results[[bf]], bai_missing = TRUE)
    } else {
      bai_size <- file.info(bai_file)$size
      bam_mtime <- file.info(bf)$mtime
      bai_mtime <- file.info(bai_file)$mtime

      if (is.na(bai_mtime) || is.na(bam_mtime)) {
        lg$warn("  Could not determine file modification times for BAM/BAI")
      } else if (bai_mtime < bam_mtime) {
        lg$warn("  BAI index file is OLDER than BAM file!")
        lg$warn("    BAM modified: {bam_mtime}")
        lg$warn("    BAI modified: {bai_mtime}")
        lg$warn("    ASCAT may produce '[W::hts_idx_load3] The index file is older than the data file' warnings")
        lg$warn("    If alleleCounter crashes, regenerate index: samtools index {bf}")
        results[[bf]] <- c(results[[bf]], bai_stale = TRUE)
      } else {
        lg$info("  BAI index file found and up-to-date")
      }
      results[[bf]] <- c(results[[bf]], bai_missing = FALSE)
    }

    # 3. Chromosome naming detection (using samtools if available)
    samtools_bin <- Sys.which("samtools")
    if (nzchar(samtools_bin)) {
      lg$info("  Detecting chromosome naming from BAM header via samtools...")
      header <- tryCatch(
        system2(samtools_bin, c("view", "-H", shQuote(bf)), stdout = TRUE, stderr = FALSE),
        error = function(e) NULL
      )
      if (!is.null(header) && length(header) > 0) {
        sq_lines <- grep("^@SQ", header, value = TRUE)
        snames <- gsub(".*SN:", "", sq_lines)
        snames <- gsub("\t.*", "", snames)

        if (length(snames) > 0) {
          has_chr_prefix <- any(grepl("^chr", snames))
          expected_prefix <- if (genome_build %in% c("hg38", "hg19")) "chr" else ""

          if (expected_prefix == "chr" && !has_chr_prefix) {
            lg$warn("  BAM uses chromosome names WITHOUT 'chr' prefix (e.g., '1', '2'),")
            lg$warn("  but genome_build='{genome_build}' expects 'chr' prefix (e.g., 'chr1', 'chr2').")
            lg$warn("  This mismatch can cause alleleCounter to fail silently or return no data.")
            lg$warn("  Consider re-aligning with chr-prefixed reference, or adjust reference files accordingly.")
            results[[bf]] <- c(results[[bf]], chr_mismatch = TRUE, chr_has_prefix = FALSE)
          } else if (expected_prefix == "chr" && has_chr_prefix) {
            lg$info("  BAM chromosome names have 'chr' prefix (compatible)")
            results[[bf]] <- c(results[[bf]], chr_mismatch = FALSE, chr_has_prefix = TRUE)
          } else {
            lg$info("  BAM chromosome names detected: {paste(utils::head(snames, 5), collapse=', ')}...")
          }

          results[[bf]] <- c(results[[bf]], chr_count = length(snames))

          # Show sample of chromosome names for debugging
          lg$info("  Found {length(snames)} sequences, examples: {paste(utils::head(snames, 5), collapse=', ')}")
        }
      } else {
        lg$info("  Could not read BAM header with samtools")
      }
    } else {
      lg$info("  samtools not found in PATH — skipping chromosome name detection")
      lg$info("  NOTE: to validate BAM chromosome naming, install samtools")
    }

    # 4. Quick BAM format check
    if (nzchar(samtools_bin)) {
      check_ok <- system2(samtools_bin, c("quickcheck", shQuote(bf)), stdout = FALSE, stderr = FALSE)
      if (check_ok != 0) {
        lg$warn("  samtools quickcheck reported issue with BAM file: {bf}")
        lg$warn("  The BAM may be truncated or corrupted — this can cause alleleCounter to crash")
        results[[bf]] <- c(results[[bf]], bam_corrupted = TRUE)
      } else {
        lg$info("  samtools quickcheck passed")
        results[[bf]] <- c(results[[bf]], bam_corrupted = FALSE)
      }
    }
    lg$info("")
  }

  invisible(results)
}
