#!/usr/bin/env Rscript
# Test GCAP workflow on two tumor-normal pairs to
# check how it works on multiple cases
data <- structure(list(
  case = c("TCGA-VD-A8K9", "TCGA-25-2400"),
  pair_id_uniq = c("TCGA-VD-A8K9-01", "TCGA-25-2400-01"),
  tumor = c("TCGA-VD-A8K9-01", "TCGA-25-2400-01"), file_name_tumor = c(
    "C1760.TCGA-VD-A8K9-01A-11D-A39W-08.1_gdc_realn.bam",
    "C239.TCGA-25-2400-01A-01W-0799-08.2_gdc_realn.bam"
  ), normal = c(
    "TCGA-VD-A8K9-10",
    "TCGA-25-2400-10"
  ), file_name_normal = c(
    "C1760.TCGA-VD-A8K9-10A-01D-A39Z-08.2_gdc_realn.bam",
    "C239.TCGA-25-2400-10A-01W-0799-08.2_gdc_realn.bam"
  ), type = c(
    "Primary Tumor",
    "Primary Tumor"
  ), project = c("TCGA-UVM", "TCGA-OV"), center = c(
    "C1760",
    "C239"
  ), file_id_tumor = c(
    "964b182c-c84c-45c8-8af6-d81c5327a46c",
    "ae39a0af-4bd1-4e95-9996-a563a6f232d5"
  ), file_id_normal = c(
    "307c987b-87e6-4bc1-9756-2af84b992101",
    "c274af2a-a462-490e-8abb-22f12d73859b"
  )
), row.names = c(NA, -2L), class = c("tbl_df", "tbl", "data.frame"), na.action = structure(c(`2` = 2L), class = "omit"))

raw_data_dir <- "~/data/gdc/wes"
tfile <- file.path(raw_data_dir, data$file_name_tumor)
nfile <- file.path(raw_data_dir, data$file_name_normal)

tn <- data$tumor
nn <- data$normal
id <- data$pair_id_uniq

library(gcap)

gcap::gcap.workflow(
  tumourseqfile = tfile,
  normalseqfile = nfile,
  tumourname = tn,
  normalname = nn,
  jobname = id,
  extra_info = data.frame(
    sample = id,
    age = c(59, 60),
    gender = c(0, 1)
  ),
  target = c("circle", "nonLinear"),
  outdir = "~/proj/gcap/test-workflow/result",
  result_file_prefix = "test_two_cases",
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
  skip_finished_ASCAT = TRUE
)
