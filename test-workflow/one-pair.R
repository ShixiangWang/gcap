#!/usr/bin/env Rscript
# Test GCAP workflow on one tumor-normal pair

data <- structure(list(
  case = "TCGA-VD-A8K9", pair_id_uniq = "TCGA-VD-A8K9-01",
  tumor = "TCGA-VD-A8K9-01", file_name_tumor = "C1760.TCGA-VD-A8K9-01A-11D-A39W-08.1_gdc_realn.bam",
  normal = "TCGA-VD-A8K9-10", file_name_normal = "C1760.TCGA-VD-A8K9-10A-01D-A39Z-08.2_gdc_realn.bam",
  type = "Primary Tumor", project = "TCGA-UVM", center = "C1760",
  file_id_tumor = "964b182c-c84c-45c8-8af6-d81c5327a46c", file_id_normal = "307c987b-87e6-4bc1-9756-2af84b992101"
), row.names = c(
  NA,
  -1L
), class = c("tbl_df", "tbl", "data.frame"), na.action = structure(c(`2` = 2L), class = "omit"))

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
    age = 60,
    gender = 1 # 1 for XY, 0 for XX. You can also use 'XY' and 'XX'
  ),
  target = c("circle", "nonLinear"),
  chrom_names = c(1:22, "X"),
  outdir = "~/proj/gcap/test-workflow/result",
  result_file_prefix = "test_one_case",
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
