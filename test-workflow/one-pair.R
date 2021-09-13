#!/usr/bin/env Rscript
# Test GCAP workflow on one tumor-normal pair

data <- na.omit(readr::read_csv("~/proj/ecDNA/data/tcga_wes_pair_info.csv"))

data <- head(data, 1)

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
    gender = 1, # 1 for XY, 0 for XX. You can also use 'XY' and 'XX'
    type = "COAD"
  ), # this is not true data
  include_type = TRUE,
  feature = "with_type",
  target = c("circle", "nonLinear"),
  use_toy = TRUE, # use builtin toy model for prediction
  chrom_names = c(1:22, "X"),
  outdir = "~/proj/gcap/test-workflow/result",
  result_file = "test_gcap_workflow_on_one_case_with_type.csv",
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
