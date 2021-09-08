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
    gender = 1 # 1 for XY, 0 for XX. You can also use 'XY' and 'XX'
  ), # this is not true data
  include_type = FALSE,
  feature = "without_type",
  target = c("circle", "nonLinear"),
  use_toy = TRUE, # use builtin toy model for prediction
  chrom_names = c(1:22, "X"),
  outdir = "~/proj/gcap/test-workflow/result",
  result_file = "test_gcap_workflow_on_one_case_without_type.csv"
)

gcap::gcap.workflow(
  tumourseqfile = tfile, 
  normalseqfile = nfile,
  tumourname = tn, 
  normalname = nn, 
  jobname = id,
  extra_info = data.frame(
    sample = id,
    age = 60,
    gender = 1,
    type = "COAD"
  ), # this is not true data
  include_type = TRUE,
  feature = c("with_type", "without_type"),
  target = c("circle", "nonLinear"),
  use_toy = TRUE, # use builtin toy model for prediction
  outdir = "~/proj/gcap/test-workflow/result",
  result_file = "test_gcap_workflow_on_one_case_all.csv",
  skip_finished_ASCAT = TRUE
)