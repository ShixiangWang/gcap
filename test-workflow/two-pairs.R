#!/usr/bin/env Rscript
# Test GCAP workflow on two tumor-normal pairs to
# check how it works on multiple cases

data <- na.omit(readr::read_csv("~/proj/ecDNA/data/tcga_wes_pair_info.csv"))

data <- head(data, 2)

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
    gender = c(0, 1),
    type = "COAD"
  ), # this is not true data
  include_type = TRUE,
  feature = "with_type",
  target = c("circle", "nonLinear"),
  use_toy = TRUE, # use builtin toy model for prediction
  outdir = "~/proj/gcap/test-workflow/result",
  result_file = "test_gcap_workflow_on_two_cases_with_type.csv",
  skip_finished_ASCAT = TRUE
)
