#!/usr/bin/env Rscript
# Test GCAP workflow on one tumor-normal pair

library(readr)

PROJ_DIR = "/data3/wsx_data/gcap-analysis/"
setwd(file.path(PROJ_DIR, "SYSUCC-GC"))
df <- read_csv("data/own-data-pairs.csv")[1:2, ]

tfile <- df$tumor_file
nfile <- df$normal_file

tn <- df$tumor_id
nn <- df$normal_id
id <- df$case

library(gcap)
gcap.workflow(
  tumourseqfile = tfile, normalseqfile = nfile, tumourname = tn, normalname = nn, jobname = id, extra_info = NULL,
  genome_build = "hg19", outdir = "/data3/wsx_data/test_ascat_latest", result_file_prefix = "test",
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter",
  alleles.prefix = file.path(
    "~/data/1000G_loci_hg19/",
    "1000genomesAlleles2012_chr"
  ),
  loci.prefix = file.path("~/data/1000G_loci_hg19/", "1000genomesloci2012chrstring_chr"),
  GCcontentfile = "~/data/GC_correction_hg19.txt",
  replictimingfile = "~/data/RT_correction_hg19.txt"
)
