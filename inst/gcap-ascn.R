#!/usr/bin/env Rscript

suppressMessages(library(GetoptLong))
GetoptLong.options(help_style = "two-column")
VERSION = as.character(packageVersion("gcap"))

model = "XGB32"
genome = "hg38"
outdir = getwd()
GetoptLong(
  "input=s", paste(
    "A file storing following columns: chromosome, start, end, total_cn,",
    "minor_cn, sample, purity, ploidy (optinal), age (optinal), gender (optinal), type (optinal)."
  ),
  "outdir=s", "Result output path.",
  "genome=s",  "Genome build version, should be hg38 or hg19.",
  "model=s", "Model name, should be one of XGB11, XGB32, XGB54."
)

suppressMessages(library(data.table))
suppressMessages(library(gcap))

gcap.ASCNworkflow(
  fread(input, header = TRUE),
  genome_build = genome,
  model = model,
  #target = "circle",
  outdir = outdir
  #result_file_prefix = paste0("gcap_", uuid::UUIDgenerate(TRUE))
)
