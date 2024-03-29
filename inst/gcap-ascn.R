#!/usr/bin/env Rscript

suppressMessages(library(GetoptLong))
GetoptLong.options(help_style = "two-column")
VERSION = as.character(packageVersion("gcap"))

model = "XGB11"
genome = "hg38"
tightness = 1L
gapCN = 3L
outdir = getwd()

GetoptLong(
  "input=s", paste(
    "A file storing following columns: chromosome, start, end, total_cn,",
    "minor_cn, sample,",
    "(optinal) purity, (optinal) ploidy, (optinal) age, (optinal) gender, (optinal) type."
  ),
  "outdir=s", "Result output path.",
  "genome=s",  "Genome build version, should be hg38 or hg19.",
  "model=s", "Trained model name, should be one of XGB11, XGB32, XGB56.",
  "tightness=i", "Control the tightness to be a circular amplicon. If the value is larger, it is more likely a fCNA assigned to 'noncircular' instead of 'circular'.",
  "gcapCN=i", "A gene with copy number above background (ploidy + gapCN in general) would be treated as focal amplicon. Smaller, more amplicons.",
  "onlyOncogenes", "Only known oncogenes are kept for circular prediction."
)

suppressMessages(library(data.table))
suppressMessages(library(gcap))

input = fread(input, header = TRUE)
if (is.null(input$purity)) input$purity = 1L

gcap.ASCNworkflow(
  input,
  genome_build = genome,
  model = model,
  tightness = tightness,
  gap_cn = gapCN,
  only_oncogenes = onlyOncogenes,
  outdir = outdir
)
