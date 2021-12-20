#!/usr/bin/env Rscript

suppressMessages(library(GetoptLong))
GetoptLong.options(help_style = "two-column")
VERSION = as.character(packageVersion("gcap"))

model = "XGB32"
genome = "hg38"
outdir = getwd()
extra = NULL

allelecounter = "~/miniconda3/envs/cancerit/bin/alleleCounter"
g1000allelesprefix = file.path("~/data/snp/1000G_loci_hg38",
                               "1kg.phase3.v5a_GRCh38nounref_allele_index_chr")
g1000lociprefix = file.path("~/data/snp/1000G_loci_hg38",
                            "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr")
GCcontentfile = "~/data/snp/GC_correction_hg38.txt"
replictimingfile = "~/data/snp/RT_correction_hg38.txt"
nthreads = min(22, parallel::detectCores())
minCounts = 10

GetoptLong(
  "input=s", paste(
    "A file storing following columns: tumorfile, normalfile, tumorname, normalname,",
    "sample (optional, it should be identical to sample in 'extra')"
  ),
  "extra=s", paste(
    "A file storing following columns: sample, age, gender and type (optional)"
  ),
  "outdir=s", "Result output path.",
  "genome=s",  "Genome build version, should be hg38 or hg19.",
  "model=s", "Model name, should be one of XGB11, XGB32, XGB54.",
  "allelecounter=s", "Path to the allele counter executable.",
  "g1000allelesprefix=s", "Prefix path to the 1000 Genomes alleles reference files.",
  "g1000lociprefix=s", "Prefix path to the 1000 Genomes SNP reference files.",
  "GCcontentfile=s", "File containing the GC content around every SNP for increasing window sizes",
  "replictimingfile=s", "File containing replication timing at every SNP for various cell lines (optional)",
  "nthreads=i", "The number of parallel processes for getting allele counts.",
  "minCounts=i", "Minimum depth required in the normal for a SNP to be considered.",
  "skip", "If skipping finished ASCAT calls to save time."
)

suppressMessages(library(data.table))
suppressMessages(library(gcap))

input = fread(input, header = TRUE)
extra = fread(extra, header = TRUE)

gcap.workflow(
  tumourseqfile = input$tumorfile,
  normalseqfile = input$normalfile,
  tumourname = input$tumorname,
  normalname = input$normalname,
  jobname = if (is.null(input$sample)) input$tumorname else input$sample,
  extra_info = extra,
  include_type = "type" %in% colnames(extra),
  genome_build = genome,
  model = model,
  #target = "circle",
  outdir = outdir,
  nthreads = nthreads,
  skip_finished_ASCAT = skip
)
