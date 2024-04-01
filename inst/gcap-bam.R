#!/usr/bin/env Rscript

suppressMessages(library(GetoptLong))
GetoptLong.options(help_style = "two-column")
VERSION = as.character(packageVersion("gcap"))

model = "XGB11"
genome = "hg38"
tightness = 1L
gapCN = 3L
outdir = getwd()
extra = NULL
bed = NULL

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
    "A file storing following columns: sample, age, gender and (optional) type"
  ),
  "bed=s", paste(
    "A BED file for only looking at SNPs within specific intervals"
  ),
  "outdir=s", "Result output path.",
  "genome=s",  "Genome build version, should be hg38 or hg19.",
  "model=s", "Trained model name, should be one of XGB11, XGB32, XGB56.",
  "tightness=i", "Control the tightness to be a circular amplicon. If the value is larger, it is more likely a fCNA assigned to 'noncircular' instead of 'circular'.",
  "gapCN=i", "A gene with copy number above background (ploidy + gapCN in general) would be treated as focal amplicon. Smaller, more amplicons.",
  "onlyOncogenes", "Only known oncogenes are kept for circular prediction.",
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
if (!is.null(extra)) extra = fread(extra, header = TRUE)
if (is.null(bed)) bed = NA

gcap.workflow(
  tumourseqfile = input$tumorfile,
  normalseqfile = input$normalfile,
  tumourname = input$tumorname,
  normalname = input$normalname,
  jobname = if (is.null(input$sample)) input$tumorname else input$sample,
  allelecounter_exe = allelecounter,
  g1000allelesprefix = g1000allelesprefix,
  g1000lociprefix = g1000lociprefix,
  GCcontentfile = GCcontentfile,
  replictimingfile = replictimingfile,
  extra_info = extra,
  include_type = "type" %in% colnames(extra),
  genome_build = genome,
  model = model,
  tightness = tightness,
  gap_cn = gapCN,
  only_oncogenes = onlyOncogenes,
  outdir = outdir,
  nthreads = nthreads,
  BED_file = bed,
  result_file_prefix = "Result",
  skip_finished_ASCAT = skip
)
