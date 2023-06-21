#!/usr/bin/env Rscript

suppressMessages(library(GetoptLong))
VERSION = as.character(packageVersion("gcap"))
Citation = "Wang et al. Machine learning-based extrachromosomal DNA identification in large-scale cohorts reveals its clinical implications in cancer."


subCommands(
  help_head = qqcat("gcap (v@{VERSION})\n"),
  help_foot = qqcat("----------\nCitation:\n  @{Citation}\nURL:\n  https://github.com/ShixiangWang/gcap"),
  "bam", "gcap-bam.R",
         "Run GCAP workflow with tumor-normal paired BAM files",
  "ascn", "gcap-ascn.R",
          "Run GCAP workflow with curated allele-specific copy number data"
)
