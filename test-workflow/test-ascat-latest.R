#!/usr/bin/env Rscript
# Test GCAP workflow on one tumor-normal pair

library(readr)

PROJ_DIR = "/data3/wsx_data/gcap-analysis/"
setwd(file.path(PROJ_DIR, "SYSUCC-GC"))
df <- read_csv("data/own-data-pairs.csv")[1, ]

tfile <- df$tumor_file
nfile <- df$normal_file

tn <- df$tumor_id
nn <- df$normal_id
id <- df$case

library(gcap)
gcap.workflow(
  tumourseqfile = tfile, normalseqfile = nfile, tumourname = tn, normalname = nn, jobname = id, extra_info = NULL,
  genome_build = "hg19", outdir = "/data3/wsx_data/test_ascat_latest", result_file_prefix = "test", chrom_names = 1:2,
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter",
  alleles.prefix = file.path(
    "~/data/1000G_loci_hg19/",
    "1000genomesAlleles2012_chr"
  ),
  loci.prefix = file.path("~/data/1000G_loci_hg19/", "1000genomesloci2012chrstring_chr"),
  GCcontentfile = "~/data/GC_correction_updated_hg19.txt",
  replictimingfile = "~/data/RT_correction_updated_hg19.txt"
)


# Update correction reference data for latest ASCAT -----------------------

# gc = data.table::fread("~/data/GC_correction_hg19.txt")
# rt = data.table::fread("~/data/RT_correction_hg19.txt")
#
# gc$V1 = paste0(gc$chr, "_", gc$Position)
# colnames(gc)[1] = "SNP_ID"
# data.table::fwrite(gc, file = "~/data/GC_correction_updated_hg19.txt", sep = "\t")
#
# rm(gc); gc()
# rt$SNP_ID = paste0(rt$chr, "_", rt$pos)
# rt_cols = colnames(rt)
# data.table::setcolorder(rt, neworder = c("SNP_ID", rt_cols[-length(rt_cols)]))
# rt
# data.table::fwrite(rt, file = "~/data/RT_correction_updated_hg19.txt", sep = "\t")

