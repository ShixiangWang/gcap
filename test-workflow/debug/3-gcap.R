# remotes::install_github("ShixiangWang/ascat@v3-for-gcap-v1", subdir = "ASCAT")
# remotes::install_github("ShixiangWang/gcap")
# install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz", repos = NULL)

library(gcap)

# id为PRJEB42904，wes_395LC是tumor，id：ERR5242993，wes_395N是normal，id：ERR5243012

# hg38 ----------------
gcap.workflow(
  tumourseqfile = "~/share/gcap_debug/bam/ERR5242993.bam", 
  normalseqfile = "~/share/gcap_debug/bam/ERR5243012.bam",
  tumourname = "wes_395LC",
  normalname = "wes_395N",
  jobname = "wes_395",
  outdir = "~/share/gcap_debug/gcap_result",
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter", 
  g1000allelesprefix = file.path(
    "~/share/gcap_reference/1000G_loci_hg38/",
    "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
  ), 
  g1000lociprefix = file.path("~/share/gcap_reference/1000G_loci_hg38/",
                              "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
  ),
  GCcontentfile = "~/share/gcap_reference/GC_correction_hg38.txt",
  replictimingfile = "~/share/gcap_reference/RT_correction_hg38.txt",
  skip_finished_ASCAT = TRUE,
  skip_ascat_call = FALSE,
  result_file_prefix = "wes_395",
  genome_build = "hg38",
  model = "XGB11"
)

