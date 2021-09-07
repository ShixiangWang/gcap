gcap.workflow <- function(
  tumourseqfile, normalseqfile,
  tumourname, normalname, jobname = tumourname,
  extra_info,
  include_type = FALSE,
  genome_build = c("hg38", "hg19"),
  feature = c("with_type", "without_type"),
  target = c("circle", "nonLinear"),
  use_toy = FALSE,
  outdir = getwd(),
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter",
  g1000allelesprefix = file.path(
    "~/data/snp/1000G_loci_hg38",
    "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
  ),
  g1000lociprefix = file.path(
    "~/data/snp/1000G_loci_hg38",
    "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
  ),
  GCcontentfile = "~/data/snp/GC_correction_hg38.txt",
  replictimingfile = "~/data/snp/RT_correction_hg38.txt",
  nthreads = 22,
  minCounts = 10,
  BED_file = NA,
  probloci_file = NA,
  chrom_names = 1:22,
  gender = "XX",
  min_base_qual = 20,
  min_map_qual = 35,
  penalty = 70
) {
  genome_build <- match.arg(genome_build)
  # support loopping
  feature <- match.arg(feature, several.ok = TRUE)
  target <- match.arg(target, several.ok = TRUE)

  lg <- set_logger()
  lg$info("=====================")
  lg$info("   GCAP WORKFLOW")
  lg$info("=====================")
  lg$info()

  lg$info("=====================")
  lg$info("   GCAP WORKFLOW")
  lg$info("=====================")

  gcap.runASCAT(
    tumourseqfile,
    normalseqfile,
    tumourname,
    normalname,
    jobname,
    outdir = getwd(),
    allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter",
    g1000allelesprefix = file.path(
      "~/data/snp/1000G_loci_hg38",
      "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
    ),
    g1000lociprefix = file.path(
      "~/data/snp/1000G_loci_hg38",
      "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
    ),
    GCcontentfile = "~/data/snp/GC_correction_hg38.txt",
    replictimingfile = "~/data/snp/RT_correction_hg38.txt",
    nthreads = 22,
    minCounts = 10,
    BED_file = NA,
    probloci_file = NA,
    chrom_names = 1:22,
    gender = "XX",
    min_base_qual = 20,
    min_map_qual = 35,
    penalty = 70
  )

}
