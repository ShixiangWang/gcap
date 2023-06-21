tumourseqfile <- "/data3/wsx/share/fq_mouse/fp_WES/bam/ERR424934.bam"
normalseqfile <- "/data3/wsx/share/fq_mouse/fp_WES/bam/ERR424937.bam"
jobname <- "ERR424934"

gcap.workflow.seqz(tumourseqfile, normalseqfile, jobname,
  genome_build = "mm10", ref_file = "/data3/wsx/aa_data_repo/mm10/mm10.fa",
  data_tmp_dir = "~/gcap_data", outdir = "~/gcap_data",
  util_exe = "~/miniconda3/bin/sequenza-utils",
  samtools_exe = "~/miniconda3/envs/circlemap/bin/samtools",
  tabix_exe = "~/miniconda3/envs/circlemap/bin/tabix"
)

gcap.workflow.seqz(tumourseqfile, normalseqfile, jobname,
  genome_build = "mm10", ref_file = "/data3/wsx/aa_data_repo/mm10/mm10.fa",
  data_tmp_dir = "~/gcap_data", outdir = "~/gcap_data",
  util_exe = "~/miniconda3/bin/sequenza-utils",
  samtools_exe = "~/miniconda3/envs/circlemap/bin/samtools",
  tabix_exe = "~/miniconda3/envs/circlemap/bin/tabix",
  only_oncogenes = TRUE
)


gcap.workflow.facets(
  tumourseqfile, normalseqfile, jobname,
  genome_build = "mm10", snp_file = "~/data/mouse-all.snp.vcf.gz",
  outdir = "~/gcap_data"
)

gcap.workflow.facets(
  tumourseqfile, normalseqfile, jobname,
  genome_build = "mm10", snp_file = "~/data/mouse-all.snp.vcf.gz",
  outdir = "~/gcap_data", only_oncogenes = TRUE
)
