# GCAP: Gene-level Circular Amplicon Prediction

<!-- badges: start -->
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fshixiangwang%2Fgcap&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)
<!-- badges: end -->

In a nutshell, **gcap** provides an end-to-end workflow for predicting
circular amplicon (also known as ecDNA, extra-chromosomal DNA ) in gene level with machine learning approach,
then classifying cancer samples into different focal amplification (fCNA) types,
based on input from WES (tumor-normal paired BAM, with corresponding `.bai` index files) data,
allele specific copy number data (e.g., results from [ASCAT](https://github.com/VanLoo-lab/ascat) or [Sequenza](https://cran.r-project.org/package=sequenza)), or even
absolute integer copy number data (e.g., results from [ABSOLUTE](https://software.broadinstitute.org/cancer/cga/absolute)). The former two data
sources are preferred as input of **gcap** .

## Installation

### Download reference files (WES bam data only)

For advanced users, you can prepare the reference files by following the instructions
from <https://github.com/shixiangwang/ascat/tree/v3.0>.

We recommend all users directly download the reference files from the links below:

- [Curated reference files for GCAP (WES)](https://zenodo.org/records/6364977)

> The prediction model was built with data on the top of hg38 genome build, so 
hg38-based BAM file input is more recommended.

### Install alleleCount (WES bam data only)

[**alleleCount**](https://github.com/cancerit/alleleCount) is required to run **ASCAT** on WES bam data,
if you haven't installed [**conda**](https://docs.conda.io/en/latest/) or [**miniconda**](https://docs.conda.io/en/latest/miniconda.html), please install firstly,
then install the **alleleCount** in terminal with:

```bash
conda create -n cancerit -c bioconda cancerit-allelecount
```

> NOTE: **gcap** set the default **alleleCount** as the `~/miniconda3/envs/cancerit/bin/alleleCounter`,
if you use **conda** or other approaches, please set the path when you use corresponding functions.

### Install ASCAT (required)

Install **ASCAT** v3.0 (modified and adapted for GCAP workflow in HPC) in R console from GitHub with:

```r
# This is a forked version ASCAT
remotes::install_github("ShixiangWang/ascat@v3-for-gcap-v1", subdir = "ASCAT")
# A ASCAT version with loose SAM flag, useful sometimes
# remotes::install_github("ShixiangWang/ascat@v3-f1", subdir = "ASCAT")
# See https://github.com/ShixiangWang/gcap/issues/27
```

> Here we used a fixed version of ASCAT for the GCAP data pre-processing, if you want to adopted the latest
> updates in processing your data, please refer to <https://github.com/VanLoo-lab/ascat> for generating the required
> allele-specific copy number data, and refer to [`?gcap.ASCNworkflow()`](https://shixiangwang.github.io/gcap/reference/gcap.ASCNworkflow.html) for downstream analysis.

### Alternatives to ASCAT

For the latest version of GCAP, [sequenza](https://shixiangwang.github.io/gcap/reference/gcap.workflow.seqz.html) or [facets](https://shixiangwang.github.io/gcap/reference/gcap.workflow.facets.html) are supported for preprocessing the bam data, please refer to the provided links for usage.

### Install GCAP (required)

Install **gcap** in R console:

```r
# r-universe
install.packages('gcap', repos = c('https://shixiangwang.r-universe.dev', 'https://cloud.r-project.org'))

# or GitHub
remotes::install_github("ShixiangWang/gcap")
```

If you would like to use CLI program in Shell terminal, run the following code in your R console after installation:

```r
gcap::deploy()
```

Two scripts `gcap-bam.R` and `gcap-ascn.R` shall be linked to your path `/usr/local/bin/`.
You can use one of them based on you input data.

**NOTE**

For users with package **GetoptLong** version `>= 1.1.0`, a main command is implemented
and also linked to `/usr/local/bin/` when calling `deploy()`. So you can type `gcap` as
a unified interface.

```sh
$ gcap
gcap (v1.0.0)
Usage: gcap [command] [options]

Commands:
  bam     Run GCAP workflow with tumor-normal paired BAM files
  ascn    Run GCAP workflow with curated allele-specific copy number data

----------
Citation:
  GCAP
URL:
  https://github.com/ShixiangWang/gcap
```

NOTE: **gcap** use XGBOOST < 1.6, if you have installed a latest version,
you can install the specified version with:

```R
install.packages("https://cran.r-project.org/src/contrib/Archive/xgboost/xgboost_1.5.2.1.tar.gz", repos = NULL)
```

## Example

Run the following code to see a quick example:

```r
library(gcap)

data("ascn")
rv <- gcap.ASCNworkflow(ascn, outdir = tempdir(), model = "XGB11")
rv
```

### Pipeline (WES bam data only)

- for one tumor-normal pair, you can refer to [one-pair.R](https://github.com/ShixiangWang/gcap/blob/master/test-workflow/one-pair.R). [test-workflow/debug](test-workflow/debug) contains a full workflow for data obtained from SRA.
- for multiple tumor-normal pairs, you can refer to [two-pair.R](https://github.com/ShixiangWang/gcap/blob/master/test-workflow/two-pairs.R).

To run **gcap** from bam files, a machine with **at least 80GB RAM** is required for
the `allelecount` process. If you set multiple threads, please note the parallel
computation is used in part of the workflow. You should balance the `nthread` setting
and the computing power your machine provides by yourself.

It generally takes `~0.5h` to finish one case (tumor-normal pair).

In our practice, when we want to process multiple cases, set `nthread = 22` and
directly let **gcap** handle multiple cases (instead of writing a loop yourself) is
good enough.

A recommended setting for Slurm is given as:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -o output-%J.o
#SBATCH -n 22
#SBATCH --mem=102400
```

Templates of practical calling command with provided hg38 and hg19 annotations are given below:

```r
# hg38 ----------------
gcap.workflow(
  tumourseqfile = tfile, normalseqfile = nfile, tumourname = tn, normalname = nn, jobname = id,
  outdir = outdir,
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter", 
  g1000allelesprefix = file.path(
    "/data/wsx/data/1000G_loci_hg38/",
    "1kg.phase3.v5a_GRCh38nounref_allele_index_chr"
  ), 
  g1000lociprefix = file.path("/data/wsx/data/1000G_loci_hg38/",
                              "1kg.phase3.v5a_GRCh38nounref_loci_chrstring_chr"
  ),
  GCcontentfile = "/data/wsx/data/GC_correction_hg38.txt",
  replictimingfile = "/data/wsx/data/RT_correction_hg38.txt",
  skip_finished_ASCAT = TRUE,
  skip_ascat_call = FALSE,
  result_file_prefix = "xxx",
  extra_info = df,
  include_type = FALSE,
  genome_build = "hg38",
  model = "XGB11"
)

# hg19 ----------------
gcap.workflow(
  tumourseqfile = tfile, normalseqfile = nfile, tumourname = tn, normalname = nn, jobname = id,
  outdir = outdir,
  allelecounter_exe = "~/miniconda3/envs/cancerit/bin/alleleCounter", g1000allelesprefix = file.path(
    "/data/wsx/data/1000G_loci_hg19/",
    "1000genomesAlleles2012_chr"
  ), g1000lociprefix = file.path("/data/wsx/data/1000G_loci_hg19/", "1000genomesloci2012chrstring_chr"),
  GCcontentfile = "/data/wsx/data/GC_correction_hg19.txt", replictimingfile = "/data/wsx/data/RT_correction_hg19.txt",
  skip_finished_ASCAT = TRUE,
  skip_ascat_call = FALSE,
  result_file_prefix = "xxx",
  extra_info = NULL,
  include_type = FALSE,
  genome_build = "hg19",
  model = "XGB11"
)

```

### Pipeline (Allele specific or absolute copy number data)

Please refer to [`?gcap.ASCNworkflow()`](https://shixiangwang.github.io/gcap/reference/gcap.ASCNworkflow.html).

### Functions

For more custom and advanced control of the analysis, you can read the structured documentation at [*package site*](https://shixiangwang.github.io/gcap/reference/index.html).

## Logging

For better debugging and rechecking.
The logging information of your operation with **gcap** would be saved into
an independent file. You can use the following commands to get the file path
and print logging message. Please note you have to use `:::` to access these
functions as they are not exported from **gcap**.

```r
> gcap:::get_log_file()
[1] "~/Library/Logs/gcap/gcap.log"
> gcap:::cat_log_file()
```

## Docker image

A docker image is available in [ghcr](https://github.com/ShixiangWang/gcap/pkgs/container/gcap) along with its corresponding [Dockerfile](https://github.com/ShixiangWang/gcap/blob/master/Dockerfile). This image comes pre-installed with all the necessary software. However, users are responsible for mapping the required reference files and input data files on their own. The Dockerfile can be customized according to the user's specific requirements, as permitted by the license we provide.

```sh
docker pull ghcr.io/shixiangwang/gcap:latest
```

## Related tools

- [gcaputils](https://github.com/ShixiangWang/gcaputils): **gcap** utils for downstream analysis and visualization.  
- [DoAbsolute](https://github.com/ShixiangWang/DoAbsolute): Automate Absolute Copy Number Calling using 'ABSOLUTE' package.
- [sigminer](https://github.com/ShixiangWang/sigminer): An easy-to-use and scalable toolkit for genomic alteration signature (a.k.a. mutational signature) analysis and visualization in R.

## Output 

 **gcap** outputs two data tables including feature table and prediction result.

## Citations 

Wang, S., Wu, C. Y., He, M. M., Yong, J. X., Chen, Y. X., Qian, L. M., ... & Zhao, Q. (2024). Machine learning-based extrachromosomal DNA identification in large-scale cohorts reveals its clinical implications in cancer. Nature Communications, 15(1), 1-17. <https://doi.org/10.1038/s41467-024-45479-6>

## LICENSE

This software and associated documentation files (the “Software”) are protected by copyright. This Software is provided “as is” (at your own risk) for internal non-commercial academic research purposes only. Please read the [Non-Commercial Academic License](LICENSE) in detail before downloading a copy. By installing or using this Software, you agree to be bound by the terms and conditions of the Non-Commercial Academic License.

All commercial use of the Software or any modification, manipulation or derivative of the Software, including but not limited to transfer, sale or licence to a commercial third party or use on behalf of a commercial third party (including but not limited to use as part of a service supplied to any third party for financial reward) is strictly prohibited and requires a commercial use licence. **This software is protected by the P. R. China patent [202211067952.6](https://www.patentguru.com/cn/search?q=202211067952.6)** For further information please email <wangsx1@sysucc.org.cn> or <zhaoqi@sysucc.org.cn>.

