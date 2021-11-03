# GCAP: Gene-level Circular Amplicon Prediction

<!-- badges: start -->
<!-- badges: end -->

The goal of **gcap** is to provide a end-to-end pipeline for predicting
circular amplicon (ecDNA) in gene level on WES data.

The input of **gcap** is bam files or allele specific copy number data and
output is a data table with features and prediction result.

## Installation

### Download reference files

For advanced users, you can prepare the reference files by following the instructions
from <https://github.com/shixiangwang/ascat/tree/v3.0>.

We recommend all users directly download the reference files from the links below:

- [Curated reference files for GCAP (WES)](https://zenodo.org/record/5533065)

> The prediction model is built with data on the top of hg38 genome build, so
hg38 based bam files as input is more recommended.

### Install softwares

[**alleleCount**](https://github.com/cancerit/alleleCount) is required to run **ASCAT** on WES bam data,
if you haven't installed [**conda**](https://docs.conda.io/en/latest/) or [**miniconda**](https://docs.conda.io/en/latest/miniconda.html), please install firstly,
then install the **alleleCount** in terminal with:

```bash
conda create -n cancerit -c bioconda cancerit-allelecount
```

> NOTE: **gcap** set the default **alleleCount** as the `~/miniconda3/envs/cancerit/bin/alleleCounter`,
if you use **conda** or other approaches, please set the path when you use corresponding functions.

Install **ASCAT** v3.0 (modified and adapted for GCAP workflow in HPC) in R console from GitHub with:

```r
# this is a fork version ASCAT
remotes::install_github("ShixiangWang/ascat@v3-for-gcap-v1", subdir = "ASCAT")
```

Install **gcap** in R console from GitHub with:

```r
remotes::install_github("ShixiangWang/gcap")
```

## Example

### Pipeline

- for one tumor-normal pair, you can refer to [one-pair.R](https://github.com/ShixiangWang/gcap/blob/master/test-workflow/one-pair.R).
- for multiple tumor-normal pairs, you can refer to [two-pair.R](https://github.com/ShixiangWang/gcap/blob/master/test-workflow/two-pairs.R).

To run **gcap** from bam files, a machine with **at least 80GB RAM** is required for
the `allelecount` process. If you set multiple threads, please note the parallel
computation is used in part of the workflow. You should balance the `nthread` setting
and the computing power your machine provides by yourself.

It generally takes `~0.5h` to finish one case (tumor-normal pair).

In our practice, when we want to process multiple cases, set `nthread = 22` and
directly let **gcap** handle multiple cases (instead of writing a loop yourself) is
good enough.

### Functions

For more custom and advanced control of the analysis, you can read the well
organized function list at [*package site*](https://shixiangwang.github.io/gcap/reference/index.html).


## LICENSE

This software and associated documentation files (the “Software”) are protected by copyright. This Software is provided “as is” (at your own risk) for internal non-commercial academic research purposes only. Please read the [Non-Commercial Academic License](LICENSE) in detail before downloading a copy. By installing or using this Software, you agree to be bound by the terms and conditions of the Non-Commercial Academic License.

All commercial use of the Software or any modification, manipulation or derivative of the Software, including but not limited to transfer, sale or licence to a commercial third party or use on behalf of a commercial third party (including but not limited to use as part of a service supplied to any third party for financial reward) is strictly prohibited and requires a commercial use licence. For further information please email <wangsx1@sysucc.org.cn> or <zhaoqi@sysucc.org.cn>.
