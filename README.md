# GCAP: Gene-level Circular Amplicon Prediction

<!-- badges: start -->
<!-- badges: end -->

In a nutshell, **gcap** provides end-to-end machine learning approache for predicting
circular amplicon (also known as ecDNA, extrachromosomal DNA ) in gene level from WES (tumor-normal paired BAM) data,
allele specific copy number data (e.g., results from [ASCAT](https://github.com/VanLoo-lab/ascat) or [Sequenza](https://cran.r-project.org/package=sequenza)), or even
absolute integer copy number data (e.g., results from [ABSOLUTE](https://software.broadinstitute.org/cancer/cga/absolute)). The former two data
sources are preferred as input of **gcap** .

## Installation

### Download reference files (WES bam data only)

For advanced users, you can prepare the reference files by following the instructions
from <https://github.com/shixiangwang/ascat/tree/v3.0>.

We recommend all users directly download the reference files from the links below:

- [Curated reference files for GCAP (WES)](https://zenodo.org/record/6364977)

> The prediction model was built with data on the top of hg38 genome build, while 
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
```

### Install GCAP (required)

Install **gcap** in R console from GitHub with:

```r
remotes::install_github("ShixiangWang/gcap")
```

## Example

### Pipeline (WES bam data only)

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

## Related tools

- [DoAbsolute](https://github.com/ShixiangWang/DoAbsolute): Automate Absolute Copy Number Calling using 'ABSOLUTE' package.
- [sigminer](https://github.com/ShixiangWang/sigminer): An easy-to-use and scalable toolkit for genomic alteration signature (a.k.a. mutational signature) analysis and visualization in R.

## Output     
 **gcap** outputs two data tables including feature table and prediction result.
 
## Citations 


## LICENSE

This software and associated documentation files (the “Software”) are protected by copyright. This Software is provided “as is” (at your own risk) for internal non-commercial academic research purposes only. Please read the [Non-Commercial Academic License](LICENSE) in detail before downloading a copy. By installing or using this Software, you agree to be bound by the terms and conditions of the Non-Commercial Academic License.

All commercial use of the Software or any modification, manipulation or derivative of the Software, including but not limited to transfer, sale or licence to a commercial third party or use on behalf of a commercial third party (including but not limited to use as part of a service supplied to any third party for financial reward) is strictly prohibited and requires a commercial use licence. For further information please email <wangsx1@sysucc.org.cn> or <zhaoqi@sysucc.org.cn>.
