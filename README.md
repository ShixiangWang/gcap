# gcap: Gene-level Circular Amplicon Prediction

<!-- badges: start -->
<!-- badges: end -->

The goal of **gcap** is to provide a end-to-end pipeline for predicting
circular amplicon (ecDNA) in gene level on WES data.

The input of **gcap** is bam files and output is a data table with features and
prediction result.

## Installation

### Download reference files

For advanced users, you can prepare the reference files by following the instructions
from <https://github.com/shixiangwang/ascat/tree/v3.0>.

We recommend all users directly download the reference files from the links below:

- [**for hg38 genome build**](https://ora.ox.ac.uk/objects/uuid:08e24957-7e76-438a-bd38-66c48008cf52) (the files are used and tested in our study)
- [for hg19 genome build](https://ora.ox.ac.uk/objects/uuid:2c1fec09-a504-49ab-9ce9-3f17bac531bc)

### Install softwares

[**alleleCount**](https://github.com/cancerit/alleleCount) is required to run **ASCAT** on WES bam data,
if you haven't installed [**conda**] or [**miniconda**], please install firstly,
then install the **alleleCount** in terminal with:

```bash
conda create -n cancerit -c bioconda cancerit-allelecount
```

Install **ASCAT** v3.0 in R console from GitHub with:

```r
# this is a fork version ASCAT
remotes::install_github("ShixiangWang/ascat@v3.0", subdir = "ASCAT")
```

Install **gcap** in R console from GitHub with:

```r
remotes::install_github("ShixiangWang/gcap")
```

## Example

### Pipeline

- for one tumor-normal pair, you can refer to [one-pair.R](test-workflow/one-pair.R).
- for multiple tumor-normal pairs, you can refer to [two-pair.R](test-workflow/two-pair.R).

### Functions

For more custom and advanced control of the analysis, you can read the well
organized function list at [*package site*](https://shixiangwang.github.io/gcap/reference/index.html).


## LICENSE

Apache License (>= 2)
