% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ASCNworkflow.R
\name{gcap.ASCNworkflow}
\alias{gcap.ASCNworkflow}
\title{GCAP workflow for gene-level amplicon prediction from ASCN input}
\usage{
gcap.ASCNworkflow(
  data,
  genome_build = c("hg38", "hg19"),
  model = "XGB11",
  tightness = 1L,
  gap_cn = 3L,
  overlap = 1,
  only_oncogenes = FALSE,
  outdir = getwd(),
  result_file_prefix = paste0("gcap_", uuid::UUIDgenerate(TRUE))
)
}
\arguments{
\item{data}{a \code{data.frame} with following columns. The key columns can be obtained
from common allele specific CNV calling software, e.g., ASCAT, Sequenza, FACETS.
\itemize{
\item chromosome: chromosome names starts with 'chr'.
\item start: start position of the segment.
\item end: end position of the segment.
\item total_cn: total integer copy number of the segment.
\item minor_cn: minor allele integer copy number of the segment. Set it
to \code{NA} if you don't have this data.
\item sample: sample identifier.
\item purity: tumor purity of the sample. Set to \code{1} if you don't know.
\item ploidy (optinal): ploidy value of the sample tumor genome.
\item age (optional): age of the case, use along with \code{gender}.
\item gender (optional): gender of the case, use along with \code{age}.
\item type (optional): cancer type of the case, use along with \code{age} and \code{gender}.
Please refer to \link{gcap.collapse2Genes} to see the supported cancer types.
This info is only used in 'XGB56' model. If you don't use this model, you
don't need to set it.
}}

\item{genome_build}{"hg38" or "hg19".}

\item{model}{model name ("XGB11", "XGB32", "XGB56") or a custom model
from input. 'toy' can be used for test.}

\item{tightness}{a coefficient to times to TCGA somatic CN to set a more strict threshold
as a circular amplicon.
If the value is larger, it is more likely a fCNA assigned to \code{noncircular}
instead of \code{circular}. \strong{When it is \code{NA}, we don't use TCGA somatic CN data as reference}.}

\item{gap_cn}{a gap copy number value.
A gene with copy number above background (\code{ploidy + gap_cn} in general) would be treated as focal amplicon.
Smaller, more amplicons.}

\item{overlap}{the overlap percentage on gene.}

\item{only_oncogenes}{if \code{TRUE}, only known oncogenes are kept for circular prediction.}

\item{outdir}{result output path.}

\item{result_file_prefix}{file name prefix (without directory path) for storing
final model prediction file in CSV format.
Default a unique file name is generated by UUID approach.}
}
\value{
a list of invisible \code{data.table} and corresponding files saved to local machine.
}
\description{
Unlike \link{gcap.workflow}, this function directly uses the allele-specific
copy number data along with some extra sample information to infer
ecDNA genes.
}
\examples{
data("ascn")
data <- ascn
rv <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
data$purity <- 1
rv2 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")
data$age <- 60
data$gender <- "XY"
rv3 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB32")
# If you want to use 'XGB56', you should include 'type' column
data$type <- "LUAD"
rv4 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB56")
# If you only have total integer copy number
data$minor_cn <- NA
rv5 <- gcap.ASCNworkflow(data, outdir = tempdir(), model = "XGB11")

# R6 class fCNA --------------------------------
print(rv)
print(rv$data)
print(rv$sample_summary)
print(rv$gene_summary)
print(rv$cytoband_summary)

# Create a subset fCNA
rv_subset <- rv$subset(total_cn > 10)
nrow(rv$data)
nrow(rv_subset$data)

rv_subset2 <- rv$subset(sample == "TCGA-02-2485-01")
nrow(rv_subset2$data)
unique(rv_subset2$data$sample)

sum_gene <- rv$getGeneSummary()
sum_gene
mat_gene <- rv$getGeneSummary(return_mat = TRUE)
mat_gene

sum_cytoband <- rv$getCytobandSummary()
sum_cytoband
mat_cytoband <- rv$getCytobandSummary(return_mat = TRUE)
mat_cytoband

}
