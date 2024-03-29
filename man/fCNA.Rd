% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fCNA-object.R
\name{fCNA}
\alias{fCNA}
\title{R6 class representing focal copy number amplification list predicted from a cohort}
\description{
Contains fields storing data and methods to get, process and visualize
fCNA information. Examples please see \code{\link[=gcap.ASCNworkflow]{gcap.ASCNworkflow()}}.
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{data}}{a \code{data.table} storing fCNA list, which typically contains following columns:
\itemize{
\item \code{sample} sample or case ID.
\item \code{band} chromosome cytoband.
\item \code{gene_id} gene ID, typically Ensembl ID. You can convert the ID with R package \code{IDConverter}.
\item \code{total_cn} total copy number value.
\item \code{minor_cn} copy number value for minor allele.
\item \code{prob} the probability the gene located in circular DNA.
\item \code{gene_class} gene level amplicon classification.
}}

\item{\code{sample_summary}}{a \code{data.table} storing sample summary data, which typically contains
at least the following columns:
\itemize{
\item \code{sample} sample or case ID. \strong{Should only include cases have been called with GCAP workflow,
otherwise the extra cases would be automatically classified as 'nofocal' (i.e. \code{NA} in \code{sample_summary} field) class}.
\item \code{purity}, \code{ploidy} for tumor purity or ploidy.
\item \code{AScore} aneuploidy score.
\item \code{pLOH} genome percentage harboring LOH events.
\item \verb{CN1 ... CN19} activity of copy number signatures.
\item \strong{\code{class}} the sample class based on amplicon type.
\item \code{ec_genes} number of genes predicted as located on circular DNA.
\item \code{ec_possibly_genes} same with \code{ec_genes} but with less confidence.
\item \code{ec_cytobands} number of cytobands predicted as located on circular DNA.
(the regions of \code{ec_possibly_genes} are not included in computation)
}}
}
\if{html}{\out{</div>}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{min_prob}}{check \verb{$new()} method for details. If you updated this value,
a function will be called to update the sample summary.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-fCNA-new}{\code{fCNA$new()}}
\item \href{#method-fCNA-subset}{\code{fCNA$subset()}}
\item \href{#method-fCNA-getSampleSummary}{\code{fCNA$getSampleSummary()}}
\item \href{#method-fCNA-getGeneSummary}{\code{fCNA$getGeneSummary()}}
\item \href{#method-fCNA-getCytobandSummary}{\code{fCNA$getCytobandSummary()}}
\item \href{#method-fCNA-saveToFiles}{\code{fCNA$saveToFiles()}}
\item \href{#method-fCNA-convertGeneID}{\code{fCNA$convertGeneID()}}
\item \href{#method-fCNA-print}{\code{fCNA$print()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-new"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-new}{}}}
\subsection{Method \code{new()}}{
Create a \code{fCNA} object.
Typically, you can obtain this object from \code{\link[=gcap.workflow]{gcap.workflow()}} or \code{\link[=gcap.ASCNworkflow]{gcap.ASCNworkflow()}}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$new(
  fcna,
  pdata = fcna[, "sample", drop = FALSE],
  min_prob = 0.6,
  only_oncogenes = FALSE,
  genome_build = c("hg38", "hg19", "mm10")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{fcna}}{a \code{data.frame} storing focal copy number amplicon list.}

\item{\code{pdata}}{a \code{data.frame} storing phenotype or sample-level related data. (Optional)}

\item{\code{min_prob}}{the minimal aggregated (in cytoband level) probability to determine a circular amplicon.}

\item{\code{only_oncogenes}}{only_oncogenes if \code{TRUE}, only known oncogenes are kept for circular prediction.}

\item{\code{genome_build}}{genome version}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-subset"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-subset}{}}}
\subsection{Method \code{subset()}}{
Return a subset \code{fCNA} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$subset(..., on = c("data", "sample_summary"))}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{subset expressions on \code{fCNA$data} or \code{fCNA$sample_summary}.}

\item{\code{on}}{if it is "data", subset operations are on data field of \code{fCNA} object,
same for "sample_summary".}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a \code{fCNA}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-getSampleSummary"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-getSampleSummary}{}}}
\subsection{Method \code{getSampleSummary()}}{
Get sample summary of fCNA
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$getSampleSummary(
  only_oncogenes = FALSE,
  genome_build = c("hg38", "hg19", "mm10")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{only_oncogenes}}{only_oncogenes if \code{TRUE}, only known oncogenes are kept for circular prediction.}

\item{\code{genome_build}}{genome version.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a \code{data.table}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-getGeneSummary"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-getGeneSummary}{}}}
\subsection{Method \code{getGeneSummary()}}{
Get gene level summary of fCNA type
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$getGeneSummary(return_mat = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{return_mat}}{if \code{TRUE}, return a cytoband by sample matrix instead of a summary.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a \code{data.table} or a \code{matrix}.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-getCytobandSummary"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-getCytobandSummary}{}}}
\subsection{Method \code{getCytobandSummary()}}{
Get cytoband level summary of fCNA type
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$getCytobandSummary(unique = FALSE, return_mat = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{unique}}{if \code{TRUE}, count sample frequency instead of gene frequency.}

\item{\code{return_mat}}{if \code{TRUE}, return a cytoband by sample matrix instead of a summary.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a \code{data.table}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-saveToFiles"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-saveToFiles}{}}}
\subsection{Method \code{saveToFiles()}}{
Save the key data to local files
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$saveToFiles(dirpath, fileprefix = "fCNA")}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dirpath}}{directory path storing output files.}

\item{\code{fileprefix}}{file prefix. Two result files shall be generated.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-convertGeneID"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-convertGeneID}{}}}
\subsection{Method \code{convertGeneID()}}{
Convert Gene IDs between Ensembl and Hugo Symbol System
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$convertGeneID(
  type = c("ensembl", "symbol"),
  genome_build = c("hg38", "hg19", "mm10")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{type}}{type of input IDs, could be 'ensembl' or 'symbol'.}

\item{\code{genome_build}}{reference genome build.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-fCNA-print"></a>}}
\if{latex}{\out{\hypertarget{method-fCNA-print}{}}}
\subsection{Method \code{print()}}{
print the fCNA object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{fCNA$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{unused.}
}
\if{html}{\out{</div>}}
}
}
}
