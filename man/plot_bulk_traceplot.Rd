% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bulk_traceplot.R
\name{plot_bulk_traceplot}
\alias{plot_bulk_traceplot}
\title{Traceplot of the posterior chain of the proportion (i.e., pi) of unspliced (U) counts in each group - bulk RNA-seq data}
\usage{
plot_bulk_traceplot(results, transcript_id)
}
\arguments{
\item{results}{a \code{list} of \code{\linkS4class{data.frame}} objects, 
computed via \code{\link{DifferentialRegulation}} (single-cell RNA-seq data), 
or \code{\link{DifferentialRegulation_bulk}} (bulk RNA-seq data).}

\item{transcript_id}{a character, indicating the transcript to plot.}
}
\value{
A \code{gtable} object.
}
\description{
\code{plot_bulk_traceplot} plots the traceplot of the posterior chain 
of the proportion (i.e., pi) of unspliced (U) counts in each group.
The vertical grey dashed line indicates the burn-in 
(the iterations on the left side of the burn-in are discarded in posterior analyses).
}
\examples{
# see the example of DifferentialRegulation_bulk function:
help(DifferentialRegulation_bulk)

}
\seealso{
\code{\link{DifferentialRegulation_bulk}}
}
\author{
Simone Tiberi \email{simone.tiberi@unibo.it}
}
