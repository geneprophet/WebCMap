% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ENSEMBLIDtoENTREZ}
\alias{ENSEMBLIDtoENTREZ}
\title{Convert ENSEMBL ID to ENTREZ ID}
\usage{
ENSEMBLIDtoENTREZ(x)
}
\arguments{
\item{a}{vector of ENSEMBL IDs}
}
\value{
a vector of ENTREZ IDs
}
\description{
Convert ENSEMBL ID to ENTREZ ID
}
\examples{
entrezid <- ENSEMBLIDtoENTREZ("ENSG00000130203")
}
