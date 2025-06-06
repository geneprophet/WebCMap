% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connectivity_analysis.R
\name{negative_query_CMap}
\alias{negative_query_CMap}
\title{To identify candidate signatures that are negatively correlated with the query signature by computing six kinds of connectivity scores}
\usage{
negative_query_CMap(signature_index)
}
\arguments{
\item{signature_index}{}
}
\value{
a \code{data.frame}
}
\description{
To identify candidate signatures that are negatively correlated with the query signature by computing six kinds of connectivity scores
}
\examples{
data(query_signature)
K=50
res = negative_query_CMap(519170)
}
