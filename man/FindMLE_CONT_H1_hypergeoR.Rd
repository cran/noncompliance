% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finitepop.R
\name{FindMLE_CONT_H1_hypergeoR}
\alias{FindMLE_CONT_H1_hypergeoR}
\title{Maximum Likelihood Estimate without assuming the sharp null for Compliers.}
\usage{
FindMLE_CONT_H1_hypergeoR(n_y0x0z0, n_y1x0z0, n_y0x0z1, n_y1x0z1, n_y0x1z1,
  n_y1x1z1)
}
\arguments{
\item{n_y0x0z0}{Number of individuals with Y=0, X=0, Z=0.}

\item{n_y1x0z0}{Number of individuals with Y=1, X=0, Z=0.}

\item{n_y0x0z1}{Number of individuals with Y=0, X=0, Z=1.}

\item{n_y1x0z1}{Number of individuals with Y=1, X=0, Z=1.}

\item{n_y0x1z1}{Number of individuals with Y=0, X=1, Z=1.}

\item{n_y1x1z1}{Number of individuals with Y=1, X=1, Z=1.}
}
\value{
The maximum likelihood, and the corresponding (possibly non-unique)
2 by 4 contingency table.
}
\description{
Find the maximum likelihood estimate of the 2 by 4 contingency table
   assuming only Compliers and Never Takers in the population,
   with the multivariate hypergeometric sampling distribution.
}
\examples{
FindMLE_CONT_H1_hypergeoR(158, 14, 52, 12, 23, 78)
}

