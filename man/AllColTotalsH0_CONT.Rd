% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finitepop.R
\name{AllColTotalsH0_CONT}
\alias{AllColTotalsH0_CONT}
\title{Finds all column totals for Compliers and Never Takers
   under the sharp null for Compliers.}
\usage{
AllColTotalsH0_CONT(n_y0x0z0.H0, n_y1x0z0.H0, n_y0x0z1.H0, n_y1x0z1.H0,
  n_y0x1z1.H0, n_y1x1z1.H0)
}
\arguments{
\item{n_y0x0z0.H0}{Number of individuals with Y=0, X=0, Z=0.}

\item{n_y1x0z0.H0}{Number of individuals with Y=1, X=0, Z=0.}

\item{n_y0x0z1.H0}{Number of individuals with Y=0, X=0, Z=1.}

\item{n_y1x0z1.H0}{Number of individuals with Y=1, X=0, Z=1.}

\item{n_y0x1z1.H0}{Number of individuals with Y=0, X=1, Z=1.}

\item{n_y1x1z1.H0}{Number of individuals with Y=1, X=1, Z=1.}
}
\value{
A data.table with all possible combinations of the
   column totals for Compliers and Never Takers
   under the sharp null for Compliers.
}
\description{
Finds all compatible column totals for Compliers and Never Takers
   under the sharp null for Compliers,
   based on an observed dataset.
}
\examples{
AllColTotalsH0_CONT(158, 14, 52, 12, 23, 78)
}

