% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{p0}
\alias{p0}
\title{Estimate m using p0 method}
\usage{
p0(data, e = 1, w = 1, d = 0, lag = 0, phi = 0, poisson = 0)
}
\arguments{
\item{data}{a vector of integer-valued colony counts.}

\item{e}{plating efficiency; a positive number not bigger than 1.}

\item{w}{mutant relative fitness; a positive number.}

\item{d}{death probability of wild-type and mutant cells; a non-negative number smaller than 0.5; relative (d)eath rate
and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).}

\item{lag}{phenotypic lag; a non-negative number.}

\item{phi}{relative size of the inoculum (N0/Nt); a non-negative number.}

\item{poisson}{average number of residual Poisson-distributed mutations on the plate; a non-negative number.}
}
\value{
a single non-negative value of m.
}
\description{
Estimate m using p0 method
}
\examples{
p0(c(0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 3, 1, 0, 10, 0, 0, 0, 1))
}
