% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{mle.rate}
\alias{mle.rate}
\title{Estimate mutation rate using Maximum Likelihood Method and calculate 95 percent CI using inverted Likelihood Ratio Test
using pairs of colony counts on selective and non-selective medium}
\usage{
mle.rate(
  data,
  e = 1,
  w = 1,
  lag = 0,
  poisson_rate = 0,
  death = 0,
  inoculum = 0,
  Nt,
  confint = TRUE,
  ci.level = 0.95,
  verbose = FALSE,
  show_logprob = FALSE
)
}
\arguments{
\item{data}{a vector of integer-valued colony counts.}

\item{e}{plating efficiency; a positive number not bigger than 1.}

\item{w}{mutant relative fitness; a positive number.}

\item{lag}{phenotypic lag; a non-negative number.}

\item{poisson_rate}{average mutation rate of Poisson-distributed mutations on the plate; a non-negative number smaller than 1.}

\item{death}{death probability of wild-type and mutant cells; a non-negative number smaller than 0.5; relative (d)eath rate
and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).}

\item{inoculum}{number of cells in the inoculum; a non-negative integer.}

\item{Nt}{a vector of integer-valued culture sizes; must be of the same length as data.}

\item{confint}{if TRUE (default), confidence intervals at ci.level will be calculated.}

\item{ci.level}{confidence interval size, default 0.95.}

\item{verbose}{if TRUE, mlemur will print messages to the console.}

\item{show_logprob}{if TRUE and confint is FALSE, the function will return the point estimate of mutation rate and
the value of the log-likelihood function. This is mostly for internal use.}
}
\value{
a single non-negative value of m, or a vector of length 3 containing MLE of m as well lower and upper limit
of 95 percent CI.
}
\description{
This function finds Maximum Likelihood Estimate of mu, the average mutation rate.
Then it proceeds to find Likelihood Ratio-based confidence limits for mu. The algorithm uses pairs of colony
counts on selective and non-selective medium for each test tube. It uses a hybrid Newton-bisection
algorithm as well as arbitrary-precision arithmetic if necessary for better stability.
}
\examples{
mle.rate(data = c(14, 20, 59, 46, 77, 38, 86, 17, 13, 25), Nt = c(1030417500,
635899434, 721044961, 1496257034, 1276590167, 1977307240, 954768054,
1192512084, 632462708, 748020573))
}
\references{
Zheng Q. Statistical and algorithmic methods for fluctuation analysis with SALVADOR as an implementation. Math Biosci.
2002;176: 237–252. doi:10.1016/S0025-5564(02)00087-1

Zheng Q. New algorithms for Luria-Delbrück fluctuation analysis. Math Biosci.
2005;196: 198–214. doi:10.1016/j.mbs.2005.03.011

Zheng Q. On Bartlett’s formulation of the Luria-Delbrück mutation model. Math Biosci.
2008;215: 48–54. doi:10.1016/j.mbs.2008.05.005

Zheng Q. rSalvador: An R package for the fluctuation experiment. G3 Genes, Genomes, Genet.
2017;7: 3849–3856. doi:10.1534/g3.117.300120
}
