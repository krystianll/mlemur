% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{lrt.mutations}
\alias{lrt.mutations}
\title{Calculate p-values using Likelihood Ratio Test}
\usage{
lrt.mutations(
  datax,
  datay,
  ex = 1,
  ey = 1,
  wx = 1,
  wy = 1,
  lagx = 0,
  lagy = 0,
  poissonx = 0,
  poissony = 0,
  deathx = 0,
  deathy = 0,
  phix = 0,
  phiy = 0,
  cvx = 0,
  cvy = 0,
  Nx = 1,
  Ny = 1,
  Mx = NA,
  My = NA,
  verbose = FALSE
)
}
\arguments{
\item{datax}{a vector of integer-valued colony counts for strain X.}

\item{datay}{a vector of integer-valued colony counts for strain Y.}

\item{ex}{plating efficiency for strain X; a positive number not bigger than 1.}

\item{ey}{plating efficiency for strain Y; a positive number not bigger than 1.}

\item{wx}{mutant relative fitness for strain X; a positive number.}

\item{wy}{mutant relative fitness for strain Y; a positive number.}

\item{lagx}{phenotypic lag for strain X; a non-negative number.}

\item{lagy}{phenotypic lag for strain Y; a non-negative number.}

\item{poissonx}{average number of residual Poisson-distributed mutations on the plate for strain X; a non-negative number.}

\item{poissony}{average number of residual Poisson-distributed mutations on the plate for strain Y; a non-negative number.}

\item{deathx}{death probability of wild-type and mutant cells for strain X; a non-negative number smaller than 0.5;
relative (d)eath rate and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).}

\item{deathy}{death probability of wild-type and mutant cells for strain Y; a non-negative number smaller than 0.5;
relative (d)eath rate and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).}

\item{phix}{relative size of the inoculum (N0/Nt) for strain X; a non-negative number.}

\item{phiy}{relative size of the inoculum (N0/Nt) for strain Y; a non-negative number.}

\item{cvx}{coefficient of variation of the final number of cells in each culture for strain X.}

\item{cvy}{coefficient of variation of the final number of cells in each culture for strain Y.}

\item{Nx}{average number of cells in culture for strain X.}

\item{Ny}{average number of cells in culture for strain Y.}

\item{Mx}{if known, MLE of m for strain X can be put here to speed up computations; mostly for internal use.}

\item{My}{if known, MLE of m for strain Y can be put here to speed up computations; mostly for internal use.}

\item{verbose}{if TRUE, mlemur will print messages to the console.}
}
\value{
a single value of p-value between 0 and 1.
}
\description{
This function calculates LRT-based p-values to assess the statistical significance of the
differences between two mutation rates (X and Y). It is loosely based on LRT.LD and LRT.LD.plating functions from rSalvador,
simplified and optimised to avoid redundant computations. It uses a hybrid Newton-bisection
algorithm as well as arbitrary-precision arithmetic if necessary for better stability.
}
\examples{
lrt.mutations(datax = c(26, 9, 16, 34, 15, 25, 77, 13, 14, 19), ex = 0.5,
datay = c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
lrt.mutations(datax = c(37, 1, 5, 43, 53, 11, 82, 2, 19, 58, 28, 70, 9, 14, 5,
9, 25, 55, 2, 41, 19), datay = c(37, 1, 5, 43, 53, 11, 82, 2, 19, 58, 28,
70, 9, 14, 5, 9, 25, 55, 2, 41, 19), lagy = 2)
}
\references{
Zheng Q. Statistical and algorithmic methods for fluctuation analysis with SALVADOR as an implementation. Math Biosci.
2002;176: 237–252. doi:10.1016/S0025-5564(02)00087-1

Zheng Q. New algorithms for Luria-Delbrück fluctuation analysis. Math Biosci.
2005;196: 198–214. doi:10.1016/j.mbs.2005.03.011

Zheng Q. On Bartlett’s formulation of the Luria-Delbrück mutation model. Math Biosci.
2008;215: 48–54. doi:10.1016/j.mbs.2008.05.005

Zheng Q. Comparing mutation rates under the Luria–Delbrück protocol. Genetica.
2016;144: 351–359. doi:10.1007/s10709-016-9904-3

Zheng Q. rSalvador: An R package for the fluctuation experiment. G3 Genes, Genomes, Genet.
2017;7: 3849–3856. doi:10.1534/g3.117.300120
}
