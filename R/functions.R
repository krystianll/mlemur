### functions for calculating ML point and interval estimates of m as well as LRT P values ###

#' Simulate fluctuation test
#'
#' @description This function simulates a fluctuation assay consisting of n cultures, starting with N0 wild-type cells
#' growing exponentially until reaching Nt wild-type cells. If cv is not 0, the number of wild-type cells per culture
#' is drawn from gamma distribution. Wild-type cell growth is taken as deterministic while mutant cell growth is taken
#' as stochastic according to a simple birth-and-death process. If plating efficiency is less than perfect, the success
#' of plating each mutant cell is simulated using binomial distribution. 
#' @param n number of parallel cultures in the experiment; a positive integer.
#' @param rate average mutation rate; a positive number smaller than 1.
#' @param N0 size of inoculum; a positive integer.
#' @param Nt final number of cells in culture; a positive integer.
#' @param mut_fit relative fitness of the mutant cells vs. wild-type cells; a positive number. Should not be used together
#' with lag.
#' @param death_prob death probability (same for mutants and wild-type cells); a non-negative number smaller
#' than 0.5; relative (d)eath rate and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @param lag mean phenotypic lag of mutant cells in generations; a non-negative number. Should not be used together
#' with mut_fit.
#' @param e plating efficiency; a positive number not bigger than 1.
#' @param cv coefficient of variation of the number of cells in culture; a non-negative number usually smaller than 1.
#' @param trim a non-negative integer; if specified, any mutant cell count bigger than this number will be replaced
#' with this number.
#' @param ret return parameter; either "list" (default) or "vector".
#' @return either a list of length two, each containing a vector of length n: "$mc" containing colony counts in parallel
#' cultures while "$nt" contains the final number of cells in each culture; or a vector identical to "$mc".
#' @export
#' @examples
#' rluria()
#' rluria(n=100, rate=1e-7, lag = 2, trim=5000, ret="v")
#' rluria(n=50, rate=1e-9, cv=0.5)
#' @references
#' Ycart B, Veziris N. Unbiased estimation of mutation rates under fluctuating final counts. PLoS One.
#' 2014;9. doi:10.1371/journal.pone.0101434
#' @references
#' Zheng Q. A second look at the final number of cells in a fluctuation experiment. J Theor Biol.
#' 2016;401: 54–63. doi:10.1016/j.jtbi.2016.04.027
#' @references
#' Mazoyer A, Drouilhet R, Despréaux S, Ycart B. Flan: An R package for inference on mutation models. R J.
#' 2017;9: 334–351. doi:10.32614/rj-2017-029
#' @references
#' Zheng Q. rSalvador: An R package for the fluctuation experiment. G3 Genes, Genomes, Genet.
#' 2017;7: 3849–3856. doi:10.1534/g3.117.300120
rluria <- function(n=10, rate=1e-8, N0=1, Nt=1e9, mut_fit=1.0, death_prob=0, lag=0, e=1, cv=0, trim=0, ret="list") {
  if (lag!=0 && mut_fit != 1.0) {warning("Phenotypic lag should not be used together with mutant fitness.")}
  ret <- match.arg(ret, c("list", "vector"))
  if (ret=="list") {
    return(rluria_list(n=n, rate=rate, N0=N0, Nt=Nt, mut_fit=mut_fit, death_prob=death_prob, lag=lag, e=e, cv=cv, trim=trim))
  } else {
    if (cv!=0) warning("Return type is set to \'vector\' but cv is not 0. The information about total culture sizes will be lost.")
    return(rluria_vec(n=n, rate=rate, N0=N0, Nt=Nt, mut_fit=mut_fit, death_prob=death_prob, lag=lag, e=e, cv=cv, trim=trim))
  }
}

# Auxiliary sequence for probability computation
# 
# @description For internal use only.
# @references
# Abramowitz M, Stegun I. Handbook of Mathematical Functions. Washington: United States Department of Commerce; 1964
# @references
# Stewart FM. Fluctuation analysis: the effect of plating efficiency. Genetica. 1991;84: 51–55. doi:10.1007/BF00123984
# @references
# Jones ME. Luria-Delbrück Fluctuation Experiments; Accounting Simultaneously for Plating Efficiency and Differential
# Growth Rate. Journal of Theoretical Biology. 1994. pp. 355–363. doi:10.1006/jtbi.1994.1032
# @references
# Angerer WP. A note on the evaluation of fluctuation experiments. Mutat Res - Fundam Mol Mech Mutagen.
# 2001;479: 207–224. doi:10.1016/S0027-5107(01)00203-2
# @references
# Zheng Q. Statistical and algorithmic methods for fluctuation analysis with SALVADOR as an implementation. Math Biosci.
# 2002;176: 237–252. doi:10.1016/S0025-5564(02)00087-1
# @references
# Zheng Q. New algorithms for Luria-Delbrück fluctuation analysis. Math Biosci.
# 2005;196: 198–214. doi:10.1016/j.mbs.2005.03.011
# @references
# Zheng Q. On Bartlett’s formulation of the Luria-Delbrück mutation model. Math Biosci.
# 2008;215: 48–54. doi:10.1016/j.mbs.2008.05.005
# @references
# Zheng Q. A note on plating efficiency in fluctuation experiments. Math Biosci.
# 2008;216: 150–153. doi:10.1016/j.mbs.2008.09.002
aux.seq <- function(e = 1, w = 1, death = 0, lag = 0, phi = 0, n = 10, integrate = FALSE) {
  
  error <- 0
  
  if (integrate == FALSE) {
    
    if (death==0 & lag==0 & phi==0) {
      
      seq <- tryCatch(aux_seq(e = e, w = w, n = n), error = function(err) {seq <- NULL})
      
      if (is.null(seq)) {error <- 1}
      else if (!all(is.finite(seq))) {error <- 1}
      else if (any(is.logical(seq))) {error <- 1}
      else if (any(seq[-1] <= 0)) {error <- 1}
      else if (length(seq) < 2) {error <- 1}
      
      if (error == 1) {seq <- tryCatch(aux_seq(e = e, w = w, n = n, option = 2), error = function(err) {seq <- NULL})}
      
    } else if (lag!=0 & w==1 & death==0 & phi==0) {
      
      seq <- tryCatch(aux_seq_lag_s_ext(lag = lag, e = e, n = n), error = function(err) {seq <- NULL})
      error <- 0
      
      if (is.null(seq)) {error <- 1}
      else if (!all(is.finite(seq))) {error <- 1}
      else if (any(is.logical(seq))) {error <- 1}
      else if (any(seq[-1] <= 0)) {error <- 1}
      else if (length(seq) < 2) {error <- 1}
      
      if (error == 1) {seq <- tryCatch(aux_seq_lag_s_ext(lag = lag, e = e, n = n, boost = TRUE), error = function(err) {seq <- NULL})}
      
    } else if (death!=0 & lag==0 & phi==0) {
      
      seq <- tryCatch(aux_seq_death_ext(e = e, w = w, d = death/(1-death), n = n), error = function(err) {seq <- NULL})
      
      if (is.null(seq)) {error <- 1}
      else if (!all(is.finite(seq))) {error <- 1}
      else if (any(is.logical(seq))) {error <- 1}
      else if (any(seq[-1] <= 0)) {error <- 1}
      else if (length(seq) < 2) {error <- 1}
      
      if (error == 1) {seq <- tryCatch(aux_seq_death_ext(e = e, w = w, d = death/(1-death), n = n, boost = TRUE), error = function(err) {seq <- NULL})}
      
    } else {
      
      return(aux.seq(e = e, w = w, death = death, lag = lag, phi = phi, n = n, integrate = TRUE))
      
    }
    
    error <- 0
    
    if (is.null(seq)) {error <- 1}
    else if (!all(is.finite(seq))) {error <- 1}
    else if (any(is.logical(seq))) {error <- 1}
    else if (any(seq[-1] <= 0)) {error <- 1}
    else if (length(seq) < 2) {error <- 1}
    
  }
  
  if (error == 1 | integrate == TRUE) {
    
    error <- 0
    
    seq <- aux_seq_integrate_s(e=e, w=w, d=death/(1-death), lag=lag, phi=phi, n=n)
    
    if (is.null(seq)) {error <- 1}
    else if (!all(is.finite(seq))) {error <- 1}
    else if (any(is.logical(seq))) {error <- 1}
    else if (any(seq[-1] <= 0)) {error <- 1}
    else if (length(seq) < 2) {error <- 1}
    
  }
  
  if (error == 1) {
    
    return(NA)
    
  } else {
    
    return(as.numeric(seq)/(1 - death/(1-death)))
    
  }
}

#' Estimate m using p0 method
#'
#' @param data a vector of integer-valued colony counts.
#' @param e plating efficiency; a positive number not bigger than 1.
#' @param w mutant relative fitness; a positive number.
#' @param d death probability of wild-type and mutant cells; a non-negative number smaller than 0.5; relative (d)eath rate
#' and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @param lag phenotypic lag; a non-negative number.
#' @param phi relative size of the inoculum (N0/Nt); a non-negative number.
#' @param poisson average number of residual Poisson-distributed mutations on the plate; a non-negative number.
#' @return a single non-negative value of m.
#' @export
#' @examples
#' p0(c(0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 3, 1, 0, 10, 0, 0, 0, 1))
p0 <- function(data, e=1, w=1, d=0, lag=0, phi=0, poisson=0) {
  if (w!=1 && lag!=0) {warning("Phenotypic lag should not be used together with mutant fitness.")}
  if (all(data == 0)) return(NA)
  m <- log(sum(data == 0)/length(data))/aux.seq(e = e, w = w, death = d, lag = lag, phi = phi, n = 1)[1] - poisson
  return(m)
}

#' Estimate m using Lea-Coulson method
#' 
#' @description The basic formula comes from a well-known paper by Lea & Coulson. Corrections for partial plating,
#' phenotypic delay (*only non-stochastic*), residual mutations, and differential growth were taken from other Angerer and Koch.
#' Correction for mutant death is based on the Angerer's rationale regarding phenotypic lag, with Newcombe correction
#' for extra cellular divisions.
#' @param data a vector of integer-valued colony counts.
#' @param e plating efficiency; a positive number not bigger than 1.
#' @param w mutant relative fitness; a positive number.
#' @param poisson average number of residual Poisson-distributed mutations on the plate; a non-negative number.
#' @param lag phenotypic lag (*only non-stochastic*); a non-negative number.
#' @param death death probability of wild-type and mutant cells; a non-negative number smaller than 0.5; relative (d)eath rate
#' and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @param confint if TRUE, 95 percent confidence intervals will be estimated by bootstrap using boot package.
#' @return a single non-negative value of m, or a vector of length 3 containing estimate of m and lower and upper limits
#' of 95 percent CI.
#' @export
#' @examples
#' lea.coulson.median(c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' @references
#' Lea DE, Coulson CA. The distribution of the numbers of mutants in bacterial populations. J Genet.
#' 1949;49: 264–285. doi:10.1007/BF02986080
#' @references
#' Newcombe HB. Delayed Phenotypic Expression of Spontaneous Mutations in Escherichia Coli. Genetics.
#' 1948;33: 447–476. doi:10.1093/genetics/33.5.447
#' @references
#' Koch AL. Mutation and growth rates from Luria-Delbrück fluctuation tests. Mutat Res - Fundam Mol Mech Mutagen.
#' 1982;95: 129–143. doi:10.1016/0027-5107(82)90252-4
#' @references
#' Angerer WP. A note on the evaluation of fluctuation experiments. Mutat Res - Fundam Mol Mech Mutagen.
#' 2001;479: 207–224. doi:10.1016/S0027-5107(01)00203-2
lea.coulson.median <- function(data, e=1, w=1, poisson=0, lag=0, death=0, confint=FALSE) {
  if (w!=1 && lag!=0) {warning("Phenotypic lag should not be used together with mutant fitness.")}
  L <- 2^lag
  d <- (1-death)/(1-2*death)
  
  mut <- median(data)/e/w-poisson
  if (mut<=0) {
    return(NA)
  } else {
    m <- (1 + sqrt(1+4*mut*d*L))/2
    for (j in 1:500) {
      f <- m*d/(mut)-1/(log(m/L)+1.24)
      df.dm <- 1*d/(mut)+(1/m)*(1/(log(m/L)+1.24)^2)
      ratio <- f/df.dm
      if (any(!is.finite(c(m, f, df.dm, ratio)))) return(NA)
      if (abs(ratio)/m<1e-6) {
        break
      } else {
        m <- m-ratio
      }
    }
    if (confint == TRUE) {
      return(c(m, c(boot::boot.ci(boot::boot(data, statistic = function(x, y) {lea.coulson.median(x[y], e = e, w = w, poisson = poisson, lag = lag, death = death, confint = FALSE)}, R=10000), type = "perc")[[4]][c(4,5)])))
    } else {
      return(m)
    }
    
  }
}

#' Estimate m using Drake method
#' 
#' @description The basic formula comes from the Drake paper. Corrections for partial plating,
#' phenotypic delay (*only non-stochastic*), residual mutations, and differential growth were taken from other Angerer and Koch.
#' Correction for mutant death is based on the Angerer's rationale regarding phenotypic lag, with Newcombe correction
#' for extra cellular divisions.
#' @param data a vector of integer-valued colony counts.
#' @param e plating efficiency; a positive number not bigger than 1.
#' @param w mutant relative fitness; a positive number.
#' @param poisson average number of residual Poisson-distributed mutations on the plate; a non-negative number.
#' @param lag phenotypic lag (*only non-stochastic*); a non-negative number.
#' @param death death probability of wild-type and mutant cells; a non-negative number smaller than 0.5; relative (d)eath rate
#' and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @return a single non-negative value of m.
#' @export
#' @examples
#' drake.formula(c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' @references
#' Newcombe HB. Delayed Phenotypic Expression of Spontaneous Mutations in Escherichia Coli. Genetics.
#' 1948;33: 447–476. doi:10.1093/genetics/33.5.447
#' @references
#' Koch AL. Mutation and growth rates from Luria-Delbrück fluctuation tests. Mutat Res - Fundam Mol Mech Mutagen.
#' 1982;95: 129–143. doi:10.1016/0027-5107(82)90252-4
#' @references
#' Drake JW. A constant rate of spontaneous mutation in DNA-based microbes. Proc Natl Acad Sci U S A.
#' 1991;88: 7160–7164. doi:10.1073/pnas.88.16.7160
#' @references
#' Angerer WP. A note on the evaluation of fluctuation experiments. Mutat Res - Fundam Mol Mech Mutagen.
#' 2001;479: 207–224. doi:10.1016/S0027-5107(01)00203-2
drake.formula <- function(data, e=1, w=1, poisson=0, lag=0, death=0) {
  if (w!=1 && lag!=0) {warning("Phenotypic lag should not be used together with mutant fitness.")}
  d <- (1-death)/(1-2*death)
  mut <- median(data)/e/w-poisson
  if (mut<=0) {
    return(NA)
  } else {
    L <- 2^lag
    m <- (1 + sqrt(1+4*mut*d*L))/2
    for (j in 1:500) {
      f <- m*d/(mut)-1/(log(m/L))
      df.dm <- 1*d/(mut)+(1/m)*(1/(log(m/L))^2)
      ratio <- f/df.dm
      if (abs(ratio)/m<1e-6) {
        break
      } else {
        m <- m-ratio
      }
    }
    return(m)
  }
}

#' Estimate m using Jones method of the median
#'
#' @description The algorithm comes from rSalvador (jones.median.plating). It is included for the sake of
#' self-containedness.
#' @param data a vector of integer-valued colony counts.
#' @param eff plating efficiency; a positive number not bigger than 1.
#' @return a single non-negative value of m.
#' @export
#' @examples
#' jones.median(c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' @references
#' Zheng Q. rSalvador: An R package for the fluctuation experiment. G3 Genes, Genomes, Genet.
#' 2017;7: 3849–3856. doi:10.1534/g3.117.300120
jones.median <- function(data, eff=1) {
  m <- rep(0, length(data))
  for (i in 1: length(data)) {
    m[i] <- (data[i]/eff-0.693)/((log(data[i])/eff)+0.367)
  }
  m <- median(m)
  return(m)
}

#' Estimate m using Koch method of quartiles
#'
#' @description The formula comes from Rosche and Foster (2000).
#' @param data a vector of integer-valued colony counts.
#' @param lag phenotypic lag; a non-negative number.
#' @return a vector of length 4 containing mean m as well as first, second and third quartile.
#' @export
#' @examples
#' koch.quartiles(c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' @references
#' Rosche WA, Foster PL. Determining Mutation Rates in Bacterial Populations. Methods.
#' 2000;20: 4–17. doi:10.1006/meth.1999.0901
koch.quartiles <- function(data, lag = 0){
  q <- as.vector(quantile(data/2^lag, c(0.25, 0.50, 0.75)))
  m1 <- (1.7335 + 0.4474*q[1] - 0.002755*q[1]*q[1])*2^lag
  m2 <- (1.1580 + 0.2730*q[2] - 0.000761*q[2]*q[2])*2^lag
  m3 <- (0.6658 + 0.1497*q[3] - 0.0001387*q[3]*q[3])*2^lag
  return(c(mean(c(m1, m2, m3)), m1, m2, m3))
}

#' Estimate m using Luria-Delbruck method of the mean
#'
#' @description The formula comes from Rosche and Foster (2000).
#' @param data a vector of integer-valued colony counts.
#' @return a single non-negative value of m.
#' @export
#' @examples
#' luria.delbruck.mean(c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' @references
#' Rosche WA, Foster PL. Determining Mutation Rates in Bacterial Populations. Methods.
#' 2000;20: 4–17. doi:10.1006/meth.1999.0901
luria.delbruck.mean <- function(data) {
  C <- length(data)
  r <- mean(data)
  if (r <= 0) {
    return(NA)
  } else {
    m <- lea.coulson.median(data)
    for (j in 1:500) {
      f <- m*log(m*C)-r
      df.dm <- log(m*C)+1
      ratio <- f/df.dm
      if (any(!is.finite(c(m, f, df.dm, ratio)))) return(NA)
      if (abs(ratio)/m<1e-6) {
        break
      } else {
        m <- m-ratio
      }
    }
  }
  return(m)
}

m_guess <- function(data, e=1, w=1, lag=0, death=0, poisson=0) {
  m <- vector(length = 0, mode = "numeric")
  if (poisson!=0) {
    data2 <- as.integer(unlist(lapply(data - poisson, max, 0)))
    if (median(data2)==0 & any(data2>0)) {
      m <- append(m, p0(data = data2, e = e, w = w, lag = lag, d = death))
      m <- append(m, m[1]/sqrt(1+poisson))
    }
  }
  if (median(data) == 0) {
    m <- append(m, p0(data = data, e = e, w = w, lag = lag, d = death))
  } else {
    m <- append(m, suppressWarnings(lea.coulson.median(data = data, e = e, w = w, poisson = poisson, lag = lag, death = death)))
  }
  return(as.vector(na.exclude(m)))
}

#' Estimate m using Maximum Likelihood Method and calculate 95 percent CI using inverted Likelihood Ratio Test
#' 
#' @description This function finds Maximum Likelihood Estimate of m, the average number of mutations in culture.
#' Then it proceeds to find Likelihood Ratio-based confidence limits for m. It is loosely based on newton.LD,
#' newton.LD.plating, newton.B0, confint.LD, confint.LD.plating, and confint.B0 functions from rSalvador,
#' simplified and optimised to avoid redundant computations. It uses a hybrid Newton-bisection
#' algorithm as well as arbitrary-precision arithmetic if necessary for better stability.
#'
#' @param data a vector of integer-valued colony counts.
#' @param e plating efficiency; a positive number not bigger than 1.
#' @param w mutant relative fitness; a positive number.
#' @param lag phenotypic lag; a non-negative number.
#' @param poisson average number of residual Poisson-distributed mutations on the plate; a non-negative number.
#' @param death death probability of wild-type and mutant cells; a non-negative number smaller than 0.5; relative (d)eath rate
#' and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @param phi relative size of the inoculum (N0/Nt); a non-negative number.
#' @param cv coefficient of variation of the final number of cells in each culture.
#' @param confint if TRUE (default), confidence intervals at ci.level will be calculated.
#' @param ci.level confidence interval size, default 0.95.
#' @param verbose if TRUE, mlemur will print messages to the console.
#' @return a single non-negative value of m, or a vector of length 3 containing MLE of m as well lower and upper limit
#' of 95 percent CI.
#' @export
#' @examples
#' mle.mutations(data = c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' mle.mutations(data = c(71, 19, 4, 32, 12, 74, 0, 35, 8, 13, 9, 5000, 31, 6, 10, 106, 0, 22, 4, 69,
#' 30, 47, 237, 15, 74, 89, 135, 11, 30, 1), lag = 2, cv = 0.3, ci.level = 0.68)
#' @references
#' Zheng Q. Statistical and algorithmic methods for fluctuation analysis with SALVADOR as an implementation. Math Biosci.
#' 2002;176: 237–252. doi:10.1016/S0025-5564(02)00087-1
#' @references
#' Zheng Q. New algorithms for Luria-Delbrück fluctuation analysis. Math Biosci.
#' 2005;196: 198–214. doi:10.1016/j.mbs.2005.03.011
#' @references
#' Zheng Q. On Bartlett’s formulation of the Luria-Delbrück mutation model. Math Biosci.
#' 2008;215: 48–54. doi:10.1016/j.mbs.2008.05.005
#' @references
#' Zheng Q. rSalvador: An R package for the fluctuation experiment. G3 Genes, Genomes, Genet.
#' 2017;7: 3849–3856. doi:10.1534/g3.117.300120
mle.mutations <- function(data, e=1, w=1, lag=0, poisson=0, death=0, phi=0, cv=0, confint=TRUE, ci.level=0.95, verbose=FALSE) {
  if (w!=1 && lag!=0) {warning("Phenotypic lag should not be used together with mutant fitness.")}
  n <- max(data)+1
  if (cv > 0) k <- 1/cv/cv else k <- 0
  seq <- aux.seq(e = e, w = w, death = death, lag = lag, phi = phi, n = n-1)
  if (any(is.na(seq))) {if (verbose) message("m: failed to get the sequence"); if (confint==TRUE) {return(c(NA,NA,NA))} else {return(c(NA))}}
  
  # mutants in the culture
  mg <- as.vector(na.exclude(m_guess(data = data, e = e, w = w, lag = lag, death = death, poisson = poisson)))
  if (length(mg)==0) mg <- 1
  
  m.est <- optim_m(mg[1], mg[1]/10, mg[1]*10, seq, n, data, k, poisson, verbose)
  m1 <- m.est[1]
  
  if (m1 < 0) {if (verbose) message("m: couldn't estimate"); if (confint==TRUE) {return(c(NA,NA,NA))} else {return(c(NA))}}
  else if (verbose) {message(paste("m: found MLE", m1))}
  
  if (confint==FALSE) {
    return(m1)
  } else {
    # confidence intervals
    if (verbose) message(paste("m: constructing confidence intervals at the level of significance \U03B1=", (ci.level)*100, "%", sep = ""))
    chi <- qchisq(ci.level,1)
    if (m1 != 0) {
      U <- m.est[2]
      J <- m.est[3]
      loglik <- m.est[4]
      root <- sqrt(chi/J)
    } else {
      loglik <- logprob(m = 1e-200, len = n, data = data, seq = seq, k = k, poisson = poisson)
      if (!is.finite(loglik)) {loglik <- logprob_boost(xm = 1e-200, len = n, data = data, seq = seq, xk = k, xpoisson = poisson)}
    }
    logprobtail <- loglik-0.5*chi
    
    # lower interval
    if (m1 == 0) {
      m1.low <- 0
    } else {
      m1.low <- root_m(current_m = m1-0.5*root, lower_m = max(m1-2*root, m1*0.1), upper_m = m1, seq = seq, len = n, data = data, k = k, poisson = poisson, lalpha = logprobtail, verbose = verbose)
      if (!is.finite(m1.low) | m1.low < 0) {if (verbose) message("m_low: culdn't estimate"); m1.low <- NA}
    }
    if (verbose) {message(paste("m_low: found MLE", m1.low))}
    
    # upper interval
    m1.up <- root_m(current_m = ifelse(m1 != 0, m1+0.5*root, 0.1), lower_m = ifelse(m1 != 0, m1, 0.01), upper_m = ifelse(m1 != 0, m1+2*root, 1), seq = seq, len = n, data = data, k = k, poisson = poisson, lalpha = logprobtail, verbose = verbose)
    if (!is.finite(m1.up) | m1.up < 0) {if (verbose) message("m_up: couldn't estimate"); m1.up <- NA}
    else if (verbose) {message(paste("m_up: found MLE", m1.up))}
    
    return(c(m1, m1.low, m1.up))
  }
}

#' Calculate p-values using Likelihood Ratio Test
#' 
#' @description This function calculates LRT-based p-values to assess the statistical significance of the
#' differences between two mutation rates (X and Y). It is loosely based on LRT.LD and LRT.LD.plating functions from rSalvador,
#' simplified and optimised to avoid redundant computations. It uses a hybrid Newton-bisection
#' algorithm as well as arbitrary-precision arithmetic if necessary for better stability.
#' @param datax a vector of integer-valued colony counts for strain X.
#' @param datay a vector of integer-valued colony counts.
#' @param ex plating efficiency for strain X; a positive number not bigger than 1.
#' @param ey plating efficiency for strain Y; a positive number not bigger than 1.
#' @param wx mutant relative fitness for strain X; a positive number.
#' @param wy mutant relative fitness for strain Y; a positive number.
#' @param lagx phenotypic lag for strain X; a non-negative number.
#' @param lagy phenotypic lag for strain Y; a non-negative number.
#' @param poissonx average number of residual Poisson-distributed mutations on the plate for strain X; a non-negative number.
#' @param poissony average number of residual Poisson-distributed mutations on the plate for strain Y; a non-negative number.
#' @param deathx death probability of wild-type and mutant cells for strain X; a non-negative number smaller than 0.5;
#' relative (d)eath rate and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @param deathy death probability of wild-type and mutant cells for strain Y; a non-negative number smaller than 0.5;
#' relative (d)eath rate and death (p)robability are connected by the relation d = p/(1-p); p = d/(1+d).
#' @param phix relative size of the inoculum (N0/Nt) for strain X; a non-negative number.
#' @param phiy relative size of the inoculum (N0/Nt) for strain Y; a non-negative number.
#' @param cvx coefficient of variation of the final number of cells in each culture for strain X.
#' @param cvy coefficient of variation of the final number of cells in each culture for strain Y.
#' @param Nx average number of cells in culture for strain X.
#' @param Ny average number of cells in culture for strain Y.
#' @param Mx if known, MLE of m for strain X can be put here to speed up computations; mostly for internal use.
#' @param My if known, MLE of m for strain Y can be put here to speed up computations; mostly for internal use.
#' @param verbose if TRUE, mlemur will print messages to the console.
#' @return a single value of p-value between 0 and 1.
#' @export
#' @examples
#' lrt.mutations(datax = c(26, 9, 16, 34, 15, 25, 77, 13, 14, 19), ex = 0.5,
#' datay = c(67, 12, 112, 24, 2151, 159, 102, 60, 32, 26))
#' @examples
#' lrt.mutations(datax = c(37, 1, 5, 43, 53, 11, 82, 2, 19, 58, 28, 70, 9, 14, 5,
#' 9, 25, 55, 2, 41, 19), datay = c(37, 1, 5, 43, 53, 11, 82, 2, 19, 58, 28,
#' 70, 9, 14, 5, 9, 25, 55, 2, 41, 19), lagy = 2)
#' @references
#' Zheng Q. Statistical and algorithmic methods for fluctuation analysis with SALVADOR as an implementation. Math Biosci.
#' 2002;176: 237–252. doi:10.1016/S0025-5564(02)00087-1
#' @references
#' Zheng Q. New algorithms for Luria-Delbrück fluctuation analysis. Math Biosci.
#' 2005;196: 198–214. doi:10.1016/j.mbs.2005.03.011
#' @references
#' Zheng Q. On Bartlett’s formulation of the Luria-Delbrück mutation model. Math Biosci.
#' 2008;215: 48–54. doi:10.1016/j.mbs.2008.05.005
#' @references
#' Zheng Q. Comparing mutation rates under the Luria–Delbrück protocol. Genetica.
#' 2016;144: 351–359. doi:10.1007/s10709-016-9904-3
#' @references
#' Zheng Q. rSalvador: An R package for the fluctuation experiment. G3 Genes, Genomes, Genet.
#' 2017;7: 3849–3856. doi:10.1534/g3.117.300120
lrt.mutations <- function(datax, datay, ex=1, ey=1, wx=1, wy=1, lagx=0, lagy=0, poissonx=0, poissony=0, deathx=0, deathy=0,
                          phix=0, phiy=0, cvx=0, cvy=0, Nx=1, Ny=1, Mx=NA, My=NA, verbose=FALSE) {
  
  if (!is.finite(Mx)) {Mx <- mle.mutations(data = datax, e = ex, w=wx, lag=lagx, poisson=poissonx, death=deathx, phi=phix, cv = cvx, confint = FALSE, verbose=verbose)} else {Mx <- Mx}
  if (!is.finite(My)) {My <- mle.mutations(data = datay, e = ey, w=wy, lag=lagy, poisson=poissony, death=deathy, phi=phiy, cv = cvy, confint = FALSE, verbose=verbose)} else {My <- My}
  if (!is.finite(Mx) | !is.finite(My)) {if (verbose) message("Couldn't estimate mutation rates"); return(NA)}
  if (verbose) message(paste("Estimated numbers of mutations are", Mx, "and", My))
  if (verbose) message(paste("Estimated mutation rates are", Mx/Nx, "and", My/Ny))
  
  if (Nx>=Ny){
    data1 <- datax; e1 <- ex; w1 <- wx; lag1 <- lagx; poisson1 <- poissonx; death1 <- deathx; phi1 <- phix; cv1 <- cvx; N1 <- Nx; M1 <- Mx
    data2 <- datay; e2 <- ey; w2 <- wy; lag2 <- lagy; poisson2 <- poissony; death2 <- deathy; phi2 <- phiy; cv2 <- cvy; N2 <- Ny; M2 <- My
  } else {
    data1 <- datay; e1 <- ey; w1 <- wy; lag1 <- lagy; poisson1 <- poissony; death1 <- deathy; phi1 <- phiy; cv1 <- cvy; N1 <- Ny; M1 <- My
    data2 <- datax; e2 <- ex; w2 <- wx; lag2 <- lagx; poisson2 <- poissonx; death2 <- deathx; phi2 <- phix; cv2 <- cvx; N2 <- Nx; M2 <- Mx
  }
  
  # mutation rates
  n1 <- max(data1)+1
  n2 <- max(data2)+1
  if (cv1>0) {k1 <- 1/cv1/cv1} else {k1 <- 0}
  if (cv2>0) {k2 <- 1/cv2/cv2} else {k2 <- 0}
  R <- N2/N1
  
  seq1 <- aux.seq(e = e1, w = w1, death = death1, lag = lag1, phi = phi1, n = n1-1)
  seq2 <- aux.seq(e = e2, w = w2, death = death2, lag = lag2, phi = phi2, n = n2-1)
  
  # combined mutation rate
  m.est <- combo_optim_m(current_m = (M1+M2)/2, lower_m = M1, upper_m = M2, R = R, seq1 = seq1, seq2 = seq2, len1 = n1, len2 = n2,
                         data1 = data1, data2 = data2, k1 = k1, k2 = k2, poisson1 = poisson1, poisson2 = poisson2, verbose = verbose)
  if (m.est[1] < 0) {if (verbose) message("m_combined: couldn't estimate m"); return(NA)}
  
  MC <- m.est[1]
  
  if (verbose) message(paste("m_combined: obtained estimate", MC))
  
  # log likelihood functions
  if (MC != 0) {
    l0 <- m.est[2]
  } else {
    la <- logprob(m = 1e-200, len = n1, data = data1, seq = seq1, k = k1, poisson = poisson1)
    if (!is.finite(la)) {la <- logprob_boost(xm = 1e-200, len = n1, data = data1, seq = seq1, xk = k1, xpoisson = poisson1)}
    lb <- logprob(m = 1e-200, len = n2, data = data2, seq = seq2, k = k2, poisson = poisson2)
    if (!is.finite(lb)) {lb <- logprob_boost(xm = 1e-200, len = n2, data = data2, seq = seq2, xk = k2, xpoisson = poisson2)}
    l0 <- la+lb
    if(!(all(is.finite(c(la,lb))))) {message("loglikelihood: non-numeric value encountered"); return(NA)}
  }
  
  lc <- logprob(m = M1, len = n1, data = data1, seq = seq1, k = k1, poisson = poisson1)
  if (!is.finite(lc)) {lc <- logprob_boost(xm = M1, len = n1, data = data1, seq = seq1, xk = k1, xpoisson = poisson1)}
  ld <- logprob(m = M2, len = n2, data = data2, seq = seq2, k = k2, poisson = poisson2)
  if (!is.finite(ld)) {ld <- logprob_boost(xm = M2, len = n2, data = data2, seq = seq2, xk = k2, xpoisson = poisson2)}
  if(!(all(is.finite(c(lc,ld))))) {message("loglikelihood: non-numeric value encountered"); return(NA)}
  l1 <- lc+ld
  
  chi <- -2*(l0-l1)
  pval <- pchisq(chi,1,lower.tail = F)
  return(pval)
}

#' Calculate profile likelihood confidence intervals for an arbitrary function of mutation rates
#'
#' @description Inspired by Zheng 2021, this function calculates MLE and profile likelihood confidence intervals
#' of an arbitrary function of data, such as fold (X1/X2), subtraction (X1-X2), fold with background
#' subtraction ((X1-X2)/(X3-X2)), or a user-defined function utilising basic mathematical operations: addition,
#' subtraction, multiplication or division. Up to 6 datasets can be used.
#' @param data a list of length 2 to 6 containing numeric vectors with colony counts.
#' @param e a list or vector of length 2 to 6 (same as data) containing respective values of plating efficiency.
#' @param w a list or vector of length 2 to 6 (same as data) containing respective values of relative mutant fitness.
#' @param cv a list or vector of length 2 to 6 (same as data) containing respective values of coefficient of variation.
#' @param Nt a list or vector of length 2 to 6 (same as data) containing respective values of average number of cells in culture.
#' @param fun describes the function of mutation rates for which confidence intervals should be calculated. Must be a character
#' vector of length 1 containing either one of default options ("fold X1/X2", "subtraction X1-X2",
#' "double fold (X1/X2)/(X3/X4)", "background subtraction fold (X1-X2)/(X3-X2)") or a user-defined function containing
#' phrases X1, X2, X3, X4, X5, X6, which denote mutation rates for Strain 1, Strain 2, ... . Only addition +, subtraction -,
#' multiplication *, and division /, are currently supported.
#' @return a vector of length 3 containg MLE of an arbitrary function of mutation rates as well as lower and upper
#' profile likelihood confidence limits.
#' @export
#' @examples
#' mle.fold(list(c(26, 9, 16, 34, 15, 25, 77, 13, 14, 19), c(67, 12, 112, 24,
#' 2151, 159, 102, 60, 32, 26)), e = c(0.5, 1))
#' @examples
#' mle.fold(list(c(33,17,15), c(1,4,10), c(45,86,156)), fun="X1*X2/X3")
#' @references
#' Venzon DJ, Moolgavkar SH. A Method for Computing Profile-Likelihood-Based Confidence Intervals. Appl Stat.
#' 1988;37: 87. doi:10.2307/2347496
#' @references
#' Zheng Q. New approaches to mutation rate fold change in Luria–Delbrück fluctuation experiments. Math Biosci.
#' 2021;335: 108572. doi:10.1016/j.mbs.2021.108572
mle.fold <- function(data, e=NULL, w=NULL, cv=NULL, Nt=NULL, fun="fold X1/X2") {
  
  if (!is.list(data) | !length(data) > 1) stop("Colony counts must be a list of length 2 to 6.")
  
  # checking if fold function is correct
  fun <- tryCatch(match.arg(fun, c("fold X1/X2", "subtraction X1-X2", "double fold (X1/X2)/(X3/X4)", "background subtraction fold (X1-X2)/(X3-X2)")), error = function(err) fun)
  if (fun == "fold X1/X2") {function.Y <- "X1/X2"}
  else if (fun == "subtraction X1-X2") {function.Y <- "X1-X2"}
  else if (fun == "double fold (X1/X2)/(X3/X4)") {function.Y <- "(X1/X2)/(X3/X4)"}
  else if (fun == "background subtraction fold (X1-X2)/(X3-X2)") {function.Y <- "(X1-X3)/(X2-X3)"}
  else {function.Y <- fun}
  allowed.characters <- gsub("(X1|X2|X3|X4|X5|X6)(?=$|\\+|\\-|\\*|\\/|\\(|\\)| )", "", function.Y, perl = TRUE)
  allowed.characters <- gsub("\\(|\\)|\\+|\\-|\\/|\\*| ", "", allowed.characters)
  if (allowed.characters != "") stop("There are invalid characters in your fun argument.")
  if (Ryacas::yac_str(paste("Simplify(", function.Y, ")", sep = "")) %in% c("{}", "Undefined", "0")) stop("Invalid fun argument.")
  used.variables <- rep(0, 6)
  for (i in 1:6) {
    used.variables[i] <- ifelse(length(grep(paste("X", i, sep = ""), function.Y)) > 0, 1, 0)
  }
  param.num <- sum(used.variables)
  checked.variables <- 0
  for (i in 1:6) {
    if (checked.variables >= param.num) {break}
    if (used.variables[i] == 1) {checked.variables <- checked.variables + 1}
    else {
      variable.to.drop <- which(used.variables == 1)[checked.variables + 1]
      function.Y <- gsub(paste("X", variable.to.drop, sep = ""), paste("X", i, sep = ""), function.Y)
      used.variables[i] <- 1
      used.variables[variable.to.drop] <- 0
      checked.variables <- checked.variables + 1
    }
  }
  
  # checking if the lengths of provided datasets match fun
  if (length(data) < param.num) {stop(paste("List of colony counts (", length(data), ") is shorter than the number of variables in fun (", sum(used.variables), ").", sep = ""))}
  if (length(data) > param.num) {warning(paste("List of colony counts (", length(data), ") is longer than the number of variables in fun (", sum(used.variables), "). Only first ", sum(used.variables), " positions will be used.", sep = "")); data <- data[1:sum(used.variables)]}
  for (i in 1:param.num) {
    if (length(data[[i]]) < 2) stop(paste("Less than two colony counts in dataset ", i, ".", sep = ""))
  }
  if (is.null(e)) {e <- rep(1, param.num)}
  e <- unlist(e)
  if (length(e) < param.num) {stop("The vector of plating efficiencies is shorter than colony counts list.")}
  else if (length(e) > param.num) {e <- e[1:param.num]}
  if (is.null(cv)) {cv <- rep(0, param.num)}
  cv <- unlist(cv)
  if (length(cv) < param.num) {stop("The vector of coefficients of variation is shorter than colony counts list.")}
  else if (length(cv) > param.num) {cv <- cv[1:param.num]}
  if (is.null(w)) {w <- rep(1, param.num)}
  w <- unlist(w)
  if (length(w) < param.num) {stop("The vector of mutant relative growth rates is shorter than colony counts list.")}
  else if (length(w) > param.num) {w <- w[1:param.num]}
  if (is.null(Nt)) {Nt <- rep(1, param.num)}
  Nt <- unlist(Nt)
  if (length(Nt) < param.num) {stop("The vector of average numbers of cells in cultures is shorter than colony counts list.")}
  else if (length(Nt) > param.num) {Nt <- Nt[1:param.num]}
  
  param.names <- c("Y", "M2", "M3", "M4", "M5", "M6")
  culture.names <- c("N1", "N2", "N3", "N4", "N5", "N6")
  
  # substituting parameters in fold function Y=f(M1,M2,...,N1,N2,...)
  function.Y <- gsub("X1", "(M1/N1)", gsub("X2", "(M2/N2)", gsub("X3", "(M3/N3)", gsub("X4", "(M4/N4)", gsub("X5", "(M5/N5)", gsub("X6", "(M6/N6)", function.Y))))))
  
  # finding mles
  M.hat <- mapply(function(v,x,y,z) {mle.mutations(data = x, e = y, w = z, cv = v, confint = FALSE)}, x = data, y = e, z = w, v = cv)
  Y.hat <- eval(parse(text = gsub("(N)([1-6])", "Nt\\[\\2\\]", gsub("(M)([1-6])", "M.hat\\[\\2\\]", function.Y))))
  
  # expressing M1 as a function of fold M1=f(Y,M2,...,N1,N2,...)
  function.M1 <- gsub(x = gsub(x = Ryacas::yac_str(paste("Solve(", paste("Y==", function.Y, sep = ""), ", M1)")), pattern = "M1==", replacement = ""), pattern = '^.|.$', replacement = '')
  
  # creating function that calculates first derivatives of M1=f(Y,M2,...,N1,N2,...)
  deriv <- mapply(function(A) {Ryacas::yac_str(paste("D(", A, ") (", function.M1, ")", sep=""))}, A = param.names[1:param.num])
  deriv.parsed <- parse(text = gsub("Y", "M[1]", gsub("(N)([1-6])", "N\\[\\2\\]", gsub("(M)([2-6])", "M\\[\\2\\]", deriv))))
  calc.deriv <- function(M,N) {
    res <- vector(length = length(M))
    for (i in 1:length(M)) {
      res[i] <- eval(deriv.parsed[i])
    }
    return(res)
  }
  
  # creating function that calculates second derivatives of M1=f(Y,M2,...,N1,N2,...)
  deriv2 <- matrix(rep(0,param.num^2), nrow = param.num, ncol = param.num, dimnames = list(param.names[1:param.num], param.names[1:param.num]))
  for (i in 1:param.num) {
    deriv2[i,] <- mapply(function(A) {Ryacas::yac_str(paste("D(", A, ") (", deriv[i], ")", sep=""))}, A = param.names[1:param.num])
  }
  deriv2.parsed <- parse(text = gsub("Y", "M[1]", gsub("(N)([1-6])", "N\\[\\2\\]", gsub("(M)([2-6])", "M\\[\\2\\]", deriv2))))
  calc.deriv2 <- function(M,N) {
    res <- matrix(NA, nrow = length(M), ncol = length(M))
    for (i in 1:length(res)) {
      res[i] <- eval(deriv2.parsed[i])
    }
    return(res)
  }
  
  # creating function that calculates M1=f(Y,M2,...,N1,N2,...)
  function.M1.parsed <- parse(text = paste(gsub("Y", "M[1]", gsub("(N)([1-6])", "N\\[\\2\\]", gsub("(M)([2-6])", "M\\[\\2\\]", function.M1)))))
  calc.M1 <-  function(M,N) {eval(function.M1.parsed)}
  
  # calculating pmfs and their derivatives under MLEs
  n <- sapply(data, max) + 1
  
  seq <- mapply(function(e, w, n) aux.seq(e = e, w = w, n = n - 1), e, w, n, SIMPLIFY = FALSE)
  
  probs.hat <- mapply(function(m, cv, seq, n) calc_probs(m, cv, seq, n), m = M.hat, cv = cv, seq = seq, n = n, SIMPLIFY = FALSE)
  
  probs1.ratio <- lapply(1:param.num, function(i) {probs.hat[[i]]$prob1[data[[i]]+1]/probs.hat[[i]]$prob[data[[i]]+1]})
  probs2.ratio <- lapply(1:param.num, function(i) {probs.hat[[i]]$prob2[data[[i]]+1]/probs.hat[[i]]$prob[data[[i]]+1]})
  
  # calculating loglikelihood, score and observed fisher information matrices under MLEs
  deriv.value <- calc.deriv(c(Y.hat, M.hat[-1]), Nt)
  U.hat <- sum(probs1.ratio[[1]])*deriv.value
  for (i in 2:param.num) {
    U.hat[i] <- U.hat[i] + sum(probs1.ratio[[i]])
  }
  
  deriv2.value <- calc.deriv2(c(Y.hat, M.hat[-1]), Nt)
  J.hat <- matrix(rep(0, param.num^2), param.num, param.num)
  for (i in 1:param.num) {
    for (j in 1:param.num) {
      J.hat[i,j] <- J.hat[i,j] + sum(probs1.ratio[[1]]^2) * deriv.value[i] * deriv.value[j] - sum(probs1.ratio[[1]]) * deriv2.value[i,j] - sum(probs2.ratio[[1]]) * deriv.value[i] * deriv.value[j]
    }
  }
  for (i in 2:param.num) {
    J.hat[i,i] <- J.hat[i,i] + sum(probs1.ratio[[i]]^2) - sum(probs2.ratio[[i]])
  }
  
  loglikelihood.hat <- sum(unlist(lapply(1:param.num, function(i) {sum(log(probs.hat[[i]]$prob[data[[i]] + 1]))})))
  
  # calculating delta
  chi <- qchisq(0.95, 1)
  l.alpha <- loglikelihood.hat - 0.5*chi
  
  h <- sqrt(chi / (J.hat[1, 1] - crossprod(J.hat[1, -1], solve(J.hat[-1, -1]) %*% J.hat[-1, 1])))
  delta <- as.vector(h) * rbind(1, -crossprod(solve(-J.hat[-1, -1]), -J.hat[-1, 1]))
  
  # lower limit
  Y.low <- NA
  lower.1 <- rep(-1, param.num)
  grid <- as.matrix(expand.grid(rep(list(seq(0.1, 1.5, length.out=31)), param.num)))
  for (i in 1:nrow(grid)) {
    if (all(lower.1[-1] > 0)) break
    
    lower.0 <- c(Y.hat, M.hat[-1]) - as.vector(grid[i,]) * as.vector(delta)
    
    iter <- 0
    repeat {
      iter <- iter + 1
      if (iter > 30) {lower.1 <- rep(-1, param.num); break}
      probs.low <- mapply(function(m, cv, seq, n) calc_probs(m, cv, seq, n), m = c(calc.M1(lower.0,Nt),lower.0[-1]), cv = cv, seq = seq, n = n, SIMPLIFY = FALSE)
      if (!all(is.finite(unlist(unlist(probs.low))))) {lower.1 <- rep(-1, param.num); break}
      probs1.ratio.low <- lapply(1:param.num, function(i) {probs.low[[i]]$prob1[data[[i]]+1]/probs.low[[i]]$prob[data[[i]]+1]})
      probs2.ratio.low <- lapply(1:param.num, function(i) {probs.low[[i]]$prob2[data[[i]]+1]/probs.low[[i]]$prob[data[[i]]+1]})
      
      deriv.value.low <- calc.deriv(lower.0, Nt)
      U.low <- rep(0, param.num)
      U.low <- sum(probs1.ratio.low[[1]])*deriv.value.low
      for (i in 2:param.num) {
        U.low[i] <- U.low[i] + sum(probs1.ratio.low[[i]])
      }
      
      deriv2.value.low <- calc.deriv2(lower.0, Nt)
      J.low <- matrix(rep(0, param.num^2), param.num, param.num)
      for (i in 1:param.num) {
        for (j in 1:param.num) {
          J.low[i,j] <- J.low[i,j] + sum(probs1.ratio.low[[1]]^2) * deriv.value.low[i] * deriv.value.low[j] - sum(probs1.ratio.low[[1]]) * deriv2.value.low[i,j] - sum(probs2.ratio.low[[1]]) * deriv.value.low[i] * deriv.value.low[j]
        }
      }
      for (i in 2:param.num) {
        J.low[i,i] <- J.low[i,i] + sum(probs1.ratio.low[[i]]^2) - sum(probs2.ratio.low[[i]])
      }
      
      if (!all(is.finite(U.low))) {lower.1 <- rep(-1, param.num); break}
      if (!all(is.finite(J.low))) {lower.1 <- rep(-1, param.num); break}
      
      loglikelihood.low <- sum(unlist(lapply(1:param.num, function(i) {sum(suppressWarnings(log(probs.low[[i]]$prob[data[[i]] + 1])))})))
      if (!all(is.finite(loglikelihood.low))) {lower.1 <- rep(-1, param.num); break}
      
      J.low[1,] <- -U.low
      U.low[1] <- loglikelihood.low - l.alpha
      U.low <- matrix(U.low, ncol = 1)
      
      delta.low <- as.vector(solve(J.low) %*% U.low)
      lower.1 <- lower.0 + delta.low
      
      if (!all(is.finite(lower.1))) {lower.1 <- rep(-1, param.num); break}
      else if (any(lower.1[-1] < 0)) {lower.1 <- rep(-1, param.num); break}
      else if (lower.1[1]>=Y.hat) {lower.1 <- rep(-1, param.num); break}
      if (sqrt(sum((delta.low)^2 / sum(lower.0^2))) < 1e-6) {Y.low <- lower.1[1]; break}
      else {lower.0 <- lower.1}
      
    }
  }
  
  # upper limit
  Y.up <- NA
  upper.1 <- rep(-1, param.num)
  grid <- as.matrix(expand.grid(rep(list(seq(0.2, 4, length.out=21)), param.num)))
  for (i in 1:nrow(grid)) {
    if (all(upper.1[-1] > 0)) break
    
    upper.0 <- c(Y.hat, M.hat[-1]) + as.vector(grid[i,]) * as.vector(delta)
    
    iter <- 0
    repeat {
      iter <- iter + 1
      if (iter > 30) {upper.1 <- rep(-1, param.num); break}
      probs.up <- mapply(function(m, cv, seq, n) calc_probs(m, cv, seq, n), m = c(calc.M1(upper.0,Nt),upper.0[-1]), cv = cv, seq = seq, n = n, SIMPLIFY = FALSE)
      if (!all(is.finite(unlist(unlist(probs.up))))) {upper.1 <- rep(-1, param.num); break}
      probs1.ratio.up <- lapply(1:param.num, function(i) {probs.up[[i]]$prob1[data[[i]]+1]/probs.up[[i]]$prob[data[[i]]+1]})
      probs2.ratio.up <- lapply(1:param.num, function(i) {probs.up[[i]]$prob2[data[[i]]+1]/probs.up[[i]]$prob[data[[i]]+1]})
      
      deriv.value.up <- calc.deriv(upper.0, Nt)
      U.up <- sum(probs1.ratio.up[[1]])*deriv.value.up
      for (i in 2:param.num) {
        U.up[i] <- U.up[i] + sum(probs1.ratio.up[[i]])
      }
      
      deriv2.value.up <- calc.deriv2(upper.0, Nt)
      J.up <- matrix(rep(0, param.num^2), param.num, param.num)
      for (i in 1:param.num) {
        for (j in 1:param.num) {
          J.up[i,j] <- J.up[i,j] + sum(probs1.ratio.up[[1]]^2) * deriv.value.up[i] * deriv.value.up[j] - sum(probs1.ratio.up[[1]]) * deriv2.value.up[i,j] - sum(probs2.ratio.up[[1]]) * deriv.value.up[i] * deriv.value.up[j]
        }
      }
      for (i in 2:param.num) {
        J.up[i,i] <- J.up[i,i] + sum(probs1.ratio.up[[i]]^2) - sum(probs2.ratio.up[[i]])
      }
      
      if (!all(is.finite(U.up))) {upper.1 <- rep(-1, param.num); break}
      if (!all(is.finite(J.up))) {upper.1 <- rep(-1, param.num); break}
      
      loglikelihood.up <- sum(unlist(lapply(1:param.num, function(i) {sum(suppressWarnings(log(probs.up[[i]]$prob[data[[i]] + 1])))})))
      if (!all(is.finite(loglikelihood.up))) {upper.1 <- rep(-1, param.num); break}
      
      J.up[1,] <- -U.up
      U.up[1] <- loglikelihood.up - l.alpha
      U.up <- matrix(U.up, ncol = 1)
      
      delta.up <- as.vector(solve(J.up) %*% U.up)
      upper.1 <- upper.0 + delta.up
      
      if (!all(is.finite(upper.1))) {upper.1 <- rep(-1, param.num); break}
      else if (any(upper.1[-1] < 0)) {upper.1 <- rep(-1, param.num); break}
      else if (upper.1[1]<=Y.hat) {upper.1 <- rep(-1, param.num); break}
      if (sqrt(sum((delta.up)^2 / sum(upper.0^2))) < 1e-6) {Y.up <- upper.1[1]; break}
      else {upper.0 <- upper.1}
      
    }
  }
  
  return(c(Y.hat, Y.low, Y.up))
  
}

# Mutation Rate Calculator
calc.rate.int <- function(data) {
  if (is.na(data$PlatingEfficiency)) {
    eff <- data$VolumeSelective/data$DilutionSelective/data$VolumeTotal
  } else {
    eff <- data$PlatingEfficiency
  }
  
  if(is.na(data$MeanCells)) {
    Nt <- mean(data$CountsNonselective) * data$DilutionNonselective / data$VolumeNonselective * data$VolumeTotal
  } else {
    Nt <- data$MeanCells
  }
  
  if (data$model && as.numeric(data$setCV) == 1) {
    cv <- data$CV
  } else if (data$model && as.numeric(data$setCV) == 0) {
    cv <- sd(data$CountsNonselective) * data$DilutionNonselective / data$VolumeNonselective * data$VolumeTotal / Nt
  } else {
    cv <- 0
  }
  
  if (is.na(cv)) {cv0 <- 0}
  else if (cv<1e-5) {cv0 <- 0}
  else if (cv>10.0) {cv0 <- 10.0}
  else {cv0 <- cv}
  
  if (is.na(data$Fitness)) {data$Fitness <- 1}
  if (is.na(data$Lag)) {data$Lag <- 0}
  if (is.na(data$Residual)) {data$Residual <- 0}
  if (is.na(data$Death)) {data$Death <- 0}
  if (is.na(data$Inoculum)) {data$Inoculum <- 0}
  
  est <- tryCatch(mle.mutations(data = data$CountsSelective, e = eff, w = data$Fitness, lag = data$Lag, poisson = data$Residual,
                                death = data$Death/(1 + data$Death), phi = data$Inoculum/Nt, cv = cv0, confint = T, ci.level = 0.95, verbose = F),
                  error = function(err) {c(NA,NA,NA)})
  
  m <- est[1]
  m_low <- est[2]
  m_up <- est[3]
  mu <- m / Nt
  mu_low <- m_low / Nt
  mu_up <- m_up / Nt
  
  outputData <- c(
    signif(eff, 6),
    formatC(Nt, format = "e", digits = 3),
    ifelse(cv0 != 0, signif(cv0, 4), 0),
    ifelse(data$Fitness != 1, signif(data$Fitness, 4), 1),
    ifelse(data$Lag != 0, signif(data$Lag, 4), 0),
    ifelse(data$Death != 0, signif(data$Death, 4), 0),
    ifelse(data$Residual != 0, signif(data$Residual, 4), 0),
    ifelse(data$Inoculum != 0, signif(data$Inoculum/Nt, 4), 0),
    signif(m, 6),
    signif(m_low, 6),
    signif(m_up, 6),
    signif(mu, 6),
    signif(mu_low, 6),
    signif(mu_up, 6),
    m
  )
  return(outputData)
}

# Calculator of pvalue
calc.pval.int <- function(datax, datay, Mx, My) {
  if (is.na(datax$PlatingEfficiency)) {
    ex <- datax$VolumeSelective/datax$DilutionSelective/datax$VolumeTotal
  } else {
    ex <- datax$PlatingEfficiency
  }
  if(is.na(datax$MeanCells)) {
    Ntx <- mean(datax$CountsNonselective) * datax$DilutionNonselective / datax$VolumeNonselective * datax$VolumeTotal
  } else {
    Ntx <- datax$MeanCells
  }
  if (datax$model && datax$setCV == 1) {
    cvx <- datax$CV
  } else if (datax$model && datax$setCV == 0) {
    cvx <- sd(datax$CountsNonselective) * datax$DilutionNonselective / datax$VolumeNonselective * datax$VolumeTotal / Ntx
  } else {
    cvx <- 0
  }
  if (is.na(cvx)) {cv0x <- 0}
  else if (cvx<1e-5) {cv0x <- 0}
  else if (cvx>10.0) {cv0x <- 10.0}
  else {cv0x <- cvx}
  if (is.na(datax$Fitness)) {datax$Fitness <- 1}
  if (is.na(datax$Lag)) {datax$Lag <- 0}
  if (is.na(datax$Residual)) {datax$Residual <- 0}
  if (is.na(datax$Death)) {datax$Death <- 0}
  if (is.na(datax$Inoculum)) {datax$Inoculum <- 0}
  
  if (is.na(datay$PlatingEfficiency)) {
    ey <- datay$VolumeSelective/datay$DilutionSelective/datay$VolumeTotal
  } else {
    ey <- datay$PlatingEfficiency
  }
  if(is.na(datay$MeanCells)) {
    Nty <- mean(datay$CountsNonselective) * datay$DilutionNonselective / datay$VolumeNonselective * datay$VolumeTotal
  } else {
    Nty <- datay$MeanCells
  }
  if (datay$model && datay$setCV == 1) {
    cvy <- datay$CV
  } else if (datay$model && datay$setCV == 0) {
    cvy <- sd(datay$CountsNonselective) * datay$DilutionNonselective / datay$VolumeNonselective * datay$VolumeTotal / Nty
  } else {
    cvy <- 0
  }
  if (is.na(cvy)) {cv0y <- 0}
  else if (cvy<1e-5) {cv0y <- 0}
  else if (cvy>10.0) {cv0y <- 10.0}
  else {cv0y <- cvy}
  if (is.na(datay$Fitness)) {datay$Fitness <- 1}
  if (is.na(datay$Lag)) {datay$Lag <- 0}
  if (is.na(datay$Residual)) {datay$Residual <- 0}
  if (is.na(datay$Death)) {datay$Death <- 0}
  if (is.na(datay$Inoculum)) {datay$Inoculum <- 0}
  
  result <- tryCatch(lrt.mutations(datax = datax$CountsSelective, datay = datay$CountsSelective, ex = ex, ey = ey, 
                                   wx = datax$Fitness, wy = datay$Fitness,
                                   lagx = datax$Lag, lagy = datay$Lag, poissonx = datax$Residual, poissony = datay$Residual,
                                   deathx = datax$Death/(1 + datax$Death), deathy = datay$Death/(1 + datay$Death), phix = datax$Inoculum/Ntx, phiy = datay$Inoculum/Nty,
                                   cvx = cv0x, cvy = cv0y, Nx = Ntx, Ny = Nty, Mx = Mx, My = My, verbose = F),
                     error=function(err) {c(NA)})
  
  return(result)
}