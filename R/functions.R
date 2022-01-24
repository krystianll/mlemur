### functions for calculating ML point and interval estimates of m as well as LRT P values ###

rluria <- function(n=10, rate=1e-8, N0=1, Nt=1e9, mut_fit=1.0, type=0, wt_dp=0, mut_dp=0, lag=0, e=1, cv=0, trim=0, ret="list") {
  ret <- match.arg(ret, c("list", "vector"))
  if (ret=="list") {
    return(rluria_list(n=n, rate=rate, N0=N0, Nt=Nt, mut_fit=mut_fit, type=type, wt_dp=wt_dp, mut_dp=mut_dp, lag=lag, e=e, cv=cv, trim=trim))
  } else {
    if (cv!=0) warning("Return type is set to \'vector\' but cv is not 0. The information about total culture sizes will be lost.")
    return(rluria_vec(n=n, rate=rate, N0=N0, Nt=Nt, mut_fit=mut_fit, type=type, wt_dp=wt_dp, mut_dp=mut_dp, lag=lag, e=e, cv=cv, trim=trim))
  }
}

# Auxiliary sequence for probability computation
# modified to make pmf computation independent of e and w for better re-use of algorithms
# error check implemented to avoid numerical instability when e and w are small
# Abramowitz M, Stegun I. Handbook of Mathematical Functions. Washington: United States Department of Commerce; 1964
# Stewart FM. Genetica 1991 doi:10.1007/BF00123984
# Jones ME. Journal of Theoretical Biology 1994 doi:10.1006/jtbi.1994.1032
# Zheng Q. Math Biosci 2002 doi:10.1016/S0025-5564(02)00087-1
# Zheng Q. Math Biosci 2005 doi:10.1016/j.mbs.2005.03.011
# Zheng Q. Math Biosci 2008 doi:10.1016/j.mbs.2008.05.005
# Zheng Q. Math Biosci 2008 doi:10.1016/j.mbs.2008.09.002

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
      
      seq <- tryCatch(aux_seq_lag_ext(L = 2^lag, e = e, n = n), error = function(err) {seq <- NULL})
      error <- 0
      
      if (is.null(seq)) {error <- 1}
      else if (!all(is.finite(seq))) {error <- 1}
      else if (any(is.logical(seq))) {error <- 1}
      else if (any(seq[-1] <= 0)) {error <- 1}
      else if (length(seq) < 2) {error <- 1}
      
      if (error == 1) {seq <- tryCatch(aux_seq_lag_ext(L = 2^lag, e = e, n = n, boost = TRUE), error = function(err) {seq <- NULL})}
      
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
    
    seq <- aux_seq_integrate(e=e, w=w, d=death/(1-death), lag=lag, phi=phi, n=n)
    
    if (is.null(seq)) {error <- 1}
    else if (!all(is.finite(seq))) {error <- 1}
    else if (any(is.logical(seq))) {error <- 1}
    else if (any(seq[-1] <= 0)) {error <- 1}
    else if (length(seq) < 2) {error <- 1}
    
  }
  
  if (error == 1) {
    
    return(NA)
    
  } else {
    
    return(as.numeric(seq))
    
  }
}

# Approximate methods of m estimation
# P0, Jones median, Drake median, and Lea-Coulson median
# Most useful as an initial guess for m in Newton method
# Jones ME. Mutat Res Mutagen Relat Subj 1993 doi:10.1016/0165-1161(93)90146-Q
# Crane GJ, Thomas SM, Jones ME. Mutat Res - Fundam Mol Mech Mutagen 1996 doi:10.1016/0027-5107(96)00009-7
# Rosche WA, Foster PL. Methods 2000 doi:10.1006/meth.1999.0901

p0 <- function(data, e=1, w=1, d=0, lag=0, phi=0, poisson=0) {
  if (all(data == 0)) return(NA)
  m <- log(sum(data == 0)/length(data))/aux.seq(e = e, w = w, death = d, lag = lag, phi = phi, n = 1)[1] - poisson
  return(m)
}

# The basic formula comes from a well-known paper Lea DE, Coulson CA. J Genet. 1949 doi:10.1007/BF02986080
# Corrections for partial plating, phenotypic delay and residual mutation after Angerer WP. Mutat Res - Fundam Mol Mech Mutagen. 2001 doi:10.1016/S0027-5107(01)00203-2
# Correction for differential growth taken from Eq. 12 in Koch AL. Mutat Res - Fundam Mol Mech Mutagen. 1982 doi:10.1016/0027-5107(82)90252-4
# Correction for mutant death based on the Angerer's rationale regarding phenotypic lag
lea.coulson.median <- function(data, e=1, w=1, poisson=0, lag=0, death=0, confint=FALSE) {
  mut <- median(data)/e/w-poisson
  if (mut<=0) {
    return(NA)
  } else {
    L <- 2^lag
    d <- (1-death)/(1-2*death)
    m <- (1 + sqrt(1+4*mut*d))/2
    for (j in 1:500) {
      f <- m/(mut)-1/(log(m/L/d)+1.24)
      df.dm <- 1/(mut)+(1/m)*(1/(log(m/L/d)+1.24)^2)
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

# Modifications were copied from Lea-Coulson
drake.formula <- function(data, eff=1, w=1, poisson=0, lag=0) {
  mut <- median(data)/eff/w-poisson
  if (mut<=0) {
    return(NA)
  } else {
    L <- 2^lag
    m <- (1 + sqrt(1+4*mut))/2
    for (j in 1:500) {
      f <- m/(mut)-1/(log(m/L))
      df.dm <- 1/(mut)+(1/m)*(1/(log(m/L))^2)
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

jones.median <- function(data, eff=1) {
  m <- rep(0, length(data))
  for (i in 1: length(data)) {
    m[i] <- (data[i]/eff-0.693)/((log(data[i])/eff)+0.367)
  }
  mest <- median(m)
  return(mest)
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
    # m <- append(m, jones.median(data, eff = e))
  }
  return(as.vector(na.exclude(m)))
}

### Calculate MLE, CI, and P value under the Bartlett model ###
# Combined functions for B0 (gamma-mixture) distribution indexed by e and cv
# based on functions newton.B0 and confint.B0 from rSalvador 1.8
# P values can be calculated using Likelihood Ratio Test, using an established framework
# simplified and optimised to avoid redundant computations, partial plating is possible in all cases
# arbitrary precision arithmetics has been implemented and trial-and-error iterative initial value of m approach is used when necessary to avoid under- and overflow
# Zheng Q. Math Biosci 2008 doi:10.1016/j.mbs.2008.09.002
# Zheng Q. Statistics (Ber) 2010 doi:10.1080/02331880903236868
# Zheng Q. Genetica 2016 doi:10.1007/s10709-016-9904-3
mle.mutations <- function(data, e=1, w=1, lag=0, poisson=0, death=0, phi=0, cv=0, confint=TRUE, ci.level=0.95, verbose=FALSE) {
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
  if (is.na(data$WildTypeDeath)) {data$WildTypeDeath <- 0}
  if (is.na(data$MutantDeath)) {data$MutantDeath <- 0}
  if (is.na(data$Inoculum)) {data$Inoculum <- 0}
  
  # data$Fitness <- data$Fitness/(1 - data$WildTypeDeath)
  
  est <- tryCatch(mle.mutations(data = data$CountsSelective, e = eff, w = data$Fitness/(1 - data$WildTypeDeath), lag = data$Lag, poisson = data$Residual,
                                death = data$MutantDeath/(1 + data$MutantDeath), phi = data$Inoculum/Nt, cv = cv0, confint = T, ci.level = 0.95, verbose = F),
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
    ifelse(data$WildTypeDeath != 0, signif(data$WildTypeDeath, 4), 0),
    ifelse(data$MutantDeath != 0, signif(data$MutantDeath, 4), 0),
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
    cvx <- data$CV
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
  if (is.na(datax$WildTypeDeath)) {datax$WildTypeDeath <- 0}
  if (is.na(datax$MutantDeath)) {datax$MutantDeath <- 0}
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
    cvy <- data$CV
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
  if (is.na(datay$WildTypeDeath)) {datay$WildTypeDeath <- 0}
  if (is.na(datay$MutantDeath)) {datay$MutantDeath <- 0}
  if (is.na(datay$Inoculum)) {datay$Inoculum <- 0}
  
  result <- tryCatch(lrt.mutations(datax = datax$CountsSelective, datay = datay$CountsSelective, ex = ex, ey = ey, 
                                   wx = datax$Fitness/(1 - datax$WildTypeDeath), wy = datay$Fitness/(1 - datay$WildTypeDeath),
                                   lagx = datax$Lag, lagy = datay$Lag, poissonx = datax$Residual, poissony = datay$Residual,
                                   deathx = datax$MutantDeath, deathy = datay$MutantDeath, phix = datax$Inoculum/Ntx, phiy = datay$Inoculum/Nty,
                                   cvx = cv0x, cvy = cv0y, Nx = Ntx, Ny = Nty, Mx = Mx, My = My, verbose = F),
                     error=function(err) {c(NA)})
  
  return(result)
}

