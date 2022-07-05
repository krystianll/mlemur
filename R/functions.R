### functions for calculating ML point and interval estimates of m as well as LRT P values ###

rluria <- function(n=10, rate=1e-8, N0=1, Nt=1e9, mut_fit=1.0, death_prob=0, lag=0, e=1, cv=0, trim=0, ret="list") {
  ret <- match.arg(ret, c("list", "vector"))
  if (ret=="list") {
    return(rluria_list(n=n, rate=rate, N0=N0, Nt=Nt, mut_fit=mut_fit, death_prob=death_prob, lag=lag, e=e, cv=cv, trim=trim))
  } else {
    if (cv!=0) warning("Return type is set to \'vector\' but cv is not 0. The information about total culture sizes will be lost.")
    return(rluria_vec(n=n, rate=rate, N0=N0, Nt=Nt, mut_fit=mut_fit, death_prob=death_prob, lag=lag, e=e, cv=cv, trim=trim))
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
      
      return(aux.seq(e = e, w = w, death = death, lag = lag, phi = phi, n = n, integrate = TRUE)/(1 - death/(1-death)))
      
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

# Approximate methods of m estimation
# P0, Jones median, Drake median, and Lea-Coulson median
# Most useful as an initial guess for m in Newton method
# Jones ME. Mutat Res Mutagen Relat Subj 1993 doi:10.1016/0165-1161(93)90146-Q
# Crane GJ, Thomas SM, Jones ME. Mutat Res - Fundam Mol Mech Mutagen 1996 doi:10.1016/0027-5107(96)00009-7
# Rosche WA, Foster PL. Methods 2000 doi:10.1006/meth.1999.0901
m_approx <- function(data, eff) {
  if (eff<=0) return(NA)
  if (median(data)==0) {
    m <- ifelse(eff==1, -log(sum(data==0)/length(data)), (1-eff)/eff/log(eff)*(log(sum(data==0)/length(data))))
  } else {
    m <- (median(data)/eff-0.693)/(log(median(data)/eff)+0.367)
  }
  return(m)
}

# p0 <- function(data, eff=1, lag=0) {
#   # unlist(lapply(test2-poisson, max, 0)); m/(1+poisson)
#   if (all(data==0)) return(NA)
#   m <- ifelse(eff==1, -log(sum(data==0)/length(data)), (1-eff)/eff/log(eff)*(log(sum(data==0)/length(data))))
#   m <- m*2^lag
#   return(m)
# }

p0 <- function(data, e=1, w=1, d=0, lag=0, phi=0, poisson=0) {
  if (all(data == 0)) return(NA)
  m <- log(sum(data == 0)/length(data))/aux.seq(e = e, w = w, death = d, lag = lag, phi = phi, n = 1)[1] - poisson
  return(m)
}

# drake.formula <- function(data, eff=1) {
#   mut <- data/eff
#   m <- rep(0, length(data))
#   for (i in 1:length(data)) {
#     m[i] <- (1 + sqrt(1+4*mut[i]))/2
#     for (j in 1:500) {
#       g <- m[i]/mut[i]-1/(log(m[i]))
#       dg <- 1/(mut[i])+(1/m[i])*(1/(log(m[i]))^2)
#       dx <- g/dg
#       if (abs(dx)/m[i]<1e-6) {
#         break
#       } else {
#         m[i] <- m[i]-dx
#       }
#     }
#   }
#   mest <- median(m)
#   return(mest)
# }

# lea.coulson.median <- function(data, eff=1) {
#   mut <- data/eff
#   m <- rep(0, length(data))
#   for (i in 1:length(data)) {
#     m[i] <- (1 + sqrt(1+4*mut[i]))/2
#     for (j in 1:500) {
#       g <- m[i]/(mut[i]-poisson)-1/(log(m[i]/L)+1.24)
#       dg <- 1/(mut[i]-poisson)+(1/m[i])*(1/(log(m[i]/L)+1.24)^2)
#       dx <- g/dg
#       if (abs(dx)/m[i]<1e-6) {
#         break
#       } else {
#         m[i] <- m[i]-dx
#       }
#     }
#   }
#   mest <- median(m)
#   return(mest)
# }

# The basic formula comes from a well-known paper Lea DE, Coulson CA. J Genet. 1949 doi:10.1007/BF02986080
# Corrections for partial plating, phenotypic delay and residual mutation after Angerer WP. Mutat Res - Fundam Mol Mech Mutagen. 2001 doi:10.1016/S0027-5107(01)00203-2
# Correction for differential growth taken from Eq. 12 in Koch AL. Mutat Res - Fundam Mol Mech Mutagen. 1982 doi:10.1016/0027-5107(82)90252-4
# Correction for mutant death based on the Angerer's rationale regarding phenotypic lag, with Newcombe correction for extra cellular divisions Newcombe HB. Genetics. 1948 doi:10.1093/genetics/33.5.447
lea.coulson.median <- function(data, e=1, w=1, poisson=0, lag=0, death=0, confint=FALSE) {
  L <- 2^lag
  d <- (1-death)/(1-2*death)
  
  mut <- median(data)/e/w-poisson
  if (mut<=0) {
    return(NA)
  } else {
    m <- (1 + sqrt(1+4*mut*d*L))/2
    for (j in 1:500) {
      # f <- m/(mut)-1/(log(m/L/d)+1.24)
      # df.dm <- 1/(mut)+(1/m)*(1/(log(m/L/d)+1.24)^2)
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

koch.quartiles <- function(data, lag = 0){
  q <- as.vector(quantile(data/2^lag, c(0.25, 0.50, 0.75)))
  m1 <- (1.7335 + 0.4474*q[1] - 0.002755*q[1]*q[1])*2^lag
  m2 <- (1.1580 + 0.2730*q[2] - 0.000761*q[2]*q[2])*2^lag
  m3 <- (0.6658 + 0.1497*q[3] - 0.0001387*q[3]*q[3])*2^lag
  print(c(mean(c(m1, m2, m3)), m1, m2, m3))
  
  q <- as.vector(quantile(data, c(0.25, 0.50, 0.75)))
  m1 <- 1.7335 + 0.4474*q[1] - 0.002755*q[1]*q[1]
  m2 <- 1.1580 + 0.2730*q[2] - 0.000761*q[2]*q[2]
  m3 <- 0.6658 + 0.1497*q[3] - 0.0001387*q[3]*q[3]
  return(c(mean(c(m1, m2, m3)), m1, m2, m3))
}

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

# death.corr <- function(m.obs, dr) {
#   z=dr/(1-dr)
#   # r <- (5*m.obs-2.236*sqrt(5*m.obs^2-6*m.obs^2*z))/(3*z)
#   r <- m.obs
#   for (i in 1:500) {
#     # print(r)
#     f.r <- log(1.0945*r^-0.0299)*r^3*z^2/m.obs^2 + log(0.6225*r^0.0551)*r^2*z/m.obs + r - m.obs
#     df.dr <- 1 + (0.0551*r^1.*z)/m.obs - (0.0299*r^2.*z^2)/m.obs^2 + (3*r^2*z^2*log(1.0945/r^0.0299))/m.obs^2 + (2*r*z*log(0.6225*r^0.0551))/m.obs
#     ratio <- f.r/df.dr
#     if (abs(ratio/r) < 1e-6) {return(r)}
#     else {r <- r-ratio}
#   }
# }
# 
# death.corr.2 <- function(o, d) {
#   z <- d/(1-d)
#   r <- suppressWarnings((5*o-2.236*sqrt(5*o^2-6*o^2*z))/(3*z))
#   if (!is.finite(r)) r <- o
#   for (i in 1:500) {
#     # print(r)
#     f.r <- -0.006859 -  0.100553/(9.434736 + r) + (-1.628842 + sqrt(1 + exp(d)) - 0.011604/(6.321826 + r))*log(o/r)
#     df.dr <- -(sqrt(1 + exp(d))/r) + 0.100553/(9.434736 + r)^2 + (2*(0.814421 + 0.005802/(6.321826 + r)))/r + (0.011604*log(o/r))/(6.321826 + r)^2
#     ratio <- f.r/df.dr
#     if (abs(ratio/r) < 1e-6) {return(r)}
#     else {r <- r-ratio}
#   }
# }

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
  # seq <- aux_seq_lag_pois_n_500[[as.integer(lag*10)]][1:n]
  # seq <- aux_seq_lag_gamma(lag = lag, n = n)
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

# lrt.mutations <- function(datax, datay, ex=1, ey=1, wx=1, wy=1, lagx=0, lagy=0, poissonx=0, poissony=0, deathx=0, deathy=0,
#                           phix=0, phiy=0, cvx=0, cvy=0, Nx=1, Ny=1, Mx=NA, My=NA, verbose=FALSE) {
#   
#   if (!is.finite(Mx)) {Mx <- mle.mutations(data = datax, e = ex, w=wx, lag=lagx, poisson=poissonx, death=deathx, phi=phix, cv = cvx, confint = FALSE, verbose=verbose)} else {Mx <- Mx}
#   if (!is.finite(My)) {My <- mle.mutations(data = datay, e = ey, w=wy, lag=lagy, poisson=poissony, death=deathy, phi=phiy, cv = cvy, confint = FALSE, verbose=verbose)} else {My <- My}
#   if (!is.finite(Mx) | !is.finite(My)) {if (verbose) message("Couldn't estimate mutation rates"); return(NA)}
#   if (verbose) message(paste("Estimated numbers of mutations are", Mx, "and", My))
#   if (verbose) message(paste("Estimated mutation rates are", Mx/Nx, "and", My/Ny))
#   
#   if (Nx>=Ny){
#     data1 <- datax; e1 <- ex; w1 <- wx; lag1 <- lagx; poisson1 <- poissonx; death1 <- deathx; phi1 <- phix; cv1 <- cvx; N1 <- Nx; M1 <- Mx
#     data2 <- datay; e2 <- ey; w2 <- wy; lag2 <- lagy; poisson2 <- poissony; death2 <- deathy; phi2 <- phiy; cv2 <- cvy; N2 <- Ny; M2 <- My
#   } else {
#     data1 <- datay; e1 <- ey; w1 <- wy; lag1 <- lagy; poisson1 <- poissony; death1 <- deathy; phi1 <- phiy; cv1 <- cvy; N1 <- Ny; M1 <- My
#     data2 <- datax; e2 <- ex; w2 <- wx; lag2 <- lagx; poisson2 <- poissonx; death2 <- deathx; phi2 <- phix; cv2 <- cvx; N2 <- Nx; M2 <- Mx
#   }
#   
#   # mutation rates
#   n1 <- max(data1)+1
#   n2 <- max(data2)+1
#   if (cv1>0) {k1 <- 1/cv1/cv1} else {k1 <- 0}
#   if (cv2>0) {k2 <- 1/cv2/cv2} else {k2 <- 0}
#   R <- N2/N1
#   
#   seq1 <- aux.seq(e = e1, w = w1, death = death1, lag = lag1, phi = phi1, n = n1-1)
#   seq2 <- aux.seq(e = e2, w = w2, death = death2, lag = lag2, phi = phi2, n = n2-1)
#   
#   # combined mutation rate
#   mg <- seq(min(M1,M2), max(M1,M2), length.out=10)
#   m1 <- -1
#   for (value in mg) {
#     if (m1>0) {break}
#     if (verbose) message(paste("m_combined: using guess", value))
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (verbose) message(paste("m_combined: iteration", iter))
#       if (iter>30) {if (verbose) message("m_combined: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- combo_score_fisher_ratio(m = m0, n_1 = n1, n_2 = n2, R = R, data_1 = data1, seq_1 = seq1, data_2 = data2, seq_2 = seq2, k_1 = k1, k_2 = k2, poisson_1 = poisson1, poisson_2 = poisson2)
#       if (!is.finite(score.fisher.ratio)) {if (verbose) message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {if (verbose) message("m_combined: negative estimate"); m1 <- -1; break}
#       if (abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   if (m1<0) if (verbose) message("m_combined: trying MPFR")
#   # MPFR
#   for (value in mg) {
#     if (m1>0) {break}
#     if (verbose) message(paste("m_combined: using guess", value))
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (verbose) message(paste("m_combined: iteration", iter))
#       if (iter>30) {if (verbose) message("m_combined: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- combo_score_fisher_ratio_boost(xm = m0, n_1 = n1, n_2 = n2, xR = R, data_1 = data1, seq_1 = seq1, data_2 = data2, seq_2 = seq2, xk_1 = k1, xk_2 = k2, xpoisson_1 = poisson1, xpoisson_2 = poisson2)
#       if (!is.finite(score.fisher.ratio)) {if (verbose) message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {if (verbose) message("m_combined: negative estimate"); m1 <- -1; break}
#       if (abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   if (m1<0) {
#     mg <- seq(min(M1,M2), max(M1,M2,M2/R), length.out=200)
#     if (verbose) message("m_combined: no success, looping")
#     # for (i in 0:100) {
#     for (value in mg) {
#       if (m1>0) {break}
#       # m0 <- min(M1,M2)+i*(abs(M1-M2)/100)
#       m0 <- value
#       if (verbose) message(paste("m_combined: using guess", m0))
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (verbose) message(paste("m_combined: iteration", iter))
#         if (iter>30) {if (verbose) message("m_combined: reached iteration limit"); m1 <- -1; break}
#         score.fisher.ratio <- combo_score_fisher_ratio(m = m0, n_1 = n1, n_2 = n2, R = R, data_1 = data1, seq_1 = seq1, data_2 = data2, seq_2 = seq2, k_1 = k1, k_2 = k2, poisson_1 = poisson1, poisson_2 = poisson2)
#         if (!is.finite(score.fisher.ratio)) {score.fisher.ratio <- combo_score_fisher_ratio_boost(xm = m0, n_1 = n1, n_2 = n2, xR = R, data_1 = data1, seq_1 = seq1, data_2 = data2, seq_2 = seq2, xk_1 = k1, xk_2 = k2, xpoisson_1 = poisson1, xpoisson_2 = poisson2)}
#         if (!is.finite(score.fisher.ratio)) {if (verbose) message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#         m1 <- m0+score.fisher.ratio
#         if (m1 <= 0) {if (verbose) message("m_combined: negative estimate"); m1 <- -1; break}
#         if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#       }
#     }
#   }
#   if (m1 < 0) {if (verbose) message("m_combined: couldn't estimate m"); return(NA)}
#   
#   MC <- m1
#   
#   if (verbose) message(paste("m_combined: obtained estimate", MC))
#   
#   # log likelihood functions
#   la <- logprob(m = MC, len = n1, data = data1, seq = seq1, k = k1, poisson = poisson1)
#   if (!is.finite(la)) {la <- logprob_boost(xm = MC, len = n1, data = data1, seq = seq1, xk = k1, xpoisson = poisson1)}
#   lb <- logprob(m = R*MC, len = n2, data = data2, seq = seq2, k = k2, poisson = poisson2)
#   if (!is.finite(lb)) {lb <- logprob_boost(xm = R*MC, len = n2, data = data2, seq = seq2, xk = k2, xpoisson = poisson2)}
#   lc <- logprob(m = M1, len = n1, data = data1, seq = seq1, k = k1, poisson = poisson1)
#   if (!is.finite(lc)) {lc <- logprob_boost(xm = M1, len = n1, data = data1, seq = seq1, xk = k1, xpoisson = poisson1)}
#   ld <- logprob(m = M2, len = n2, data = data2, seq = seq2, k = k2, poisson = poisson2)
#   if (!is.finite(ld)) {ld <- logprob_boost(xm = M2, len = n2, data = data2, seq = seq2, xk = k2, xpoisson = poisson2)}
#   if(!(all(is.finite(c(la,lb,lc,ld))))) {message("loglikelihood: non-numeric value encountered"); return(NA)}
#   l0 <- la+lb
#   l1 <- lc+ld
#   chi <- -2*(l0-l1)
#   pval <- 1-pchisq(chi,1)
#   return(pval)
# }

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
  # pval <- 1-pchisq(chi,1)
  pval <- pchisq(chi,1,lower.tail = F)
  return(pval)
}

# ### Calculate MLE, CI, and P value under Luria-Delbruck model ###
# # Combined functions for both Lea-Coulson (phi=1, w=1) and Mandelbrot-Koch (phi=1, w!=1) distributions
# # based on functions newton.LD, newton.MK, newton.LD.plating, confint.LD, confint.MK, confint.LD.plating from rSalvador 1.8
# # simplified and optimised to avoid redundant computations, partial plating is possible in all cases
# # arbitrary precision arithmetics has been implemented and trial-and-error iterative initial value of m approach is used when necessary to avoid under- and overflow
# # Zheng Q. Math Biosci 2002 doi:10.1016/S0025-5564(02)00087-1
# # Zheng Q. Math Biosci 2005 doi:10.1016/j.mbs.2005.03.011
# # Zheng Q. Math Biosci 2008 doi:10.1016/j.mbs.2008.09.002
# # Zheng Q. Genetica 2016 doi:10.1007/s10709-016-9904-3
# mle.ld <- function(data, e=1, w=1, lag=0, poisson=0, death=0, confint=TRUE, verbose=FALSE) {
#   n <- max(data)+1
#   
#   r <- 1/w
#   seq <- aux.seq(e = e, w = w, death = death, lag = lag, phi = 0, n = n-1)
#   # if (lag!=0) {
#   #   if (w != 1) message("Differential growth rate cannot be used together with phenotypic lag correction.")
#   #   L <- 2^lag
#   #   seq <- aux_seq_lag_ext(L, e, n-1)
#   # } else if (death!=0) {
#   #   # if (e != 1) message("Partial plating cannot be used together with death correction.")
#   #   # seq <- aux_seq_death(death, w, n-1)
#   #   seq <- aux_seq_death(e, w, death/(1-death), n-1)
#   # } else {
#   #   seq <- aux.seq(e = e, w = w, n = n-1)
#   # }
#   
#   if (length(seq)<2) {if (verbose) message("Failed to get the sequence"); return(NA)}
#   
#   # mg1 <- m_approx(data, e*w)
#   # mg2 <- m_approx((min(data)+1), e*w)
#   # if (w!=1) {mg3 <- ifelse(w<1, m_approx(data, e)*(1.2*log(r)+1.2), m_approx(data, e)*(1.12*r^2-0.14*r+0.02))} else {mg3 <- NA}
#   # if (w!=1) {mg4 <- ifelse(w<1, m_approx((min(data)+1), e)*(1.2*log(r)+1.2), m_approx((min(data)+1), e)*(1.12*r^2-0.14*r+0.02))} else {mg4 <- NA}
#   # if (poisson!=0) {data2 <- data-poisson; data2[which(data2<0)] <- 0; mg5 <- m_approx(data2, e)} else {mg5 <- NA}
#   # mg <- as.vector(na.exclude(c(mg5, mg5*0.5, mg1, mg1*0.5, mg2, mg2*0.5, mg3, mg3*0.5, mg4, mg4*0.5)))
#   # mg <- mg[mg<1e6]
#   # mg <- mg[mg>0]
#   mg <- m_guess(data = data, e = e, w = w, lag = lag, death = death, poisson = poisson)
#   mg <- c(mg, 0.5*mg, 0.7*mg, 1.3*mg, 0.1*mg, 0.001, 1)
#   # print(mg)
#   m1 <- -1
#   
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter > 30) {if (verbose) message("m: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- score_fisher_ratio(m0, seq, n, data, poisson = poisson)
#       # stop()
#       if (is.finite(score.fisher.ratio)==FALSE) {if (verbose) message("m: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   # MPFR
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter > 30) {if (verbose) message("m: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- score_fisher_ratio(m0, seq, n, data, poisson = poisson)
#       if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- score_fisher_ratio_boost(m0, seq, n, data, xpoisson = poisson)}
#       if (is.finite(score.fisher.ratio)==FALSE) {if (verbose) message("m: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   if (m1 < 0) {
#     if (verbose) message("no success, looping")
#     mdown <- m_approx(max(1, min(data)), e*w)*0.1
#     if (w <= 1) {mtop <- m_approx(max(data), e*w)/0.5} 
#     else {mtop <- m_approx(max(data), e*w)^(1/w)/0.5}
#     for (i in 0:100) {
#       if (m1>0) {break}
#       if (verbose) message(paste("loop", i))
#       m0 <- mdown+i*(abs(mtop-mdown)/100)
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {if (verbose) message("m: reached iteration limit"); m1 <- -1; break}
#         score.fisher.ratio <- score_fisher_ratio(m0, seq, n, data, poisson = poisson)
#         if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- score_fisher_ratio_boost(m0, seq, n, data, xpoisson = poisson)}
#         if (is.finite(score.fisher.ratio)==FALSE) {if (verbose) message("m: non-numeric value encountered"); m1 <- -1; break}
#         m1 <- m0+score.fisher.ratio
#         if (m1 <= 0) {m1 <- -1; break}
#         if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#       }
#     }
#   }
#   
#   if (m1 < 0) {if (verbose) message("Couldn't estimate m"); if (confint==TRUE) {return(c(NA,NA,NA))} else {return(c(NA))}}
#   if (confint==FALSE) {
#     return(m1)
#   } else {
#     
#     # confidence intervals
#     chi <- qchisq(1-0.05,1)
#     root.logprobtail <- root_logprobtail(m1, seq, n, data, chi, poisson = poisson)
#     if (all(is.finite(root.logprobtail))==FALSE) {root.logprobtail <- root_logprobtail_boost(m1, seq, n, data, chi, xpoisson = poisson)}
#     if (all(is.finite(root.logprobtail))==FALSE) {if (verbose) message("Couldn't set the boundaries for CI"); return(c(m1,NA,NA))}
#     root <- root.logprobtail[1]
#     logprobtail <- root.logprobtail[2]
#     
#     # lower interval
#     mg.low <- as.vector(na.exclude(c(m1-0.5*root, m1-0.9*root, (m1-0.5*root)*0.5)))
#     m1.low <- -1
#     for (value in mg.low) {
#       if (m1.low>0) {break}
#       m0.low <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {if (verbose) message("m_low: reached iteration limit"); m1.low <- -1; break}
#         score.fisher.ratio.low <- log_score_ratio(m0.low, seq, n, data, logprobtail, poisson = poisson)
#         if (is.finite(score.fisher.ratio.low)==FALSE) {if (verbose) message("m_low: non-numeric value encountered"); m1.low <- -1; break}
#         m1.low <- m0.low-score.fisher.ratio.low
#         if (m1.low <= 0) {m1.low <- -1; break}
#         if(abs(score.fisher.ratio.low)/m0.low < 1e-6) {if(m1.low>m1) {m1.low <- -1; break} else {break}} else {m0.low <- m1.low}
#       }
#     }
#     # MPFR
#     for (value in mg.low) {
#       if (m1.low>0) {break}
#       m0.low <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {if (verbose) message("m_low: reached iteration limit"); m1.low <- -1; break}
#         score.fisher.ratio.low <- log_score_ratio(m0.low, seq, n, data, logprobtail, poisson = poisson)
#         if (is.finite(score.fisher.ratio.low)==FALSE) {score.fisher.ratio.low <- log_score_ratio_boost(m0.low, seq, n, data, logprobtail, xpoisson = poisson)}
#         if (is.finite(score.fisher.ratio.low)==FALSE) {if (verbose) message("m_low: non-numeric value encountered"); m1.low <- -1; break}
#         m1.low <- m0.low-score.fisher.ratio.low
#         if (m1.low <= 0) {m1.low <- -1; break}
#         if(abs(score.fisher.ratio.low)/m0.low < 1e-6) {if(m1.low>m1) {m1.low <- -1; break} else {break}} else {m0.low <- m1.low}
#       }
#     }
#     if (m1.low < 0) {
#       for (i in 0:100) {
#         if (m1.low>0) {break}
#         if (verbose) message(paste("loop", i))
#         m0.low <- (1+i)*m1*0.01
#         iter <- 0
#         repeat {
#           iter <- iter+1
#           if (iter>30) {if (verbose) message("m_low: reached iteration limit"); m1.low <- -1; break}
#           score.fisher.ratio.low <- log_score_ratio(m0.low, seq, n, data, logprobtail, poisson = poisson)
#           if (is.finite(score.fisher.ratio.low)==FALSE) {score.fisher.ratio.low <- log_score_ratio_boost(m0.low, seq, n, data, logprobtail, xpoisson = poisson)}
#           if (is.finite(score.fisher.ratio.low)==FALSE) {if (verbose) message("m: non-numeric value encountered"); m1.low <- -1; break}
#           m1.low <- m0.low-score.fisher.ratio.low
#           if (m1.low <= 0) {m1.low <- -1; break}
#           if(abs(score.fisher.ratio.low)/m0.low < 1e-6) {if(m1.low>m1) {m1.low <- -1; break} else {break}} else {m0.low <- m1.low}
#         }
#       }
#     }
#     if (m1.low < 0) {if (verbose) message("Couldn't estimate m_low"); m1.low <- NA}
#     
#     # upper interval
#     mg.up <- as.vector(na.exclude(c(m1+0.5*root, m1+0.9*root, (m1+0.5*root)*0.5)))
#     m1.up <- -1
#     for (value in mg.up) {
#       if (m1.up>0) {break}
#       m0.up <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {if (verbose) message("m_up: reached iteration limit"); m1.up <- -1; break}
#         score.fisher.ratio.up <- log_score_ratio(m0.up, seq, n, data, logprobtail, poisson = poisson)
#         if (is.finite(score.fisher.ratio.up)==FALSE) {if (verbose) message("m_up: non-numeric value encountered"); m1.up <- -1; break}
#         m1.up <- m0.up-score.fisher.ratio.up
#         if (m1.up <= 0) {m1.up <- -1; break}
#         if(abs(score.fisher.ratio.up)/m0.up < 1e-6) {if(m1.up<m1) {m1.up <- -1; break} else {break}} else {m0.up <- m1.up}
#       }
#     }
#     # MPFR
#     for (value in mg.up) {
#       if (m1.up>0) {break}
#       m0.up <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {if (verbose) message("m_up: reached iteration limit"); m1.up <- -1; break}
#         score.fisher.ratio.up <- log_score_ratio(m0.up, seq, n, data, logprobtail, poisson = poisson)
#         if (is.finite(score.fisher.ratio.up)==FALSE) {score.fisher.ratio.up <- log_score_ratio_boost(m0.up, seq, n, data, logprobtail, xpoisson = poisson)}
#         if (is.finite(score.fisher.ratio.up)==FALSE) {if (verbose) message("m_up: non-numeric value encountered"); m1.up <- -1; break}
#         m1.up <- m0.up-score.fisher.ratio.up
#         if (m1.up <= 0) {m1.up <- -1; break}
#         if(abs(score.fisher.ratio.up)/m0.up < 1e-6) {if(m1.up<m1) {m1.up <- -1; break} else {break}} else {m0.up <- m1.up}
#       }
#     }
#     if (m1.up < 0) {
#       for (i in 0:100) {
#         if (m1.up>0) {break}
#         if (verbose) message(paste("loop", i))
#         m0.up <- (1+i/10)*m1
#         iter <- 0
#         repeat {
#           iter <- iter+1
#           if (iter>30) {if (verbose) message("m_up: reached iteration limit"); m1.up <- -1; break}
#           score.fisher.ratio.up <- log_score_ratio_ld(m0.up, seq, n, data, logprobtail, poisson = poisson)
#           if (is.finite(score.fisher.ratio.up)==FALSE) {score.fisher.ratio.up <- log_score_ratio_ld_boost(m0.up, seq, n, data, logprobtail, xpoisson = poisson)}
#           if (is.finite(score.fisher.ratio.up)==FALSE) {if (verbose) message("m: non-numeric value encountered"); m1.up <- -1; break}
#           m1.up <- m0.up-score.fisher.ratio.up
#           if (m1.up <= 0) {m1.up <- -1; break}
#           if(abs(score.fisher.ratio.up)/m0.up < 1e-6) {if(m1.up<m1) {m1.up <- -1; break} else {break}} else {m0.up <- m1.up}
#         }
#       }
#     }
#     if (m1.up < 0) {if (verbose) message("Couldn't estimate m_up"); m1.up <- NA}
#     
#     return(c(m1, m1.low, m1.up))
#   }
#   
# }
# 
# lrt.ld <- function(datax, datay, ex=1, ey=1, wx=1, wy=1, Nx=1, Ny=1, Mx=NA, My=NA) {
#   
#   if (is.finite(Mx)==FALSE) {Mx <- mle.ld(datax, ex, wx, confint = FALSE)}
#   if (is.finite(My)==FALSE) {My <- mle.ld(datay, ey, wy, confint = FALSE)}
#   if (is.finite(Mx)==FALSE | is.finite(My)==FALSE) {message("Couldn't estimate mutation rates"); return(c(NA,NA))}
#   
#   if (Nx>Ny){
#     data1 <- datax
#     data2 <- datay
#     e1 <- ex
#     e2 <- ey
#     w1 <- wx
#     w2 <- wy
#     N1 <- Nx
#     N2 <- Ny
#     M1 <- Mx
#     M2 <- My
#   } else {
#     data1 <- datay
#     data2 <- datax
#     e1 <- ey
#     e2 <- ex
#     w1 <- wy
#     w2 <- wx
#     N1 <- Ny
#     N2 <- Nx
#     M1 <- My
#     M2 <- Mx
#   }
#   
#   n1 <- max(data1)+1
#   n2 <- max(data2)+1
#   R <- N2/N1
#   dataC <- c(data1, data2)
#   eC <- (e1*length(data1)+e2*length(data2))/(length(dataC))
#   wC <- (w1*length(data1)+w2*length(data2))/(length(dataC))
#   seq1 <- aux.seq(e = e1, w = w1, n = n1-1)
#   seq2 <- aux.seq(e = e2, w = w2, n = n2-1)
#   
#   # combined mutation rate
#   mg1 <- m_approx(dataC, eC*wC)
#   mg2 <- m_approx((min(dataC)+1), eC*wC)
#   if (wC!=1) {mg3 <- ifelse(wC<1, m_approx(dataC, eC)*(1.2*log(1/wC)+1.2), m_approx(dataC, eC)*(1.12*(1/wC)^2-0.14*(1/wC)+0.02))} else {mg3 <- NA}
#   if (wC!=1) {mg4 <- ifelse(wC<1, m_approx((min(dataC)+1), eC)*(1.2*log(1/wC)+1.2), m_approx((min(dataC)+1), eC)*(1.12*(1/wC)^2-0.14*(1/wC)+0.02))} else {mg4 <- NA}
#   mg <- as.vector(na.exclude(c(mg1, mg1*0.5, mg2, mg2*0.5, mg3, mg3*0.5, mg4, mg4*0.5, min(M1,M2), min(M1,M2)*0.5)))
#   mg <- mg[mg<1e6]
#   m1 <- -1
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter>30) {message("m_combined: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- combo_score_fisher_ratio(m0, n1, n2, R, data1, seq1, data2, seq2)
#       if (is.finite(score.fisher.ratio)==FALSE) {message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if (abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   # MPFR
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter>30) {message("m_combined: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- combo_score_fisher_ratio(m0, n1, n2, R, data1, seq1, data2, seq2)
#       if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- combo_score_fisher_ratio_boost(m0, n1, n2, R, data1, seq1, data2, seq2)}
#       if (is.finite(score.fisher.ratio)==FALSE) {message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if (abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   if (m1<0) {
#     for (i in 0:100) {
#       if (m1>0) {break}
#       message(paste("loop", i))
#       m0 <- min(M1,M2)+i*(abs(M1-M2)/100)
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m_combined: reached iteration limit #loop)"); m1 <- -1; break}
#         score.fisher.ratio <- combo_score_fisher_ratio(m0, n1, n2, R, data1, seq1, data2, seq2)
#         if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- combo_score_fisher_ratio_boost(m0, n1, n2, R, data1, seq1, data2, seq2)}
#         if (is.finite(score.fisher.ratio)==FALSE) {message("m_combined: non-numeric value encountered #loop)"); m1 <- -1; break}
#         m1 <- m0+score.fisher.ratio
#         if (m1 <= 0) {m1 <- -1; break}
#         if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#       }
#     }
#   }
#   if (m1 < 0) {message("Couldn't estimate combined m"); return(c(NA,NA))}
#   MC <- m1
#   
#   # log likelihood functions
#   a <- logprob(MC, n1, data1, seq1)
#   if (is.finite(a)==FALSE) {a <- logprob_boost(MC, n1, data1, seq1)}
#   b <- logprob(R*MC, n2, data2, seq2)
#   if (is.finite(b)==FALSE) {b <- logprob_boost(R*MC, n2, data2, seq2)}
#   c <- logprob(M1, n1, data1, seq1)
#   if (is.finite(c)==FALSE) {c <- logprob_boost(M1, n1, data1, seq1)}
#   d <- logprob(M2, n2, data2, seq2)
#   if (is.finite(d)==FALSE) {d <- logprob_boost(M2, n2, data2, seq2)}
#   if(!(all(is.finite(c(a,b,c,d))))) {message("loglikelihood: non-numeric value encountered"); return(c(NA, NA))}
#   l0 <- a+b
#   l1 <- c+d
#   chi <- -2*(l0-l1)
#   pval <- 1-pchisq(chi,1)
#   effsize <- chi/(-2*l0)
#   return(c(pval, effsize))
#   
# }
# 
# ### Calculate MLE, CI, and P value under the Bartlett model ###
# # Combined functions for B0 (gamma-mixture) distribution indexed by e and cv
# # based on functions newton.B0 and confint.B0 from rSalvador 1.8
# # P values can be calculated using Likelihood Ratio Test, using an established framework
# # simplified and optimised to avoid redundant computations, partial plating is possible in all cases
# # arbitrary precision arithmetics has been implemented and trial-and-error iterative initial value of m approach is used when necessary to avoid under- and overflow
# # Zheng Q. Math Biosci 2008 doi:10.1016/j.mbs.2008.09.002
# # Zheng Q. Statistics (Ber) 2010 doi:10.1080/02331880903236868
# # Zheng Q. Genetica 2016 doi:10.1007/s10709-016-9904-3
# mle.b0 <- function(data, e=1, w=1, cv=1e-3, confint=TRUE) {
#   n <- max(data)+1
#   k <- 1/cv/cv
#   
#   seq <- aux.seq(e = e, w = w, n = n-1)
#   
#   # mutants in the culture
#   mg1 <- m_approx(data, e)
#   mg2 <- m_approx((min(data)+1), e)
#   mg <- c(mg1, mg2, mg1*0.5)
#   m1 <- -1
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter > 30) {message("m: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- score_fisher_ratio(m0, seq, n, data, k)
#       # stop()
#       if (is.finite(score.fisher.ratio)==FALSE) {message("m: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   # MPFR
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter > 30) {message("m: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- score_fisher_ratio(m0, seq, n, data, k)
#       if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- score_fisher_ratio_boost(m0, seq, n, data, k)}
#       if (is.finite(score.fisher.ratio)==FALSE) {message("m: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   if (m1 < 0) {
#     message("no success, looping")
#     mdown <- jones.median(max(1, min(data)), e)*0.5
#     mtop <- jones.median(max(data), e)/0.5
#     for (i in 0:100) {
#       if (m1>0) {break}
#       message(paste("loop", i))
#       m0 <- mdown+i*(abs(mtop-mdown)/100)
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m: reached iteration limit"); m1 <- -1; break}
#         score.fisher.ratio <- score_fisher_ratio(m0, seq, n, data, k)
#         if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- score_fisher_ratio_boost(m0, seq, n, data, k)}
#         if (is.finite(score.fisher.ratio)==FALSE) {message("m: non-numeric value encountered"); m1 <- -1; break}
#         m1 <- m0+score.fisher.ratio
#         if (m1 <= 0) {m1 <- -1; break}
#         if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#       }
#     }
#   }
#   if (m1 < 0) {message("Couldn't estimate m"); if (confint==TRUE) {return(c(NA,NA,NA))} else {return(c(NA))}}
#   if (confint==FALSE) {
#     return(m1)
#   } else {
#     
#     # confidence intervals
#     chi <- qchisq(1-0.05,1)
#     root.logprobtail <- root_logprobtail(m1, seq, n, data, chi, k)
#     if (all(is.finite(root.logprobtail))==FALSE) {root.logprobtail <- root_logprobtail_boost(m1, seq, n, data, chi, k)}
#     if (all(is.finite(root.logprobtail))==FALSE) {message("Couldn't set the boundaries for CI"); return(c(m1,NA,NA))}
#     root <- root.logprobtail[1]
#     logprobtail <- root.logprobtail[2]
#     
#     # lower interval
#     mg.low <- as.vector(na.exclude(c(m1-0.5*root, m1-0.9*root, (m1-0.5*root)*0.5)))
#     m1.low <- -1
#     for (value in mg.low) {
#       if (m1.low>0) {break}
#       m0.low <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m_low: reached iteration limit"); m1.low <- -1; break}
#         score.fisher.ratio.low <- log_score_ratio(m0.low, seq, n, data, logprobtail, k)
#         if (is.finite(score.fisher.ratio.low)==FALSE) {message("m_low: non-numeric value encountered"); m1.low <- -1; break}
#         m1.low <- m0.low-score.fisher.ratio.low
#         if (m1.low <= 0) {m1.low <- -1; break}
#         if(abs(score.fisher.ratio.low)/m0.low < 1e-6) {if(m1.low>m1) {m1.low <- -1; break} else {break}} else {m0.low <- m1.low}
#       }
#     }
#     # MPFR
#     for (value in mg.low) {
#       if (m1.low>0) {break}
#       m0.low <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m_low: reached iteration limit"); m1.low <- -1; break}
#         score.fisher.ratio.low <- log_score_ratio(m0.low, seq, n, data, logprobtail, k)
#         if (is.finite(score.fisher.ratio.low)==FALSE) {score.fisher.ratio.low <- log_score_ratio_boost(m0.low, seq, n, data, logprobtail, k)}
#         if (is.finite(score.fisher.ratio.low)==FALSE) {message("m_low: non-numeric value encountered"); m1.low <- -1; break}
#         m1.low <- m0.low-score.fisher.ratio.low
#         if (m1.low <= 0) {m1.low <- -1; break}
#         if(abs(score.fisher.ratio.low)/m0.low < 1e-6) {if(m1.low>m1) {m1.low <- -1; break} else {break}} else {m0.low <- m1.low}
#       }
#     }
#     if (m1.low < 0) {
#       for (i in 0:100) {
#         if (m1.low>0) {break}
#         message(paste("loop", i))
#         m0.low <- (1+i)*m1*0.01/(ifelse(cv<=1, 1, log(cv)))
#         iter <- 0
#         repeat {
#           iter <- iter+1
#           if (iter>30) {message("m_low: reached iteration limit"); m1.low <- -1; break}
#           score.fisher.ratio.low <- log_score_ratio(m0.low, seq, n, data, logprobtail, k)
#           if (is.finite(score.fisher.ratio.low)==FALSE) {score.fisher.ratio.low <- log_score_ratio_boost(m0.low, seq, n, data, logprobtail, k)}
#           if (is.finite(score.fisher.ratio.low)==FALSE) {message("m_low: non-numeric value encountered"); m1.low <- -1; break}
#           m1.low <- m0.low-score.fisher.ratio.low
#           if (m1.low <= 0) {m1.low <- -1; break}
#           if(abs(score.fisher.ratio.low)/m0.low < 1e-6) {if(m1.low>m1) {m1.low <- -1; break} else {break}} else {m0.low <- m1.low}
#         }
#       }
#     }
#     if (m1.low < 0) {message("Couldn't estimate m_low"); m1.low <- NA}
#     
#     # upper interval
#     mg.up <- as.vector(na.exclude(c(m1+0.5*root, m1+0.9*root, (m1+0.5*root)*0.5, (m1+0.5*root)/0.5)))
#     m1.up <- -1
#     for (value in mg.up) {
#       if (m1.up>0) {break}
#       m0.up <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m_up: reached iteration limit"); m1.up <- -1; break}
#         score.fisher.ratio.up <- log_score_ratio(m0.up, seq, n, data, logprobtail, k)
#         if (is.finite(score.fisher.ratio.up)==FALSE) {message("m_up: non-numeric value encountered"); m1.up <- -1; break}
#         m1.up <- m0.up-score.fisher.ratio.up
#         if (m1.up <= 0) {m1.up <- -1; break}
#         if(abs(score.fisher.ratio.up)/m0.up < 1e-6) {if(m1.up<m1) {m1.up <- -1; break} else {break}} else {m0.up <- m1.up}
#       }
#     }
#     # MPFR
#     for (value in mg.up) {
#       if (m1.up>0) {break}
#       m0.up <- value
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m_up: reached iteration limit"); m1.up <- -1; break}
#         score.fisher.ratio.up <- log_score_ratio(m0.up, seq, n, data, logprobtail, k)
#         if (is.finite(score.fisher.ratio.up)==FALSE) {score.fisher.ratio.up <- log_score_ratio_boost(m0.up, seq, n, data, logprobtail, k)}
#         if (is.finite(score.fisher.ratio.up)==FALSE) {message("m_up: non-numeric value encountered"); m1.up <- -1; break}
#         m1.up <- m0.up-score.fisher.ratio.up
#         if (m1.up <= 0) {m1.up <- -1; break}
#         if(abs(score.fisher.ratio.up)/m0.up < 1e-6) {if(m1.up<m1) {m1.up <- -1; break} else {break}} else {m0.up <- m1.up}
#       }
#     }
#     if (m1.up < 0) {
#       for (i in 0:100) {
#         if (m1.up>0) {break}
#         message(paste("loop", i))
#         m0.up <- (1+i)*m1*exp(cv)*ifelse(cv<=1, 0.1, 1)
#         iter <- 0
#         repeat {
#           iter <- iter+1
#           if (iter>30) {message("m_up: reached iteration limit"); m1.up <- -1; break}
#           score.fisher.ratio.up <- log_score_ratio(m0.up, seq, n, data, logprobtail, k)
#           if (is.finite(score.fisher.ratio.up)==FALSE) {score.fisher.ratio.up <- log_score_ratio_boost(m0.up, seq, n, data, logprobtail, k)}
#           if (is.finite(score.fisher.ratio.up)==FALSE) {message("m_up: non-numeric value encountered"); m1.up <- -1; break}
#           m1.up <- m0.up-score.fisher.ratio.up
#           if (m1.up <= 0) {m1.up <- -1; break}
#           if(abs(score.fisher.ratio.up)/m0.up < 1e-6) {if(m1.up<m1) {m1.up <- -1; break} else {break}} else {m0.up <- m1.up}
#         }
#       }
#     }
#     if (m1.up < 0) {message("Couldn't estimate m_up"); m1.up <- NA}
#     
#     return(c(m1, m1.low, m1.up))
#   }
# }
# 
# lrt.b0 <- function(datax, datay, ex=1, ey=1, cvx=1e-3, cvy=1e-3, Nx=1, Ny=1, Mx=NA, My=NA) {
#   
#   if (is.finite(Mx)==FALSE) {Mx <- mle.b0(data = datax, e = ex, cv = cvx, confint = FALSE)} else {Mx <- Mx}
#   if (is.finite(My)==FALSE) {My <- mle.b0(data = datay, e = ey, cv = cvy, confint = FALSE)} else {My <- My}
#   if (is.finite(Mx)==FALSE | is.finite(My)==FALSE) {message("Couldn't estimate mutation rates"); return(c(NA,NA))}
#   
#   if (Nx>=Ny){
#     data1 <- datax
#     data2 <- datay
#     e1 <- ex
#     e2 <- ey
#     cv1 <- cvx
#     cv2 <- cvy
#     N1 <- Nx
#     N2 <- Ny
#     M1 <- Mx
#     M2 <- My
#   } else {
#     data1 <- datay
#     data2 <- datax
#     e1 <- ey
#     e2 <- ex
#     cv1 <- cvy
#     cv2 <- cvx
#     N1 <- Ny
#     N2 <- Nx
#     M1 <- My
#     M2 <- Mx
#   }
#   
#   # mutation rates
#   n1 <- max(data1)+1
#   n2 <- max(data2)+1
#   dataC <- c(data1, data2)
#   eC <- (e1*length(data1)+e2*length(data2))/(length(dataC))
#   k1 <- 1/cv1/cv1
#   k2 <- 1/cv2/cv2
#   R <- N2/N1
#   
#   seq1 <- aux.seq(e = e1, w = 1, n = n1-1)
#   seq2 <- aux.seq(e = e2, w = 1, n = n2-1)
#   
#   # combined mutation rate
#   mg1 <- m_approx(dataC, eC)
#   mg2 <- m_approx((min(dataC)+1), eC)
#   mg <- as.vector(na.exclude(c(mg1, mg2, mg1*0.5, mg2*0.5, min(M1, M2), min(M1, M2)*0.5, min(M1,M2)/R)))
#   m1 <- -1
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter>30) {message("m_combined: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- combo_score_fisher_ratio(m0, n1, n2, R, data1, seq1, data2, seq2, k1, k2)
#       if (is.finite(score.fisher.ratio)==FALSE) {message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if (abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   # MPFR
#   for (value in mg) {
#     if (m1>0) {break}
#     m0 <- value
#     iter <- 0
#     repeat {
#       iter <- iter+1
#       if (iter>30) {message("m_combined: reached iteration limit"); m1 <- -1; break}
#       score.fisher.ratio <- combo_score_fisher_ratio(m0, n1, n2, R, data1, seq1, data2, seq2, k1, k2)
#       if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- combo_score_fisher_ratio_boost(m0, n1, n2, R, data1, seq1, data2, seq2, k1, k2)}
#       if (is.finite(score.fisher.ratio)==FALSE) {message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#       m1 <- m0+score.fisher.ratio
#       if (m1 <= 0) {m1 <- -1; break}
#       if (abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#     }
#   }
#   if (m1<0) {
#     for (i in 0:100) {
#       if (m1>0) {break}
#       message(paste("loop ", i))
#       m0 <- min(M1,M2)+i*(abs(M1-M2)/100)
#       iter <- 0
#       repeat {
#         iter <- iter+1
#         if (iter>30) {message("m_combined: reached iteration limit"); m1 <- -1; break}
#         score.fisher.ratio <- combo_score_fisher_ratio(m0, n1, n2, R, data1, seq1, data2, seq2, k1, k2)
#         if (is.finite(score.fisher.ratio)==FALSE) {score.fisher.ratio <- combo_score_fisher_ratio_boost(m0, n1, n2, R, data1, seq1, data2, seq2, k1, k2)}
#         if (is.finite(score.fisher.ratio)==FALSE) {message("m_combined: non-numeric value encountered"); m1 <- -1; break}
#         m1 <- m0+score.fisher.ratio
#         if (m1 <= 0) {m1 <- -1; break}
#         if(abs(score.fisher.ratio)/m0 < 1e-6) {break} else {m0 <- m1}
#       }
#     }
#   }
#   if (m1 < 0) {message("m_combined: no convergence (negative estimate)"); return(c(NA,NA))}
#   
#   MC <- m1
#   
#   # log likelihood functions
#   la <- logprob(MC, n1, data1, seq1, k1)
#   if (is.finite(la)==FALSE) {la <- logprob_boost(MC, n1, data1, seq1, k1)}
#   lb <- logprob(R*MC, n2, data2, seq2, k2)
#   if (is.finite(lb)==FALSE) {lb <- logprob_boost(R*MC, n2, data2, seq2, k2)}
#   lc <- logprob(M1, n1, data1, seq1, k1)
#   if (is.finite(lc)==FALSE) {lc <- logprob_boost(M1, n1, data1, seq1, k1)}
#   ld <- logprob(M2, n2, data2, seq2, k2)
#   if (is.finite(ld)==FALSE) {ld <- logprob_boost(M2, n2, data2, seq2, k2)}
#   if(!(all(is.finite(c(la,lb,lc,ld))))) {message("loglikelihood: non-numeric value encountered"); return(c(NA, NA))}
#   l0 <- la+lb
#   l1 <- lc+ld
#   chi <- -2*(l0-l1)
#   pval <- 1-pchisq(chi,1)
#   effsize <- chi/(-2*l0)
#   return(c(pval, effsize))
# }

# lrt.joint.lag.3 <- function(data, e=1, plot.CR=TRUE, resolution=15, verbose=FALSE, statistic=FALSE, bartlett=FALSE) {
#   
#   n <- max(data) + 1
#   seq <- aux.seq(e = e, w = 1, n = n - 1)
#   
#   # checking if the model is single-parameter only
#   if (verbose) {message("Finding MLE under LC model")}
#   m.one <- mle.ld(data = data, e = e)
#   logprob.one <- logprob(m = m.one[1], len = n, data = data, seq = aux.seq(e = e, w = 1, n = n - 1), k = 0, poisson = 0)
#   # m.minus1 <- suppressMessages(mle.ld(data = data, e = e, lag = -1, confint = FALSE))
#   # logprob.minus1 <- logprob_joint(m = m.minus1, e = e, param = 2^-1, data = data, len = n, option = "lag")
#   # m.minus0.999 <- suppressMessages(mle.ld(data = data, e = e, lag = -0.999, confint = FALSE))
#   # logprob.minus0.999 <- logprob_joint(m = m.minus0.999, e = e, param = 2^-0.999, data = data, len = n, option = "lag")
#   # m.01 <- suppressMessages(mle.ld(data = data, e = e, lag = 0.015, confint = FALSE))
#   # logprob.01 <- logprob_joint(m = m.01, e = e, param = sqrt(2^0.015-1), data = data, len = n, option = "lag")
#   # if (!is.finite(logprob.01)) {logprob.01 <- logprob_joint_boost(xm = m.01, xe = e, xparam = sqrt(2^0.015-1), data = data, len = n, option = "lag")}
#   
#   # print(c(logprob.minus1, logprob.minus0.999, logprob.one, logprob.01))
#   # stop()
#   
#   # if (logprob.one > logprob.01) {
#   # if (logprob.minus1 > logprob.minus0.999) {
#   # lag.max <- 0
#   # L.max <- 1
#   # m.max <- m.one[1]
#   # } else {
#   t1 <- Sys.time()
#   # finding joint MLEs using Newton-Raphson method
#   if (verbose) {message("Finding MLE under LC-lag model")}
#   mg1 <- m.one[1]
#   mg2 <- m_approx(data, e)
#   mg3 <- m_approx((min(data) + 1), e)
#   mg <- as.vector(na.exclude(c(as.vector(outer(c(mg1, mg2, mg3), c(1, 0.7, 0.5, 0.3, 0.1, 1.3, 1.5, 1.7, 2, 2.5, 3, 3.5, 4, 4.5, 5))), 0.1, 1)))
#   mg <- mg[mg < 1e6]
#   
#   # Lg <- c(0.5, 0.8, 1.0001, 1.01, seq(1.1, 5, length.out=10))
#   Lg <- c(1.0001, 1.01, seq(1.1, 5, length.out=10))
#   Lg <- sqrt(Lg-1)
#   
#   matrix1 <- c(-1, -1)
#   matrix0 <- c(0, 0)
#   
#   for (value.L in Lg) {
#     # print(value.L)
#     if (all((matrix1) > 0)) {break}
#     for (value.m in mg) {
#       # print(value.m)
#       if (all((matrix1) > 0)) {break}
#       matrix0[2] <- value.L
#       matrix0[1] <- value.m
#       iter <- 0
#       repeat {
#         iter <- iter + 1
#         if (iter > 30) {if (verbose) {message("joint lag: reached iteration limit")}; matrix1 <- c(-1, -1); break}
#         m0 <- matrix0[1]
#         L0 <- matrix0[2]
#         product <- score_fisher_product_joint(m = m0, e = e, param = L0^2+1, data = data, len = n, option = "lag")
#         # product <- score_fisher_product_joint(m = m0, e = e, param = L0, data = data, len = n, option = "lag")
#         product <- product[1:2]
#         # if (any(is.finite(product) == FALSE)) {if (verbose) {message("joint lag: boost")}; product <- score_fisher_product_joint_boost(m0, e, L0, data, n); product <- product[1:2]}
#         if (any(is.finite(product) == FALSE)) {if (verbose) {message("joint lag: non-numeric value encountered")}; matrix1 <- c(-1,-1); break}
#         product <- product[1:2]
#         matrix1 = matrix0 + product
#         # if (matrix1[1] <= 0 | matrix1[2] <= 1 | any(matrix1 > 1e6)) {if (verbose) {message("joint lag: value oob")}; matrix1 <- c(-1, -1); break}
#         if (any(matrix1[1] <= 0) | any(matrix1 > 1e6)) {if (verbose) {message("joint lag: value oob")}; matrix1 <- c(-1, -1); break}
#         # if (matrix1[1] <= 0 | matrix1[2] < 0.5 | any(matrix1 > 1e6)) {if (verbose) {message("joint lag: value oob")}; matrix1 <- c(-1, -1); break}
#         if (sqrt(sum((product)^2 / sum(matrix0^2))) < 1e-6) {break} else {matrix0 <- matrix1}
#       }
#     }
#   }
#   
#   # if (any(matrix1 < 0)) {
#   #   mg.ld <- seq(0.5*(min(data) + 1), max(data), length.out = 50)
#   #   if (e != 1) mg.ld <- c(mg.ld/e, mg.ld)
#   #   for (value.p in mg.p) {
#   #     if (all((matrix1) > 0)) {break}
#   #     for (value.ld in mg.ld) {
#   #       if (all((matrix1) > 0)) {break}
#   #       matrix0[2] <- value.p
#   #       matrix0[1] <- value.ld
#   #       iter <- 0
#   #       repeat {
#   #         iter <- iter + 1
#   #         if (iter > 40) {if (verbose) {message("joint residual: reached iteration limit")}; matrix1 <- c(-1, -1); break}
#   #         m0.ld <- matrix0[1]
#   #         m0.p <- matrix0[2]
#   #         product <- score_fisher_product_joint(m0.ld, e, m0.p, data, n, "mixed")
#   #         product <- product[1:2]
#   #         if (any(is.finite(product) == FALSE)) {message(if (verbose) "joint residual: boost"); product <- score_fisher_product_joint_boost(m0.ld, e, m0.p, data, n, "mixed")}
#   #         if (any(is.finite(product) == FALSE)) {if (verbose) {message("joint residual: non-numeric value encountered")}; matrix1 <- c(-1,-1); break}
#   #         product <- product[1:2]
#   #         matrix1 = matrix0 + product
#   #         if (any(matrix1 <= 0)) {if (verbose) {message("joint residual: value oob")}; matrix1 <- c(-1, -1); break}
#   #         if (sqrt(sum((product)^2 / sum(matrix0^2))) < 1e-4) {break} else {matrix0 <- matrix1}
#   #       }
#   #     }
#   #   }
#   # }
#   if (any(matrix1 < 0)) {
#     if (verbose) {message("Couldn't estimate joint estimates of m and residual m")}
#     if (plot.CR==FALSE) {
#       return(list("pval" = NA, "mutations" = m.one[1], "adj_mutations" = NA, "phenotypic_lag" = NA))
#     } else {
#       return(list("pval" = NA, "mutations" = m.one[1], "adj_mutations" = NA, "phenotypic_lag" = NA, "plot" = NA))
#     }
#   } else {
#     L.max <- matrix1[2]
#     lag.max <- log(L.max^2+1, 2)
#     # lag.max <- log(L.max, 2)
#     m.max <- matrix1[1]
#   }
#   # }
#   
#   # t2 <- Sys.time()
#   # # print(eigen(matrix(c(J11, J12, J12, J22), ncol=2))$value)
#   # # print(optim_m_param_joint(current_m = sqrt(m_guess(data = data, e = e, w = 1, lag = 0, death = 0, poisson = 1)[1]), lower_m = sqrt(m.one[1]), upper_m = sqrt(m.one[1]/10),
#   # # current_param = 1, lower_param = sqrt(1e-18), upper_param = sqrt(mean(data)), e = e, data = data, len = n, option = "mixed", verbose = verbose))
#   # ans1 <- (optim_m_param_joint(current_m = sqrt(m_guess(data = data, e = e, w = 1, lag = 1, death = 0, poisson = 0)[1]), lower_m = sqrt(m.one[1]), upper_m = sqrt(m_guess(data = data, e = e, w = 1, lag = 3, death = 0, poisson = 0)[1]),
#   #                              current_param = 1, lower_param = sqrt(1e-18), upper_param = sqrt(3), e = e, data = data, len = n, option = "lag", verbose = verbose))
#   # t3 <- Sys.time()
#   # ans2 <- (optim_m_param_joint_2(current_m = sqrt(m_guess(data = data, e = e, w = 1, lag = 0.1, death = 0, poisson = 0)[1]), current_param = sqrt(0.1), e = e, data = data, len = n, option = "lag", verbose = verbose))
#   # t4 <- Sys.time()
#   # print(c(m.max, lag.max))
#   # print(ans1)
#   # print(ans2)
#   # print(t2-t1)
#   # print(t3-t2)
#   # print(t4-t3)
#   # # stop()
#   # return(0)
#   
#   if (verbose) {message("Calculating p value")}
#   logprob.max <- logprob_joint(m = m.max, e = e, param = L.max^2+1, data = data, len = n, option = "lag")
#   
#   # pval <- 1 - pchisq(-2*(logprob.one - logprob.max),1)
#   if (m.max==m.one[1]) {
#     pval <- 0.5
#     if (statistic==TRUE) pval <- c(pval, 0)
#   } else {
#     lambda <- -2*(logprob.one - logprob.max)
#     if (bartlett==TRUE) {
#       X.95 <- which(cumsum(prob_joint(m = m.max, e = e, param = 1, len = 100000, option = "lag"))>0.999)[1]
#       if (is.na(X.95)) X.95 <- 100000
#       lambda.star <- lapply(1:200, function(i) {lrt.joint.lag(floor(simu.cultures.lag(length(data), m.max*1e-9, 1, 1, 1, 1e9, trim = X.95)*e), e, plot.CR = F, verbose = F, statistic = T, bartlett = F)[[1]][2]})
#       # lambda.star <- lapply(1:200, function(i) {lrt.joint.lag(floor(simu.cultures.lag(length(data), m.max*1e-9, 1, 1, 1, 1e9)*e), e, plot.CR = F, verbose = F, statistic = T, bartlett = F)[[1]][2]})
#       lambda.star <- unlist(lambda.star)
#       while (any(lambda.star==0)) {
#         lambda.star[which(lambda.star==0)] <- lrt.joint.lag(floor(simu.cultures.lag(length(data), m.max*1e-9, 1, 1, 1, 1e9, trim = X.95)*e), e, plot.CR = F, verbose = F, statistic = T, bartlett = F)[[1]][2]
#       }
#       lambda <- lambda/mean(unlist(lambda.star))
#     }
#     # lrt.stat <- -2*(logprob.one - logprob.max)
#     pval <- 1 - emdbook::pchibarsq(lambda,1)
#     if (statistic==TRUE) pval <- c(pval, lambda)
#   }
#   
#   if (plot.CR==FALSE) {
#     # if (m.max != m.one[1]) {
#     #   logprob.95 <- logprob.max-0.5*qchisq(0.95,2)
#     # } else {
#     #   logprob.95 <- logprob.max-0.5*emdbook::qchibarsq(0.95,2)
#     # }
#     # logprob.true <- logprob_joint(20, e, 2^0, data, n, "lag")
#     # if (logprob.true>=logprob.95) {confint=1} else {confint=0}
#     # return(list("pval" = pval, "m_one" = m.one[1], "m_lag" = m.max, "lag" = lag.max, "interval" = confint))
#     return(list("pval" = pval, "mutations" = m.one[1], "adj_mutations" = m.max, "phenotypic_lag" = lag.max))
#   } else {
#     
#     # plotting confidence region using grid-searching-type algorithm
#     if (resolution < 10) {resolution <- 10}
#     # lag.map <- seq(0.1*lag.max, ifelse(lag.max==0, 2, 1.9*lag.max), length.out=resolution)
#     lag.map <- seq(ifelse(lag.max<=0, -1, 0.1*lag.max), ifelse(lag.max<=0, 2, 1.9*lag.max), length.out=resolution)
#     m.map <- seq(0.5*m.max, 1.5*m.max, length.out=resolution)
#     
#     if (verbose) {message("Creating plot")}
#     logprob.95 <- logprob.max-0.5*qchisq(0.95,2)
#     table <- matrix(nrow=resolution, ncol=resolution)
#     rownames(table) <- m.map
#     colnames(table) <- lag.map
#     for (i in 1:resolution) {
#       for (j in 1:resolution) {
#         table[i,j] <- logprob_joint(m = m.map[i], e = e, param = 2^lag.map[j], data = data, len = n, option = "lag")
#       }
#     }
#     for (i in 1:3) {
#       while(!all(table[1,]<logprob.95)) {
#         if (m.map[1]==0) {break}
#         m.down <- (1-1/resolution)*m.map[1]
#         if (m.down<0.01) {m.down <- 0}
#         m.map <- c(m.down, m.map)
#         appendix <- mapply(function(x) logprob_joint(m = m.down, e = e, param = 2^lag.map[x], data = data, len = n, option = "lag"), x=1:length(lag.map))
#         table <- rbind(appendix, table)
#         rownames(table)[1] <- m.down
#       }
#       while(!all(table[length(rownames(table)),]<logprob.95)) {
#         m.up <- (1+1/resolution)*m.map[length(m.map)]
#         m.map <- c(m.map, m.up)
#         appendix <- mapply(function(x) logprob_joint(m = m.up, e = e, param = 2^lag.map[x], data = data, len = n, option = "lag"), x=1:length(lag.map))
#         table <- rbind(table, appendix)
#         rownames(table)[length(m.map)] <- m.up
#       }
#       while(!all(table[,1]<logprob.95)) {
#         if (lag.map[1]==0) {break}
#         lag.down <- (1-1/resolution)*lag.map[1]
#         if (lag.down<0.01) {lag.down <- 0}
#         # if (2^(lag.map[1])<=0.5) {break}
#         # lag.down <- log((1-1/resolution)*2^(lag.map[1]), 2)
#         # if (2^(lag.down)<=0.5) {lag.down <- -1}
#         lag.map <- c(lag.down, lag.map)
#         appendix <- mapply(function(x) logprob_joint(m = m.map[x], e = e, param = 2^lag.down, data = data, len = n, option = "lag"), x=1:length(m.map))
#         table <- cbind(appendix, table)
#         colnames(table)[1] <- lag.down
#       }
#       while(!all(table[,length(colnames(table))]<logprob.95)) {
#         lag.up <- (1+1/resolution)*lag.map[length(lag.map)]
#         lag.map <- c(lag.map, lag.up)
#         appendix <- mapply(function(x) logprob_joint(m = m.map[x], e = e, param = 2^lag.up, data = data, len = n, option = "lag"), x=1:length(m.map))
#         table <- cbind(table, appendix)
#         colnames(table)[length(lag.map)] <- lag.up
#       }
#     }
#     plot <- ggplot2::ggplot(reshape2::melt(table),ggplot2::aes(x=Var1, y=Var2)) + ggplot2::geom_contour(ggplot2::aes(z=value), breaks=logprob.95, colour="blue") + ggplot2::geom_point(ggplot2::aes(x=m.max, y=lag.max), colour="blue")  + ggplot2::geom_point(ggplot2::aes(x=m.one[1], y=0), colour="black") + ggplot2::geom_segment(ggplot2::aes(x = m.one[2], y = 0, xend = m.one[3], yend = 0), colour="black", size=0.25) + ggplot2::labs(x="m", y="lag", title="MLE and 95% confidence region") + ggplot2::theme_bw()
#     # plotly::plot_ly(z = as.data.frame(table)) %>% add_surface(z=table)
#     # plotly::plot_ly(type="surface", x=rownames(table), y=colnames(table), z=table)
#     return(list("pval" = pval, "mutations" = m.one[1], "adj_mutations" = m.max, "phenotypic_lag" = lag.max, "plot" = plot))
#   }
# }

mle.fold <- function(data, e=NULL, w=NULL, cv=NULL, Nt=NULL, fun="fold X1/X2") {
  
  if (!is.list(data) | !length(data) > 1) stop("Colony counts must be a list of length 2 to 6.")
  
  # checking if fold function is correct
  fun <- tryCatch(match.arg(fun, c("fold X1/X2", "substraction X1-X2", "double fold (X1/X2)/(X3/X4)", "background substraction fold (X1-X2)/(X3-X2)")), error = function(err) fun)
  if (fun == "fold X1/X2") {function.Y <- "X1/X2"}
  else if (fun == "substraction X1-X2") {function.Y <- "X1-X2"}
  else if (fun == "double fold (X1/X2)/(X3/X4)") {function.Y <- "(X1/X2)/(X3/X4)"}
  else if (fun == "background substraction fold (X1-X2)/(X3-X2)") {function.Y <- "(X1-X3)/(X2-X3)"}
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

calc.rate <- function(cs=NULL, eff=NULL, nt=NULL, fit=NULL, cv=NULL, vt=NULL, vs=NULL,
                      ds=NULL, vn=NULL, dn=NULL, cn=NULL, model="LD", pval=TRUE) {
  
  model <- match.arg(arg = model, choices = c("LD", "B"), several.ok = FALSE)
  if (model != "LD" & model != "B") {stop("Incorrect model. Model must be either 'LD' or 'B'.")}
  
  if (missing(cs)) {stop("Counts on selective medium were not provided.")}
  if (!is.vector(cs)) {stop("Counts on selective medium have invalid data structure.")}
  if (!is.list(cs)) {cs <- list(cs)}
  datasets <- length(cs)
  
  for (i in 1:datasets) {
    cs[[i]] <- suppressWarnings(as.double(as.vector(cs[[i]])))
    if (any(!is.finite(cs[[i]]))) {stop(paste("There are non-numeric characters in counts on selective medium in dataset ", names(cs)[i], ".", sep = ""))}
    else if (length(cs[[i]]) < 2) {stop(paste("There are less than 2 colony counts on selective medium in dataset ", names(cs)[i], ".", sep = ""))}
    else if (any(cs[[i]] < 0)) {stop(paste("There are negative values in colony counts on selective medium in dataset ", names(cs)[i], ".", sep = ""))}
    else if (sum(cs[[i]]) == 0) {stop(paste("At least one colony count on selective medium in datatset ", names(cs)[i], " must be bigger than 0.", sep = ""))}
  }
  
  if (((missing(eff) & missing(vs) & missing(ds) & missing(nt)) | (missing(eff) & missing(vt) & missing(vs) & missing(ds) & !missing(nt))) & !(missing(eff) & missing(vs) & missing(ds) & missing(nt) & missing(vn) & missing(dn) & missing(cn) & !missing(vt))) {
    message("Plating efficiency will be set to 1.")
    eff <- rep(1, datasets)
  } else if (missing(eff) & (missing(vt) | missing(vs) | missing(ds))) {
    if (missing(vt)) {stop("Cannot calculate plating efficiency because total culture volume is missing.")}
    else if (missing(vs)) {stop("Cannot calculate plating efficiency because volume plated on selective medium is missing.")}
    else if (missing(ds)) {stop("Cannot calculate plating efficiency because dilution factor on selective medium is missing.")}
  }
  
  if (((missing(nt) & missing(vn) & missing(dn) & missing(cn) & missing(eff)) | (missing(nt) & missing(vn) & missing(dn) & missing(cn) & missing(vt) & !missing(eff))) & !(missing(eff) & missing(vs) & missing(ds) & missing(nt) & missing(vn) & missing(dn) & missing(cn) & !missing(vt))) {
    message("Total culture size will be set to 1.")
    nt <- rep(1, datasets)
  } else if (missing(nt) & (missing(vt) | missing(vn) | missing(dn) | missing(cn))) {
    if (missing(vt)) {stop("Cannot calculate total culture size because total culture volume is missing.")}
    else if (missing(vn)) {stop("Cannot calculate total culture size because volume plated on non-selective medium is missing.")}
    else if (missing(dn)) {stop("Cannot calculate total culture size because dilution factor on non-selective medium is missing.")}
    else if (missing(cn)) {stop("Cannot calculate total culture size because counts on non-selective medium is missing.")}
  }
  
  if (missing(eff) | missing(nt)) {
    vt <- suppressWarnings(as.double(as.vector(vt)))
    if (any(!is.finite(vt))) {stop("There are non-numeric characters in total culture volume.")}
    else if (length(vt) < datasets) {stop(paste("There are ", length(vt), " values in total culture volume, expected ", datasets, ".", sep = ""))}
    else if (length(vt) > datasets) {warning(paste("There are ", length(vt), " values in total culture volume, only first ", datasets, " will be used.", sep = "")); vt <- vt[1:datasets]}
    if (any(vt <= 0)) {stop("Each total culture volume must be > 0.")}
  }
  
  if (missing(eff)) {
    vs <- suppressWarnings(as.double(as.vector(vs)))
    if (any(!is.finite(vs))) {stop("There are non-numeric characters in volume plated on selective medium.")}
    else if (length(vs) < datasets) {stop(paste("There are ", length(vs), " values in volume plated on selective medium, expected ", datasets, ".", sep = ""))}
    else if (length(vs) > datasets) {warning(paste("There are ", length(vs), " values in volume plated on selective medium, only first ", datasets, " will be used.", sep = "")); vs <- vs[1:datasets]}
    if (any(vs <= 0)) {stop("Each volume plated on selective medium must be > 0.")}
    
    ds <- suppressWarnings(as.double(as.vector(ds)))
    if (any(!is.finite(ds))) {stop("There are non-numeric characters in dilution factor on selective medium.")}
    else if (length(ds) < datasets) {stop(paste("There are ", length(ds), " values in dilution factor on selective medium, expected ", datasets, ".", sep = ""))}
    else if (length(ds) > datasets) {warning(paste("There are ", length(ds), " values in dilution factor on selective medium, only first ", datasets, " will be used.", sep = "")); ds <- ds[1:datasets]}
    if (any(ds <= 0)) {stop("Each dilution factor on selective medium must be > 0.")}
    
    if (any(vs / ds > vt)) {stop("The amount plated on selective medium is bigger than total culture volume.")}
    
    eff <- vs / ds / vt
    
  } else {
    eff <- suppressWarnings(as.double(as.vector(eff)))
    if (any(!is.finite(eff))) {stop("There are non-numeric characters in plating efficiency.")}
    else if (length(eff) < datasets) {stop(paste("There are ", length(eff), " values in plating efficiency, expected ", datasets, ".", sep = ""))}
    else if (length(eff) > datasets) {warning(paste("There are ", length(eff), " values in plating efficiency, only first ", datasets, " will be used.", sep = "")); eff <- eff[1:datasets]}
    if (any(eff <= 0 | eff > 1)) {stop("Each plating efficiency must be > 0 and \U2264 1.")}
  }
  
  if (missing(nt)) {
    vn <- suppressWarnings(as.double(as.vector(vn)))
    if (any(!is.finite(vn))) {stop("There are non-numeric characters in volume plated on non-selective medium.")}
    else if (length(vn) < datasets) {stop(paste("There are ", length(vn), " values in volume plated on non-selective medium, expected ", datasets, ".", sep = ""))}
    else if (length(vn) > datasets) {warning(paste("There are ", length(vn), " values in volume plated on non-selective medium, only first ", datasets, " will be used.", sep = "")); vn <- vn[1:datasets]}
    if (any(vn <= 0)) {stop("Each volume plated on non-selective medium must be > 0.")}
    
    dn <- suppressWarnings(as.double(as.vector(dn)))
    if (any(!is.finite(dn))) {stop("There are non-numeric characters in dilution factor on non-selective medium.")}
    else if (length(dn) < datasets) {stop(paste("There are ", length(dn), " values in dilution factor on non-selective medium, expected ", datasets, ".", sep = ""))}
    else if (length(dn) > datasets) {warning(paste("There are ", length(dn), " values in dilution factor on non-selective medium, only first ", datasets, " will be used.", sep = "")); dn <- dn[1:datasets]}
    if (any(dn <= 0)) {stop("Each dilution factor on non-selective medium must be > 0.")}
    
    if (any(vn / dn > vt)) {stop("The amount plated on non-selective medium is bigger than total culture volume.")}
    
    CVmissing <- missing(cv)
    if (!is.vector(cn)) {stop("Counts on non-selective medium have invalid data structure.")}
    if (!is.list(cn)) {cn <- list(cn)}  
    if (length(cn) < datasets) {stop(paste("There are ", length(cn), " sets in counts on non-selective medium, expected ", datasets, ".", sep = ""))}
    if (length(cn) > datasets) {warning(paste("There are ", length(cn), " sets in counts on non-selective medium, only first ", datasets, " will be used.", sep = "")); cn <- cn[[1:datasets]]}
    for (i in 1:datasets) {
      cn[[i]] <- suppressWarnings(as.double(as.vector(cn[[i]])))
      if (any(!is.finite(cn[[i]]))) {stop(paste("There are non-numeric characters in counts on non-selective medium in dataset ", names(cn)[i], ".", sep = ""))}
      if (any(cn[[i]] < 0)) {stop(paste("There are negative values in colony counts on non-selective medium in dataset ", names(cn)[i], ".", sep = ""))}
      if (sum(cn[[i]]) == 0) {stop(paste("At least one colony count on non-selective medium in datatset ", names(cn)[i], " must be bigger than 0.", sep = ""))}
      if (model == "B" & CVmissing == TRUE) {
        if (length(cn[[i]]) < 2) {stop(paste("There are less than 2 colony counts on non-selective medium in dataset ", names(cn)[i], ", at least 2 are required to calculate CV.", sep = ""))}
        cv[i] <- sd(cn[[i]]) / mean(cn[[i]])
      } else {
        if (length(cn[[i]]) < 1) {stop(paste("There is less than 1 colony count on non-selective medium in dataset ", names(cn)[i], ".", sep = ""))}
      }
      nt[i] <- mean(cn[[i]])*vn[[i]]/dn[[i]]/vt[[i]]
    }
    
  } else {
    nt <- suppressWarnings(as.double(as.vector(nt)))
    if (any(!is.finite(nt))) {stop("There are non-numeric characters in total culture size.")}
    else if (length(nt) < datasets) {stop(paste("There are ", length(nt), " values in total culture size, expected ", datasets, ".", sep = ""))}
    else if (length(nt) > datasets) {warning(paste("There are ", length(nt), " values in total culture size, only first ", datasets, " will be used.", sep = "")); nt <- nt[1:datasets]}
    if (any(nt <= 0)) {stop("Each total culture size must be > 0.")}
  }
  
  if (model == "B") {
    if (missing(cv)) {warning("Bartlett model was chosen but CV is missing, model will switch to Luria-Delbruck."); model <- "LD"}
    else {
      cv <- suppressWarnings(as.double(as.vector(cv)))
      if (any(!is.finite(cv))) {stop("There are non-numeric characters in coefficient of variation.")}
      else if (length(cv) < datasets) {stop(paste("There are ", length(cv), " values in coefficient of variation, expected ", datasets, ".", sep = ""))}
      else if (length(cv) > datasets) {warning(paste("There are ", length(cv), " values in coefficient of variation, only first ", datasets, " will be used.", sep = "")); cv <- cv[1:datasets]}
      if (any(cv <= 0)) {stop("Each coefficient of variation must be > 0.")}
      if (any(cv < 1e-5)) {warning("Coefficients of variation smaller than 1e-5 will be set to 1e-5."); cv[cv < 1e-5] <- 1e-5}
      if (any(cv > 10)) {warning("Coefficients of variation bigger than 10 will be set to 10."); cv[cv < 10] <- 10}
    }
  } else {
    if (!missing(cv)) {message("Coefficient of variation will be ignored because it is not a valid parameter under the Luria-Delbruck model.")}
  }
  
  if (model == "LD") {
    if (missing(fit)) {message("Mutant relative fitness will be set to 1."); fit <- rep(1, datasets)}
    else {
      fit <- suppressWarnings(as.double(as.vector(fit)))
      if (any(!is.finite(fit))) {stop("There are non-numeric characters in mutant relative fitness.")}
      else if (length(fit) < datasets) {stop(paste("There are ", length(fit), " values in mutant relative fitness, expected ", datasets, ".", sep = ""))}
      else if (length(fit) > datasets) {warning(paste("There are ", length(fit), " values in mutant relative fitness, only first ", datasets, " will be used.", sep = "")); fit <- fit[1:datasets]}
      if (any(fit <= 0)) {stop("Each mutant relative fitness must be > 0.")}
    }
  } else {
    if (!missing(fit)) {message("Mutant relative fitness will be ignored because it is not a valid parameter under the Bartlett model.")}
  }
  
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
                                   deathx = datax$Death, deathy = datay$Death, phix = datax$Inoculum/Ntx, phiy = datay$Inoculum/Nty,
                                   cvx = cv0x, cvy = cv0y, Nx = Ntx, Ny = Nty, Mx = Mx, My = My, verbose = F),
                     error=function(err) {c(NA)})
  
  return(result)
}