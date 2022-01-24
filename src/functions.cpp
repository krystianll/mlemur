// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
using namespace Rcpp;
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <chrono>
namespace mp = boost::multiprecision;
using mpfr_20 = mp::number<mp::mpfr_float_backend<20,mp::allocate_stack>,mp::et_off>;

/// EXPERIMENT SIMULATION ///

// [[Rcpp::export]]
Rcpp::List rluria_list(const int n=10, const double rate=1e-8, const double N0=1, const double Nt=1e9, const double mut_fit=1.0, const int type=0,
                       const double wt_dp=0, const double mut_dp=0, const double lag=0, const double e=1, const double cv=0, const int trim=0) {
  NumericVector res(n);
  
  double beta1s=1.0*(1.0-2.0*wt_dp)/(1.0-wt_dp);
  double beta2, beta2s, delta2;
  if (type==0){
    beta2=mut_fit;
    beta2s=beta2*(1.0-2.0*mut_dp)/(1.0-mut_dp);
  } else {
    beta2s=mut_fit*beta1s;
    beta2=mut_fit*(1.0-mut_dp)/(1.0-2.0*mut_dp);
  };
  delta2=beta2*mut_dp/(1.0-mut_dp);
  double tg=log(2)/beta2s;
  
  int i,j;
  // std::vector<double> Ttot(n), Tlag(n), m(n), mutations(n);
  NumericVector Ttot(n), Tlag(n), m(n), mutations(n);
  // std::vector<double> Nts(n);
  NumericVector Nts(n);
  if (cv==0) {
    Ttot[0]=log(Nt/N0)/beta1s;
    Tlag[0]=Ttot[0]-tg*lag;
    m[0]=rate*N0*(exp(beta1s*Tlag[0])-1.0);
    Nts[0]=Nt;
    for (i=1;i<n;++i) {
      Ttot[i]=Ttot[0];
      Tlag[i]=Tlag[0];
      m[i]=m[0];
      Nts[i]=Nt;
    };
  } else {
    // Nts=rlognorm(log(Nt)-log(1.0+cv*cv)/2.0, sqrt(log(1.0+cv*cv)), n);
    Nts=rlnorm(n, log(Nt)-log(1.0+cv*cv)/2.0, sqrt(log(1.0+cv*cv)));
    for (i=0;i<n;++i) {
      Ttot[i]=log(Nts[i]/N0)/beta1s;
      Tlag[i]=Ttot[i]-tg*lag;
      m[i]=rate*N0*(exp(beta1s*Tlag[i])-1.0);
    };
  }
  for (i=0;i<n;++i) {
    // mutations[i]=rpoisson(m[i]);
    mutations[i]=rpois(1, m[i])[0];
  }
  
  // std::vector<double> mt(n);
  NumericVector mt(n);
  double epoch, p, u, mutants;
  int R;
  
  for (i=0;i<=n-1;++i) {
    res[i]=0;
  };
  
  for (i=0;i<=n-1;++i) {
    if (mutations[i]==0) {
      res[i]+=0;
    } else {
      // mt=runiform(mutations[i], 0, 1);
      mt=runif(mutations[i], 0, 1);
      for (j=0; j<mutations[i]; j++) {
        mutants=0;
        
        epoch=log(1.0+mt[j]*(exp(beta1s*Tlag[i])-1.0))/(beta1s);
        p=exp(-beta2s*(Ttot[i]-epoch));
        u=1-delta2*(exp(beta2s*(Ttot[i]-epoch))-1.0)/(beta2*exp(beta2s*(Ttot[i]-epoch))-delta2);
        // R=rbinomial(u, 1);
        R=rbinom(1, 1, u)[0];
        
        if (R==1) {
          // mutants+=1+rgeometric(u*p);
          mutants+=1+rgeom(1, u*p)[0];
        };
        res[i]+=mutants;
      };
      if (e<1) {
        // res[i]=rbinomial(e, res[i]);
        res[i]=rbinom(1, res[i], e)[0];
      }
      if (trim>0) {
        if (res[i]>trim) {
          res[i]=trim;
        };
      };
    };
  };
  
  return Rcpp::List::create(_["mc"] = res, _["nt"] = Nts);
}

// [[Rcpp::export]]
NumericVector rluria_vec(const int n=10, const double rate=1e-8, const double N0=1, const double Nt=1e9, const double mut_fit=1.0, const int type=0,
                         const double wt_dp=0, const double mut_dp=0, const double lag=0, const double e=1, const double cv=0, const int trim=0) {
  NumericVector res(n);
  
  double beta1s=1.0*(1.0-2.0*wt_dp)/(1.0-wt_dp);
  double beta2, beta2s, delta2;
  if (type==0){
    beta2=mut_fit;
    beta2s=beta2*(1.0-2.0*mut_dp)/(1.0-mut_dp);
  } else {
    beta2s=mut_fit*beta1s;
    beta2=mut_fit*(1.0-mut_dp)/(1.0-2.0*mut_dp);
  };
  delta2=beta2*mut_dp/(1.0-mut_dp);
  double tg=log(2)/beta2s;
  
  int i,j;
  // std::vector<double> Ttot(n), Tlag(n), m(n), mutations(n);
  NumericVector Ttot(n), Tlag(n), m(n), mutations(n);
  if (cv==0) {
    Ttot[0]=log(Nt/N0)/beta1s;
    Tlag[0]=Ttot[0]-tg*lag;
    m[0]=rate*N0*(exp(beta1s*Tlag[0])-1.0);
    for (i=1;i<n;++i) {
      Ttot[i]=Ttot[0];
      Tlag[i]=Tlag[0];
      m[i]=m[0];
    };
  } else {
    // std::vector<double> Nts(n);
    NumericVector Nts(n);
    // Nts=rlognorm(log(Nt)-log(1.0+cv*cv)/2.0, sqrt(log(1.0+cv*cv)), n);
    Nts=rlnorm(n, log(Nt)-log(1.0+cv*cv)/2.0, sqrt(log(1.0+cv*cv)));
    for (i=0;i<n;++i) {
      Ttot[i]=log(Nts[i]/N0)/beta1s;
      Tlag[i]=Ttot[i]-tg*lag;
      m[i]=rate*N0*(exp(beta1s*Tlag[i])-1.0);
    };
  }
  for (i=0;i<n;++i) {
    // mutations[i]=rpoisson(m[i]);
    mutations[i]=rpois(1, m[i])[0];
  }
  
  // std::vector<double> mt(n);
  NumericVector mt(n);
  double epoch, p, u, mutants;
  int R;
  
  for (i=0;i<=n-1;++i) {
    res[i]=0;
  };
  
  for (i=0;i<=n-1;++i) {
    if (mutations[i]==0) {
      res[i]+=0;
    } else {
      // mt=runiform(mutations[i], 0, 1);
      mt=runif(mutations[i], 0, 1);
      for (j=0; j<mutations[i]; j++) {
        mutants=0;
        
        epoch=log(1.0+mt[j]*(exp(beta1s*Tlag[i])-1.0))/(beta1s);
        p=exp(-beta2s*(Ttot[i]-epoch));
        u=1-delta2*(exp(beta2s*(Ttot[i]-epoch))-1.0)/(beta2*exp(beta2s*(Ttot[i]-epoch))-delta2);
        // R=rbinomial(u, 1);
        R=rbinom(1, 1, u)[0];
        
        if (R==1) {
          // mutants+=1+rgeometric(u*p);
          mutants+=1+rgeom(1, u*p)[0];
        };
        res[i]+=mutants;
      };
      if (e<1) {
        // res[i]=rbinomial(e, res[i]);
        res[i]=rbinom(1, res[i], e)[0];
      }
      if (trim>0) {
        if (res[i]>trim) {
          res[i]=trim;
        };
      };
    };
  };
  
  return res;
}

/// AULIARY SEQUENCE FOR PROBABILITY COMPUTATION ///

/// GENERAL INTEGRAL FORMULA ///

// [[Rcpp::export]]
std::vector<double> aux_seq_integrate(double e=1, double w=1, double d=0, double lag=0, double phi=0, int n=10) {
  // std::vector<double> aux_seq_integrate(double e, double w, double d, double lag, double phi, int n) {
  using namespace boost::math::quadrature;
  tanh_sinh<double> tanh;
  double r=1./w;
  double chi=pow(2.,-lag);
  int k;
  std::vector<double> res(n+1);
  
  if (e==1 & d==0) {
    auto f = [&r, &k](double x) {
      return (r*pow(x,r)*pow(1.-x,k-1));
    };
    
    res[0]=-pow(chi,r);
    for (k=1;k<=n;++k) {
      res[k]=tanh.integrate(f, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
      // res[k]=tanh.integrate(f, pow(phi,w), chi)/(pow(chi,r)-phi);
    };
    
  } else {
    auto f1 = [&d, &r, &e](double x) {
      return ((r*pow(x,r-1)*((-1 + e)*x + d*(-e + x)))/(e*(-1 + x) + (-1 + d)*x));
    };
    auto f = [&d, &r, &k, &e](double x) {
      return (r*pow(1.-d,2)*pow(x,r)*e/pow(e+(1.-d-e)*x,2)*pow((1.-x)/(1.+(1.-d-e)/e*x),k-1));
    };
    
    res[0]=-pow(chi,r)+tanh.integrate(f1, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
    // res[0]=-pow(chi,r)+tanh.integrate(f1, pow(phi,w), chi)/(pow(chi,r)-phi);
    
    for (k=1;k<=n;++k) {
      res[k]=tanh.integrate(f, pow(phi,w), chi)/(1.-phi*pow(chi,-r));
      // res[k]=tanh.integrate(f, pow(phi,w), chi)/(pow(chi,r)-phi);
    };
  };
  
  return res;
}

/// SEQUENCE FOR LD DISTRIBUTION WITH OPTIONAL PARTIAL PLATING AND DIFFERENTIAL GROWTH ///

// not exported
std::vector<double> aux_seq_exact(double e, double w, int n) {
  std::vector<double> hSeq(n+1);
  int i;
  double r=1.0/w;
  double z = 1.0-e;
  hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    hSeq[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// not exported
std::vector<double> aux_seq_exact_boost(double ee, double w, int n) {
  std::vector<double> hSeq(n+1);
  std::vector<mpfr_20> hSeq2(n+1);
  int i;
  double r=1.0/w;
  double z = 1.0-ee;
  mpfr_20 e=ee;
  hSeq2[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    hSeq2[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  for (i=0; i<=n; ++i) {
    hSeq[i]=hSeq2[i].convert_to<double>();
  };
  return(hSeq);
}

// not exported
std::vector<double> aux_seq_boost(double e, double w, int n) {
  std::vector<double> hSeq(n+1);
  int i;
  double r=1.0/w;
  double z = 1.0-e;
  hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    hSeq[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_boost(mpfr_20 e, mpfr_20 w, int n) {
  std::vector<mpfr_20> hSeq(n+1);
  int i;
  mpfr_20 r=1.0/w;
  mpfr_20 z = 1.0-e;
  hSeq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(1.0)}, {2.0 + r}, z)/(1.0+r);
  for (i=1; i<=n; ++i) {
    hSeq[i]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+i+1.0}, z)*boost::math::beta(i,r+1.0)*pow(e,r)*r;
  };
  return(hSeq);
}

// [[Rcpp::export]]
std::vector<double> aux_seq(double e, double w, int n, int option=0) {
  // NumericVector aux_seq(double e, double w, int n, int option) {
  if (option==1) {return aux_seq_exact(e, w, n);}
  else if (option==2) {return aux_seq_exact_boost(e, w, n);};
  
  if (n<1 || e<=0 || e>1 || w<=0) {return std::vector<double> {0};};
  
  std::vector<double> seq(n+1);
  int i;
  double r=1.0/w;
  double z=1.0-e;
  
  double lowerlimit=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z)*boost::math::beta(1,r+1.0)*pow(e,r)*r;
  double upperlimit=boost::math::hypergeometric_pFq({r, r+1.0}, {r+1.0+n}, z)*boost::math::beta(n,r+1.0)*pow(e,r)*r;
  if (lowerlimit==0 || upperlimit==0) {return std::vector<double> {0};};
  
  if (e==1) {
    seq[0]=-1.0;
    seq[1]=1.0/(1.0+r)*r;
    if (n>1) {
      for (i=2; i<=n; ++i) {
        seq[i]=(i-1.0)/(i+r)*seq[i-1];
      };
    };
  } else if (w!=1 && e!=1) {
    std::vector<double> bSeq(n+1);
    bSeq[0]=-1.0;
    bSeq[1]=1.0/(1.0+r);
    if (n>1) {
      for (i=2; i<=n; ++i) {
        bSeq[i]=(i-1.0)/(i+r)*bSeq[i-1];
      };
    };
    std::vector<double> fSeq(n+1);
    fSeq[0]=0;
    if (n==2) {
      fSeq[1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z);
    } else if (n==3) {
      fSeq[1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z);
      fSeq[2]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+3.0}, z);
    } else if (n>=4) {
      if (e>=0.5) {
        fSeq[n]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+n+1.0}, z);
        fSeq[n-1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+n}, z);
        for (i=(n-1); i>=2; --i) {
          fSeq[i-1]=i*(i+1.0)*z*fSeq[i+1]/(r+i+1.0)/(r+i)/e+(r+i-2.0*i*z)*fSeq[i]/e/(r+i);
        };
      } else {
        fSeq[1]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+2.0}, z);
        fSeq[2]=boost::math::hypergeometric_pFq({r, r+1.0}, {r+3.0}, z);
        for (i=2; i<=(n-1); ++i) {
          fSeq[i+1]=(r+1.0+i)/(i*(i+1.0)*z)*(e*(r+i)*fSeq[i-1]-(r+i-2.0*i*z)*fSeq[i]);
        };
      };
    };
    seq[0]=-1.0+r*z*boost::math::hypergeometric_pFq({1.0, 1.0}, {2.0 + r}, z)/(1.0+r);
    for (i=1; i<=n; ++i) {
      seq[i]=pow(e,r)*bSeq[i]*fSeq[i]*r;
    };
  } else if (w==1 && e>=0.5) {
    seq[0]=log(e)*e/z;
    seq[n]=boost::math::hypergeometric_pFq({1.0, 2.0}, {n+2.0}, z)*e/n/(n+1.0);
    if (n>1) {
      for (i=n-1; i>=1; --i) {
        seq[i]=1.0/i/(i+1.0)-seq[i+1]*z/e;
      };
    };
  } else if (w==1 && e<0.5) {
    seq[0]=log(e)*e/z;
    seq[1]=(e/z)*(-1.0-log(e)/z);
    if (n>1) {
      for (i=2; i<=n; ++i) {
        seq[i]=e/z*(1.0/i/(i-1.0)-seq[i-1]);
      };
    };
  } else {
    return(std::vector<double> {0});
  };
  if (abs(seq[1]-lowerlimit)/lowerlimit>1e-9 || abs(seq[n]-upperlimit)/upperlimit>1e-9) {
    Rcout << "exact";
    return(aux_seq_exact(e, w, n));
  } else {
    return(seq);
  }
}

// [[Rcpp::export]]
Rcpp::List aux_seq_deriv1_deriv2(double e, double w, int n, double h=1.0e-4, bool boost=false) {
  // Rcpp::List aux_seq_deriv1_deriv2(double e, double w, int n, double h, bool boost) {
  std::vector<double> deriv1(n+1), deriv2(n+1);
  int i;
  
  if (boost) {
    std::vector<mpfr_20> seq(n+1), seqplus(n+1), seqminus(n+1), deriv11(n+1), deriv22(n+1);
    mpfr_20 xe=e;
    mpfr_20 xw=w;
    mpfr_20 whp=w*(1.0+h);
    mpfr_20 whm=w*(1.0-h);
    mpfr_20 dwp=whp-w;
    mpfr_20 dwm=w-whm;
    mpfr_20 dw=whp-whm;
    
    seq=aux_seq_boost(xe, xw, n);
    seqplus=aux_seq_boost(xe, whp, n);
    seqminus=aux_seq_boost(xe, whm, n);
    
    for (i=0;i<=n;++i) {
      deriv11[i]=(seqplus[i]-seqminus[i])/dw;
      deriv22[i]=(seqplus[i]+seqminus[i]-2*seq[i])/dwp/dwm;
      deriv1[i]=deriv11[i].convert_to<double>();
      deriv2[i]=deriv22[i].convert_to<double>();
    };
    
  } else {
    std::vector<double> seq(n+1), seqplus(n+1), seqminus(n+1);
    volatile double whp=w*(1.0+h);
    volatile double whm=w*(1.0-h);
    double dwp=whp-w;
    double dwm=w-whm;
    double dw=whp-whm;
    
    seq=aux_seq_boost(e, w, n);
    seqplus=aux_seq_boost(e, whp, n);
    seqminus=aux_seq_boost(e, whm, n);
    
    for (i=0;i<=n;++i) {
      deriv1[i]=(seqplus[i]-seqminus[i])/dw;
      deriv2[i]=(seqplus[i]+seqminus[i]-2*seq[i])/dwp/dwm;
    };
  };
  
  return Rcpp::List::create(Rcpp::Named("seq1") = deriv1,
                            Rcpp::Named("seq1") = deriv2);
}

/// SEQUENCES FOR LD WITH PHENOTYPIC LAG ///

// not exported
std::vector<double> aux_seq_lag(const double L, const double e, const int n) {
  
  std::vector<double> seq(n+1);
  int i;
  
  if (e == 1.0) {
    seq[0]=-1.0/L;
    if (n>=1) {
      for (i=1; i<=n; ++i) {
        seq[i]=(1.0-pow(1.0-1/L,i)*(1.0+i/L))/i/(i+1.0);
      };
    };
  } else {
    double odds=e/(1.0-e);
    double xi=e*(1.0-L);
    double div=xi/(xi-1.0);
    const int len=n+100; //to increase accuracy of backwards computation for small n
    std::vector<double> q(len+2);
    q[0]=odds*log(-e*L/(xi-1.0));
    if (n>=1) {
      q[1]=odds/(xi-1.0)-q[0];
      if (n>=2) {
        for (i=2; i<=len+1; ++i) {
          q[i]=odds/i/(i-1)*(1.0-(xi-i)*pow(div,i-1)/(xi-1.0));
        };
      };
    };
    seq[0]=q[0];
    if (n>=1) {
      seq[1]=-odds*q[0]+q[1];
    };
    if (n>=2) {
      if (e<0.5) {
        for (i=1; i<=n-1; ++i) {
          seq[i+1]=-odds*seq[i]+q[i+1];
        };
      } else {
        std::vector<double> seq2(len+1);
        seq2[len]=(1.0-e)*q[len+1];
        for (i=len-1; i>=2; --i) {
          seq2[i]=1/odds*(q[i+1]-seq2[i+1]);
        };
        for (i=2; i<=n; ++i) {
          seq[i]=seq2[i];
        };
      };
    };
  };
  
  return(seq);
}

// not exported
std::vector<mpfr_20> aux_seq_lag(const mpfr_20 L, const mpfr_20 e, const int n) {
  
  std::vector<mpfr_20> seq(n+1);
  int i;
  
  if (e == 1.0) {
    seq[0]=-1.0/L;
    if (n>=1) {
      for (i=1; i<=n; ++i) {
        seq[i]=(1.0-pow(1.0-1/L,i)*(1.0+i/L))/i/(i+1.0);
      };
    };
  } else {
    mpfr_20 odds=e/(1.0-e);
    mpfr_20 xi=e*(1.0-L);
    mpfr_20 div=xi/(xi-1.0);
    const int len=n+100; //to increase accuracy of backwards computation for small n
    std::vector<mpfr_20> q(len+2);
    q[0]=odds*log(-e*L/(xi-1.0));
    if (n>=1) {
      q[1]=odds/(xi-1.0)-q[0];
      if (n>=2) {
        for (i=2; i<=len+1; ++i) {
          q[i]=odds/i/(i-1)*(1.0-(xi-i)*pow(div,i-1)/(xi-1.0));
        };
      };
    };
    seq[0]=q[0];
    if (n>=1) {
      seq[1]=-odds*q[0]+q[1];
    };
    if (n>=2) {
      if (e<0.5) {
        for (i=1; i<=n-1; ++i) {
          seq[i+1]=-odds*seq[i]+q[i+1];
        };
      } else {
        std::vector<mpfr_20> seq2(len+1);
        seq2[len]=(1.0-e)*q[len+1];
        for (i=len-1; i>=2; --i) {
          seq2[i]=1/odds*(q[i+1]-seq2[i+1]);
        };
        for (i=2; i<=n; ++i) {
          seq[i]=seq2[i];
        };
      };
    };
  };
  
  return(seq);
}

// not exported
std::vector<double> aux_seq_lag_deriv1(const double L, const double e, const int n) {
  
  std::vector<double> seq(n+1);
  int i;
  
  if (e == 1.0) {
    seq[0]=1.0/L/L;
    double div=(L-1.0)/L;
    if (n>=1) {
      for (i=1; i<=n; ++i) {
        seq[i]=-pow(div,i-1)/pow(L,3);
      };
    };
  } else {
    double odds=e/(1.0-e);
    double xi=e*(1.0-L);
    double div=xi/(xi-1.0);
    const int len=n+100; //to increase accuracy of backwards computation for small n
    std::vector<double> q(len+2);
    q[0]=-e/L/(xi-1.0);
    if (n>=1) {
      q[1]=odds/L*(2*e-1.0-e*xi)/pow(xi-1,2);
      if (n>=2) {
        for (i=2; i<=len+1; ++i) {
          q[i]=odds*e*pow(div,i-2)/pow(xi-1,3);
        };
      };
    };
    seq[0]=q[0];
    if (n>=1) {
      seq[1]=-odds*q[0]+q[1];
    };
    if (n>=2) {
      if (e<0.5) {
        for (i=1; i<=n-1; ++i) {
          seq[i+1]=-odds*seq[i]+q[i+1];
        };
      } else {
        std::vector<double> seq2(len+1);
        seq2[len]=(1.0-e)*q[len+1];
        for (i=len-1; i>=2; --i) {
          seq2[i]=1/odds*(q[i+1]-seq2[i+1]);
        };
        for (i=2; i<=n; ++i) {
          seq[i]=seq2[i];
        };
      };
    };
  };
  
  return(seq);
}

// not exported
std::vector<mpfr_20> aux_seq_lag_deriv1(const mpfr_20 L, const mpfr_20 e, const int n) {
  
  std::vector<mpfr_20> seq(n+1);
  int i;
  
  if (e == 1.0) {
    seq[0]=1.0/L/L;
    mpfr_20 div=(L-1.0)/L;
    if (n>=1) {
      for (i=1; i<=n; ++i) {
        seq[i]=-pow(div,i-1)/pow(L,3);
      };
    };
  } else {
    mpfr_20 odds=e/(1.0-e);
    mpfr_20 xi=e*(1.0-L);
    mpfr_20 div=xi/(xi-1.0);
    const int len=n+100; //to increase accuracy of backwards computation for small n
    std::vector<mpfr_20> q(len+2);
    q[0]=-e/L/(xi-1.0);
    if (n>=1) {
      q[1]=odds/L*(2*e-1.0-e*xi)/pow(xi-1,2);
      if (n>=2) {
        for (i=2; i<=len+1; ++i) {
          q[i]=odds*e*pow(div,i-2)/pow(xi-1,3);
        };
      };
    };
    seq[0]=q[0];
    if (n>=1) {
      seq[1]=-odds*q[0]+q[1];
    };
    if (n>=2) {
      if (e<0.5) {
        for (i=1; i<=n-1; ++i) {
          seq[i+1]=-odds*seq[i]+q[i+1];
        };
      } else {
        std::vector<mpfr_20> seq2(len+1);
        seq2[len]=(1.0-e)*q[len+1];
        for (i=len-1; i>=2; --i) {
          seq2[i]=1/odds*(q[i+1]-seq2[i+1]);
        };
        for (i=2; i<=n; ++i) {
          seq[i]=seq2[i];
        };
      };
    };
  };
  
  return(seq);
}

// not exported
std::vector<double> aux_seq_lag_deriv2(const double L, const double e, const int n) {
  
  std::vector<double> seq(n+1);
  int i;
  
  if (e == 1.0) {
    seq[0]=-2.0/L/L/L;
    double div=(L-1.0)/L;
    if (n>=1) {
      for (i=1; i<=n; ++i) {
        seq[i]=(3*L-i-2.0)*pow(div,i-2)/pow(L,5);
      };
    };
  } else {
    double odds=e/(1.0-e);
    double xi=e*(1.0-L);
    double div=xi/(xi-1.0);
    const int len=n+100; //to increase accuracy of backwards computation for small n
    std::vector<double> q(len+2);
    q[0]=-e*((2*L-1.0)*e+1.0)/L/L/pow(xi-1.0,2);
    if (n>=1) {
      q[1]=-odds*(-1.0+3.0*xi+3.0*e*e*(2.0*L-1.0)+e*e*e*(1.0-3.0*L+2.0*L*L))/L/L/pow(1.0-xi,3);
      if (n>=2) {
        for (i=2; i<=len+1; ++i) {
          q[i]=odds*pow(div,i)*(2.0-3*xi-i)/e/pow(1.0-xi,2)/pow(L-1.0,3);
        };
      };
    };
    seq[0]=q[0];
    if (n>=1) {
      seq[1]=-odds*q[0]+q[1];
    };
    if (n>=2) {
      if (e<0.5) {
        for (i=1; i<=n-1; ++i) {
          seq[i+1]=-odds*seq[i]+q[i+1];
        };
      } else {
        std::vector<double> seq2(len+1);
        seq2[len]=(1.0-e)*q[len+1];
        for (i=len-1; i>=2; --i) {
          seq2[i]=1/odds*(q[i+1]-seq2[i+1]);
        };
        for (i=2; i<=n; ++i) {
          seq[i]=seq2[i];
        };
      };
    };
  };
  
  return(seq);
}

// not exported
std::vector<mpfr_20> aux_seq_lag_deriv2(const mpfr_20 L, const mpfr_20 e, const int n) {
  
  std::vector<mpfr_20> seq(n+1);
  int i;
  
  if (e == 1.0) {
    seq[0]=-2.0/L/L/L;
    mpfr_20 div=(L-1.0)/L;
    if (n>=1) {
      for (i=1; i<=n; ++i) {
        seq[i]=(3*L-i-2.0)*pow(div,i)/pow(L,3)/pow(L-1.0,2);
      };
    };
  } else {
    mpfr_20 odds=e/(1.0-e);
    mpfr_20 xi=e*(1.0-L);
    mpfr_20 div=xi/(xi-1.0);
    const int len=n+100; //to increase accuracy of backwards computation for small n
    std::vector<mpfr_20> q(len+2);
    q[0]=-e*((2*L-1.0)*e+1.0)/L/L/pow(xi-1.0,2);
    if (n>=1) {
      q[1]=-odds*(-1.0+3.0*xi+3.0*e*e*(2.0*L-1.0)+e*e*e*(1.0-3.0*L+2.0*L*L))/L/L/pow(1.0-xi,3);
      if (n>=2) {
        for (i=2; i<=len+1; ++i) {
          q[i]=odds*pow(div,i)*(2.0-3*xi-i)/e/pow(1.0-xi,2)/pow(L-1.0,3);
        };
      };
    };
    seq[0]=q[0];
    if (n>=1) {
      seq[1]=-odds*q[0]+q[1];
    };
    if (n>=2) {
      if (e<0.5) {
        for (i=1; i<=n-1; ++i) {
          seq[i+1]=-odds*seq[i]+q[i+1];
        };
      } else {
        std::vector<mpfr_20> seq2(len+1);
        seq2[len]=(1.0-e)*q[len+1];
        for (i=len-1; i>=2; --i) {
          seq2[i]=1/odds*(q[i+1]-seq2[i+1]);
        };
        for (i=2; i<=n; ++i) {
          seq[i]=seq2[i];
        };
      };
    };
  };
  
  return(seq);
}

// [[Rcpp::export]]
std::vector<double> aux_seq_lag_ext(double L, double e, int n, bool boost=false, int deriv=0) {
  // NumericVector aux_seq_lag_ext(double L, double e, int n, bool boost, int deriv) {
  
  std::vector<double> seq(n+1);
  int i;
  
  if (boost==false) {
    std::vector<double> seq2(n+1);
    
    if (deriv==0) {
      seq2=aux_seq_lag(L, e, n);
    } else if (deriv==1) {
      seq2=aux_seq_lag_deriv1(L, e, n);
    } else {
      seq2=aux_seq_lag_deriv2(L, e, n);
    };
    
    for (i=0; i<=n; ++i) {
      seq[i]=seq2[i];
    };
  } else {
    mpfr_20 xL=L;
    mpfr_20 xe=e;
    std::vector<mpfr_20> seq2(n+1);
    
    if (deriv==0) {
      seq2=aux_seq_lag(xL, xe, n);
    } else if (deriv==1) {
      seq2=aux_seq_lag_deriv1(xL, xe, n);
    } else {
      seq2=aux_seq_lag_deriv2(xL, xe, n);
    };
    
    for (i=0; i<=n; ++i) {
      seq[i]=seq2[i].convert_to<double>();
    };
  };
  
  return(seq);
}

/// SEQUENCE FOR LD WITH CELL DEATH ///

// not exported
std::vector<double> aux_seq_death(double e, double w, double d, int n) {
  std::vector<double> hSeq(n+1);
  int i;
  double r=1.0/w;
  if (e==1) {
    hSeq[0]=-1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0);
    for (i=1; i<=n; ++i) {
      hSeq[i]=pow(1.0-d,2)*r*boost::math::hypergeometric_pFq({i+1.0, r+1.0}, {r+i+1.0}, d)*boost::math::beta(i,r+1.0);
    };
  } else {
    if (d>=(1.0-2*e)) {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, d)+(1.0-e-d)/e*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, ((e+d-1.0)/e)));
      for (i=1; i<=n; ++i) {
        // Rcout << i << " ";
        try {
          hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*boost::math::hypergeometric_pFq({i+1.0, r+1.0}, {i+r+1.0}, ((e+d-1.0)/e));
        } catch (...) {
          Rcpp::stop("Error");
          // hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*static_cast<double>(boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(i+r+1.0)}, static_cast<mpfr_20>((e+d-1.0)/e)));
        }
        // Rcout << hSeq[i] << " | ";
      }
    } else {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({1.0, r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({1.0, r+1.0}, {r+2.0}, d)+(1.0-e-d)/e*(e/(1-d))*boost::math::hypergeometric_pFq({1.0, 1.0}, {r+2.0}, ((e+d-1.0)/(d-1.0))));
      for (i=1; i<=n; ++i) {
        try {
          hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*boost::math::hypergeometric_pFq({r+1.0, r}, {i+r+1.0}, ((e+d-1.0)/(d-1.0)));
        } catch (...) {
          Rcpp::stop("Error");
          // hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*static_cast<double>(boost::math::hypergeometric_pFq({static_cast<mpfr_20>(r+1.0), static_cast<mpfr_20>(r)}, {static_cast<mpfr_20>(i+r+1.0)}, static_cast<mpfr_20>((e+d-1.0)/(d-1.0))));
        }
      }
    }
  }
  
  return(hSeq);
}

// not exported
std::vector<mpfr_20> aux_seq_death(mpfr_20 e, mpfr_20 w, mpfr_20 d, int n) {
  std::vector<mpfr_20> hSeq(n+1);
  int i;
  mpfr_20 r=1.0/w;
  if (e==1) {
    hSeq[0]=-1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {r+static_cast<mpfr_20>(2.0)}, d)*boost::math::beta(r,2.0);
    for (i=1; i<=n; ++i) {
      hSeq[i]=pow(1.0-d,2)*r*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+i+1.0)}, d)*boost::math::beta(i,r+1.0);
    };
  } else {
    if (d>=(1.0-2*e)) {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {r+2.0}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, d)+(1.0-e-d)/e*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, ((e+d-1.0)/e)));
      for (i=1; i<=n; ++i) {
        hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(i+1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(i+r+1.0)}, ((e+d-1.0)/e));
      }
    } else {
      hSeq[0] = -1.0+r*d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), r}, {static_cast<mpfr_20>(r+2.0)}, d)*boost::math::beta(r,2.0)+r*(1.0-d)/(r+1.0)*(d*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(r+1.0)}, {static_cast<mpfr_20>(r+2.0)}, d)+(1.0-e-d)/e*(e/(1-d))*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(1.0), static_cast<mpfr_20>(1.0)}, {static_cast<mpfr_20>(r+2.0)}, ((e+d-1.0)/(d-1.0))));
      for (i=1; i<=n; ++i) {
        hSeq[i] = r*pow((1-d),2)/e*boost::math::beta(i,r+1.0)*pow(e/(1.0-d),r+1)*boost::math::hypergeometric_pFq({static_cast<mpfr_20>(r+1.0), r}, {static_cast<mpfr_20>(i+r+1.0)}, ((e+d-1.0)/(d-1.0)));
      }
    }
  }
  
  return(hSeq);
}

// [[Rcpp::export]]
std::vector<double> aux_seq_death_ext(double e, double w, double d, int n, bool boost=false) {
  // NumericVector aux_seq_death_ext(double e, double w, double d, int n, bool boost) {
  std::vector<double> hSeq(n+1);
  int i;
  if (boost==false) {
    std::vector<double> seq(n+1);
    seq=aux_seq_death(e, w, d, n);
    for (i=0;i<=n;++i) {
      hSeq[i]=seq[i];
    };
    return hSeq;
  } else {
    std::vector<mpfr_20> seq(n+1);
    seq=aux_seq_death(static_cast<mpfr_20>(e), static_cast<mpfr_20>(w), static_cast<mpfr_20>(d), n);
    for (i=0;i<=n;++i) {
      hSeq[i]=seq[i].convert_to<double>();
    };
    return hSeq;
  }
}

/// CALCULATION OF PMFS ///

/// PMF FOR LURIA-DELBRUCK APPROXIMATE DISTRIBUTION ///

// [[Rcpp::export]]
std::vector<double> prob_ld(const double m, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(m*seq[0]); 
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n]+=j*seq[j]*prob[n-j];
    };
    prob[n]*=m/n;
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> prob_ld(const mpfr_20 m, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(m*seq[0]); 
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n]+=j*seq[j]*prob[n-j];
    };
    prob[n]*=m/n;
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> prob_ld(const mpfr_20 m, std::vector<mpfr_20> &seq, const int len){
  int n,j;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(m*seq[0]); 
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n]+=j*seq[j]*prob[n-j];
    };
    prob[n]*=m/n;
  };
  
  return prob;
}

// not exported
std::vector<double> prob_ld_deriv(std::vector<double> &seq, std::vector<double> &prob, const int len){
  int n,j;
  std::vector<double> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

// not exported
std::vector<mpfr_20> prob_ld_deriv(std::vector<mpfr_20> &seq, std::vector<double> &prob, const int len){
  int n,j;
  std::vector<mpfr_20> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

// not exported
std::vector<mpfr_20> prob_ld_deriv(std::vector<double> &seq, std::vector<mpfr_20> &prob, const int len){
  int n,j;
  std::vector<mpfr_20> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

// not exported
std::vector<mpfr_20> prob_ld_deriv(std::vector<mpfr_20> &seq, std::vector<mpfr_20> &prob, const int len){
  int n,j;
  std::vector<mpfr_20> prob_deriv(len);
  
  for(n=1;n<=len-1;++n) {
    prob_deriv[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob_deriv[n-1]+=seq[j-1]*prob[n-j];
    };
  };
  
  return prob_deriv;
}

/// PMF FOR LURIA-DELBRUCK & GAMMA MIXTURE DISTRIBUTION ///

// [[Rcpp::export]]
std::vector<double> xi_seq(const double A, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<double> xi(len);
  
  for(n=1;n<=len-1;++n) {
    xi[n]=0;
  };
  
  xi[0]=log(1.0-A*seq[0]);
  xi[1]=-A*seq[1]/(1.0-A*seq[0]);
  for(n=2;n<=len-1;++n) {
    for(j=1;j<=n-1;++j) {
      xi[n]+=j*xi[j]*seq[n-j];
    };
    xi[n]=-A*(n*seq[n]-xi[n])/(n*(1.0-A*seq[0]));
  };
  
  return xi;
}

// not exported
std::vector<mpfr_20> xi_seq(const mpfr_20 A, std::vector<double> &seq, const int len){
  int n,j;
  std::vector<mpfr_20> xi(len);
  
  for(n=1;n<=len-1;++n) {
    xi[n]=0;
  };
  
  xi[0]=log(1.0-A*seq[0]);
  xi[1]=-A*seq[1]/(1.0-A*seq[0]);
  for(n=2;n<=len-1;++n) {
    for(j=1;j<=n-1;++j) {
      xi[n]+=j*xi[j]*seq[n-j];
    };
    xi[n]=-A*(n*seq[n]-xi[n])/(n*(1.0-A*seq[0]));
  };
  
  return xi;
}

// [[Rcpp::export]]
std::vector<double> prob_b0(const double A, const double k, const double seq0, std::vector<double> &xi, const int len){
  int n,j;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=1.0/pow((1.0-A*seq0),k);
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n] += j*xi[j]*prob[n-j];
    };
    prob[n]*=(-k/n);
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> prob_b0(const mpfr_20 A, const mpfr_20 k, const mpfr_20 seq0, std::vector<mpfr_20> &xi, const int len){
  int n,j;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=1.0/pow((1.0-A*seq0),k);
  for(n=1;n<=len-1;++n) {
    for(j=1;j<=n;++j){
      prob[n] += j*xi[j]*prob[n-j];
    };
    prob[n]*=(-k/n);
  };
  
  return prob;
}

// not exported
std::vector<double> prob_b0_deriv1(const double A, const double k, std::vector<double> &seq, std::vector<double> &xi, const int len){
  int n,j;
  std::vector<double> probk1(len),prob1(len);
  double seq0=seq[0];
  
  probk1=prob_b0(A, k+1.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    prob1[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob1[n-1]+=seq[j-1]*probk1[n-j];
    };
  };
  
  return prob1;
}

// not exported
std::vector<mpfr_20> prob_b0_deriv1(const mpfr_20 A, const mpfr_20 k, std::vector<double> &seq, std::vector<mpfr_20> &xi, const int len){
  int n,j;
  std::vector<mpfr_20> probk1(len),prob1(len);
  mpfr_20 seq0=seq[0];
  
  probk1=prob_b0(A, k+1.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    prob1[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob1[n-1]+=seq[j-1]*probk1[n-j];
    };
  };
  
  return prob1;
}

// not exported
std::vector<double> prob_b0_deriv2(const double A, const double k, std::vector<double> &seq, std::vector<double> &xi, const int len){
  int n,j;
  std::vector<double> probk2(len),conv(len),prob2(len);
  double seq0=seq[0];
  
  probk2=prob_b0(A, k+2.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    conv[n]=0;
    prob2[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      conv[n-1]+=seq[j-1]*seq[n-j];
    };
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob2[n-1]+=conv[j-1]*probk2[n-j];
    };
    prob2[n-1]*=((k+1)/k);
  };
  
  return prob2;
}

// not exported
std::vector<mpfr_20> prob_b0_deriv2(const mpfr_20 A, const mpfr_20 k, std::vector<double> &seq, std::vector<mpfr_20> &xi, const int len){
  int n,j;
  std::vector<mpfr_20> probk2(len),conv(len),prob2(len);
  mpfr_20 seq0=seq[0];
  
  probk2=prob_b0(A, k+2.0, seq0, xi, len);
  
  for(n=1;n<=len-1;++n) {
    conv[n]=0;
    prob2[n]=0;
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      conv[n-1]+=seq[j-1]*seq[n-j];
    };
  };
  
  for(n=1;n<=len;++n){
    for(j=1;j<=n;++j){
      prob2[n-1]+=conv[j-1]*probk2[n-j];
    };
    prob2[n-1]*=(k+1)/k;
  };
  
  return prob2;
}

/// PROBABILITY FOR POISSON DISTRIBUTION ///
// not exported
std::vector<double> prob_pois(const double m, const int len){
  int n;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(-m); 
  for(n=1;n<=len-1;++n) {
    prob[n]=prob[n-1]*m/n;
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> prob_pois(const mpfr_20 m, const int len){
  int n;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(-m); 
  for(n=1;n<=len-1;++n) {
    prob[n]=prob[n-1]*m/n;
  };
  
  return prob;
}

// not exported
std::vector<double> prob_pois_2(const double lambda, const int len){
  int n;
  std::vector<double> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(-lambda*lambda); 
  for(n=1;n<=len-1;++n) {
    prob[n]=prob[n-1]*lambda*lambda/n;
  };
  
  return prob;
}

// not exported
std::vector<mpfr_20> prob_pois_2(const mpfr_20 lambda, const int len){
  int n;
  std::vector<mpfr_20> prob(len);
  
  for(n=1;n<=len-1;++n) {
    prob[n]=0;
  };
  
  prob[0]=exp(-lambda*lambda); 
  for(n=1;n<=len-1;++n) {
    prob[n]=prob[n-1]*lambda*lambda/n;
  };
  
  return prob;
}

// not exported
std::vector<double> prob_pois_deriv1(const double m, std::vector<double> &prob, const int len){
  int n;
  std::vector<double> prob1(len);
  
  for(n=0;n<=len-1;++n){
    prob1[n]=prob[n]*(n/m-1.0);
  };
  
  return prob1;
}

// not exported
std::vector<mpfr_20> prob_pois_deriv1(const mpfr_20 m, std::vector<mpfr_20> &prob, const int len){
  int n;
  std::vector<mpfr_20> prob1(len);
  
  for(n=0;n<=len-1;++n){
    prob1[n]=prob[n]*(n/m-1.0);
  };
  
  return prob1;
}

// not exported
std::vector<double> prob_pois_deriv1_2(const double lambda, std::vector<double> &prob, const int len){
  int n;
  std::vector<double> prob1(len);
  
  for(n=0;n<=len-1;++n){
    prob1[n]=prob[n]*(2*n/lambda-2*lambda);
  };
  
  return prob1;
}

// not exported
std::vector<mpfr_20> prob_pois_deriv1_2(const mpfr_20 lambda, std::vector<mpfr_20> &prob, const int len){
  int n;
  std::vector<mpfr_20> prob1(len);
  
  for(n=0;n<=len-1;++n){
    prob1[n]=prob[n]*(2*n/lambda-2*lambda);
  };
  
  return prob1;
}

// not exported
std::vector<double> prob_pois_deriv2(const double m, std::vector<double> &prob, const int len){
  int n;
  std::vector<double> prob2(len);
  
  for(n=0;n<=len-1;++n){
    prob2[n]=prob[n]*(n*(n-1)/pow(m, 2.0)-2.0*n/m+1.0);
  };
  
  return prob2;
}

// not exported
std::vector<mpfr_20> prob_pois_deriv2(const mpfr_20 m, std::vector<mpfr_20> &prob, const int len){
  int n;
  std::vector<mpfr_20> prob2(len);
  
  for(n=0;n<=len-1;++n){
    prob2[n]=prob[n]*(n*(n-1)/pow(m, 2.0)-2.0*n/m+1.0);
  };
  
  return prob2;
}

// not exported
std::vector<double> prob_pois_deriv2_2(const double lambda, std::vector<double> &prob, const int len){
  int n;
  std::vector<double> prob2(len);
  
  for(n=0;n<=len-1;++n){
    prob2[n]=prob[n]*(4*n*n/lambda/lambda-2.0+4*lambda*lambda-2*n/lambda/lambda-8*n);
  };
  
  return prob2;
}

// not exported
std::vector<mpfr_20> prob_pois_deriv2_2(const mpfr_20 lambda, std::vector<mpfr_20> &prob, const int len){
  int n;
  std::vector<mpfr_20> prob2(len);
  
  for(n=0;n<=len-1;++n){
    prob2[n]=prob[n]*(4*n*n/lambda/lambda-2.0+4*lambda*lambda-2*n/lambda/lambda-8*n);
  };
  
  return prob2;
}

/// GENERAL FUNCTION FOR CALCULATING PMF - NO DERIVATIVES ///

// [[Rcpp::export]]
std::vector<double> prob_mutations(double m, int n, double e=1, double w=1, double cv=0, double death=0, double lag=0, double phi=0, double poisson=0){
  std::vector<double> seq(n);
  std::vector<double> prob(n);

  if (death==0 & lag==0 & phi==0) {
    seq=aux_seq(e, w, n-1);
  } else if (lag!=0 & w==1 & death==0 & phi==0) {
    seq=aux_seq_lag_ext(pow(2, lag), e, n-1);
  } else if (death!=0 & lag==0 & phi==0) {
    seq=aux_seq_death_ext(e, w, death/(1.-death), n-1);
  } else {
    seq=aux_seq_integrate(e, w, death/(1.-death), lag, phi, n-1);
  };
  
  if (cv==0) {
    prob=prob_ld(m, seq, n);
  } else {
    double k=1/cv/cv;
    double A=m/k;
    std::vector<double> xi(n);
    xi=xi_seq(A, seq, n);
    prob=prob_b0(A, k, seq[0], xi, n);
  };
  
  if (poisson==0) {
    return prob;
  } else {
    std::vector<double> prob_p(n);
    prob_p=prob_pois(poisson, n);
    std::vector<double> prob_ldp(n);
    prob_ldp=prob_ld_deriv(prob, prob_p, n);
    return prob_ldp;
  };
}

/// PMF AND DERIVATIVES FOR FOLD ///

// [[Rcpp::export]]
Rcpp::List calc_probs(const double m, const double cv, std::vector<double> &seq, const int len, const double poisson=0, bool boost=false){
  std::vector<double> prob(len),prob1(len),prob2(len);
  
  if (cv==0) {
    if (boost==false){
      if (poisson==0){
        prob=prob_ld(m, seq, len);
      } else {
        std::vector<double> prob_lc(len), prob_p(len);
        prob_lc=prob_ld(m, seq, len);
        prob_p=prob_pois(poisson, len);
        prob=prob_ld_deriv(prob_lc, prob_p, len);
      };
      prob1=prob_ld_deriv(seq, prob, len);
      prob2=prob_ld_deriv(seq, prob1, len);
    } else {
      std::vector<mpfr_20> xprob(len),xprob1(len),xprob2(len);
      mpfr_20 xm=m;
      int i;
      
      if (poisson==0){
        xprob=prob_ld(xm, seq, len);
      } else {
        std::vector<mpfr_20> xprob_lc(len), xprob_p(len);
        mpfr_20 xpoisson=poisson;
        xprob_lc=prob_ld(xm, seq, len);
        xprob_p=prob_pois(xpoisson, len);
        xprob=prob_ld_deriv(xprob_lc, xprob_p, len);
      };
      xprob1=prob_ld_deriv(seq, xprob, len);
      xprob2=prob_ld_deriv(seq, xprob1, len);
      
      for (i=0;i<=len-1;++i){
        prob[i]=xprob[i].convert_to<double>();
        prob1[i]=xprob1[i].convert_to<double>();
        prob2[i]=xprob2[i].convert_to<double>();
      };
    };
  } else {
    if (boost==false){
      double k=1.0/cv/cv;
      std::vector<double> xi(len);
      double seq0=seq[0];
      double A=m/k;
      xi=xi_seq(A, seq, len);
      if (poisson==0){
        prob=prob_b0(A, k, seq0, xi, len);
        prob1=prob_b0_deriv1(A, k, seq, xi, len);
        prob2=prob_b0_deriv2(A, k, seq, xi, len);
      } else {
        std::vector<double> prob_lc(len), prob_lc_1(len), prob_lc_2(len), prob_p(len);
        prob_lc=prob_b0(A, k, seq0, xi, len);
        prob_p=prob_pois(poisson, len);
        prob=prob_ld_deriv(prob_lc, prob_p, len);
        prob_lc_1=prob_b0_deriv1(A, k, seq, xi, len);
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, len);
        prob1=prob_ld_deriv(prob_lc_1, prob_p, len);
        prob2=prob_ld_deriv(prob_lc_2, prob_p, len);
      };
    } else {
      std::vector<mpfr_20> xprob(len),xprob1(len),xprob2(len);
      mpfr_20 k=1.0/cv/cv;
      std::vector<mpfr_20> xi(len);
      mpfr_20 seq0=seq[0];
      mpfr_20 A=m/k;
      int i;
      xi=xi_seq(A, seq, len);
      if (poisson==0){
        xprob=prob_b0(A, k, seq0, xi, len);
        xprob1=prob_b0_deriv1(A, k, seq, xi, len);
        xprob2=prob_b0_deriv2(A, k, seq, xi, len);
      } else {
        std::vector<mpfr_20> xprob_lc(len), xprob_lc_1(len), xprob_lc_2(len), xprob_p(len);
        mpfr_20 xpoisson=poisson;
        xprob_lc=prob_b0(A, k, seq0, xi, len);
        xprob_p=prob_pois(xpoisson, len);
        xprob=prob_ld_deriv(xprob_lc, xprob_p, len);
        xprob_lc_1=prob_b0_deriv1(A, k, seq, xi, len);
        xprob_lc_2=prob_b0_deriv2(A, k, seq, xi, len);
        xprob1=prob_ld_deriv(xprob_lc_1, xprob_p, len);
        xprob2=prob_ld_deriv(xprob_lc_2, xprob_p, len);
      };
      for (i=0;i<=len-1;++i){
        prob[i]=xprob[i].convert_to<double>();
        prob1[i]=xprob1[i].convert_to<double>();
        prob2[i]=xprob2[i].convert_to<double>();
      };
    };
  };
  // Rcout << std::nextafter(1.0, std::numeric_limits<float>::max) << " " << std::nextafter(1.0e-5, std::numeric_limits<float>::max) << " " << std::nextafter(10.0, std::numeric_limits<float>::max) << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("prob") = prob,
                            Rcpp::Named("prob1") = prob1,
                            Rcpp::Named("prob2") = prob2);
}

// [[Rcpp::export]]
Rcpp::List calc_probs_2(const double m, const double cv, std::vector<double> &seq, const int len, bool boost=false){
  std::vector<double> prob(len),prob1(len),prob2(len);
  
  if (cv==0) {
    if (boost==false){
      prob=prob_ld(m, seq, len);
      prob1=prob_ld_deriv(seq, prob, len);
      prob2=prob_ld_deriv(seq, prob1, len);
    } else {
      std::vector<mpfr_20> xprob(len),xprob1(len),xprob2(len);
      mpfr_20 xm=m;
      int i;
      
      xprob=prob_ld(xm, seq, len);
      xprob1=prob_ld_deriv(seq, xprob, len);
      xprob2=prob_ld_deriv(seq, xprob1, len);
      
      for (i=0;i<=len-1;++i){
        prob[i]=xprob[i].convert_to<double>();
        prob1[i]=xprob1[i].convert_to<double>();
        prob2[i]=xprob2[i].convert_to<double>();
      };
    };
  } else {
    if (boost==false){
      double k=1.0/cv/cv;
      std::vector<double> xi(len);
      double seq0=seq[0];
      double A=m/k;
      xi=xi_seq(A, seq, len);
      prob=prob_b0(A, k, seq0, xi, len);
      prob1=prob_b0_deriv1(A, k, seq, xi, len);
      prob2=prob_b0_deriv2(A, k, seq, xi, len);
    } else {
      std::vector<mpfr_20> xprob(len),xprob1(len),xprob2(len);
      mpfr_20 k=1.0/cv/cv;
      std::vector<mpfr_20> xi(len);
      mpfr_20 seq0=seq[0];
      mpfr_20 A=m/k;
      int i;
      xi=xi_seq(A, seq, len);
      xprob=prob_b0(A, k, seq0, xi, len);
      xprob1=prob_b0_deriv1(A, k, seq, xi, len);
      xprob2=prob_b0_deriv2(A, k, seq, xi, len);
      for (i=0;i<=len-1;++i){
        prob[i]=xprob[i].convert_to<double>();
        prob1[i]=xprob1[i].convert_to<double>();
        prob2[i]=xprob2[i].convert_to<double>();
      };
    };
  };
  
  return Rcpp::List::create(Rcpp::Named("prob") = prob,
                            Rcpp::Named("prob1") = prob1,
                            Rcpp::Named("prob2") = prob2);
}

// [[Rcpp::export]]
double logprob(double m, int len, std::vector<int>& data, std::vector<double>& seq, double k=0, double poisson=0){
  int i,number;
  int samplesize=data.size();
  double sumlogprob=0;
  std::vector<double> prob(len);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    sumlogprob+=log(prob[number]);
  };
  
  return sumlogprob;
}

// [[Rcpp::export]]
double logprob_boost(double xm, int len, std::vector<int>& data, std::vector<double>& seq, double xk=0, double xpoisson=0){
  int i,number;
  int samplesize=data.size();
  mpfr_20 m=xm;
  mpfr_20 sumlogprob=0;
  std::vector<mpfr_20> prob(len);
  
  if (xk==0){
    if (xpoisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      mpfr_20 poisson=xpoisson;
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 k=xk;
    mpfr_20 A=m/k;
    xi=xi_seq(A, seq, len);
    if (xpoisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      mpfr_20 poisson=xpoisson;
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    sumlogprob+=log(prob[number]);
  };
  
  return sumlogprob.convert_to<double>();
}

/// HELPER FUNCTIONS FOR ONE-PARAMETER ESTIMATES ///

void derivatives(double &U, double &J, double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, const bool &fisher=true){
  int i,number;
  int samplesize=data.size();
  std::vector<double> prob(len),prob1(len),prob2(len),subscore(samplesize);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
    prob1=prob_ld_deriv(seq, prob, len);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, len);};
  } else {
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
      prob1=prob_b0_deriv1(A, k, seq, xi, len);
      if (fisher) {prob2=prob_b0_deriv2(A, k, seq, xi, len);};
    } else {
      std::vector<double> prob_p(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
      prob_lc_1=prob_b0_deriv1(A, k, seq, xi, len);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,len);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(A, k, seq, xi, len);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,len);
      };
    };
  };
  
  U=0;
  J=0;
  loglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    subscore[i]=(prob1[number]/prob[number]);
    U+=subscore[i];
    if (fisher) {J+=subscore[i]*subscore[i]-prob2[number]/prob[number];};
    loglik+=log(prob[number]);
  };
}

void derivatives_boost(double &U, double &J, double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, const bool &fisher=true){
  int i,number;
  int samplesize=data.size();
  std::vector<mpfr_20> prob(len),prob1(len),prob2(len),subscore(samplesize);
  mpfr_20 xm=m,xk=k,xpoisson=poisson,xU=0,xJ=0,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, len);
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_ld(xm, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
    prob1=prob_ld_deriv(seq, prob, len);
    if (fisher) {prob2=prob_ld_deriv(seq, prob1, len);};
  } else {
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 xA=xm/xk;
    xi=xi_seq(xA, seq, len);
    if (poisson==0){
      prob=prob_b0(xA, xk, seq0, xi, len);
      prob1=prob_b0_deriv1(xA, xk, seq, xi, len);
      if (fisher) {
        prob2=prob_b0_deriv2(xA, xk, seq, xi, len);
      };
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_b0(xA, xk, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
      prob_lc_1=prob_b0_deriv1(xA, xk, seq, xi, len);
      prob1=prob_ld_deriv(prob_p,prob_lc_1,len);
      if (fisher) {
        prob_lc_2=prob_b0_deriv2(xA, xk, seq, xi, len);
        prob2=prob_ld_deriv(prob_p,prob_lc_2,len);
      };
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    subscore[i]=(prob1[number]/prob[number]);
    xU+=subscore[i];
    if (fisher) {xJ+=subscore[i]*subscore[i]-prob2[number]/prob[number];};
    xloglik+=log(prob[number]);
  };
  
  U=static_cast<double>(xU);
  J=static_cast<double>(xJ);
  loglik=static_cast<double>(xloglik);
}

void loglik(double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson){
  int i,number;
  int samplesize=data.size();
  std::vector<double> prob(len);
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(m, seq, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_ld(m, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<double> xi(len);
    double seq0=seq[0];
    double A=m/k;
    xi=xi_seq(A, seq, len);
    if (poisson==0){
      prob=prob_b0(A, k, seq0, xi, len);
    } else {
      std::vector<double> prob_p(len),prob_lc(len);
      prob_p=prob_pois(poisson, len);
      prob_lc=prob_b0(A, k, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  loglik=0;
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    loglik+=log(prob[number]);
  };
}

void loglik_boost(double &loglik, const double &m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson){
  int i,number;
  int samplesize=data.size();
  std::vector<mpfr_20> prob(len);
  mpfr_20 xk=k,xpoisson=poisson,xm=m,xloglik=0;
  
  if (k==0){
    if (poisson==0){
      prob=prob_ld(xm, seq, len);
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_ld(xm, seq, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  } else {
    std::vector<mpfr_20> xi(len);
    mpfr_20 seq0=seq[0];
    mpfr_20 xA=xm/xk;
    xi=xi_seq(xA, seq, len);
    if (poisson==0){
      prob=prob_b0(xA, xk, seq0, xi, len);
    } else {
      std::vector<mpfr_20> prob_p(len),prob_lc(len);
      prob_p=prob_pois(xpoisson, len);
      prob_lc=prob_b0(xA, xk, seq0, xi, len);
      prob=prob_ld_deriv(prob_p,prob_lc,len);
    };
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    xloglik+=log(prob[number]);
  };
  
  loglik=static_cast<double>(xloglik);
}

// [[Rcpp::export]]
std::vector<double> optim_m(double current_m, double lower_m, double upper_m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_U,lower_J,lower_loglik,upper_U,upper_J,upper_loglik;
  
  derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);
  if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);};
  
  derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);
  if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);};
  
  if(std::abs(upper_U)<1e-9){
    derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);
    if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);};
    return std::vector<double>{upper_m,upper_U,upper_J,upper_loglik};
  };
  
  if(std::abs(lower_U)<1e-9){
    derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);
    if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);};
    return std::vector<double>{lower_m,lower_U,lower_J,lower_loglik};
  };
  
  while(upper_U>0){
    lower_m=upper_m;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_m*=5;
    derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);
    if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, false);};
    if(std::abs(upper_U)<1e-9){
      derivatives(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);
      if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_boost(upper_U, upper_J, upper_loglik, upper_m, seq, len, data, k, poisson, true);};
      return std::vector<double>{upper_m,upper_U,upper_J,upper_loglik};
    };
  };
  while(lower_U<0){
    upper_m=lower_m;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_m/=10;
    if (lower_m<1e-24) {return {1e-20};};
    derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);
    if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, false);};
    if(std::abs(lower_U)<1e-9){
      derivatives(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);
      if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_boost(lower_U, lower_J, lower_loglik, lower_m, seq, len, data, k, poisson, true);};
      return std::vector<double>{lower_m,lower_U,lower_J,lower_loglik};
    };
  };
  
  if(verbose) {Rcout << "boundaries:: m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "boundaries:: U: " << lower_U << " " << upper_U << "\n";};
  if(verbose) {Rcout << "boundaries:: loglik: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m<lower_m) || (current_m>upper_m)) {current_m=(lower_m+upper_m)/2;};
  
  derivatives(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {
    derivatives_boost(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  };
  if(verbose) {Rcout<< "U: " << current_U << " J: " << current_J << " loglik: " << current_loglik << " m: " << current_m << "\n\n";};
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {return {-1.0};};
  
  while(std::abs(current_U)>1.0e-9){
    iter++;
    if (iter>50) {return {-1.0};};
    
    new_m=current_m+current_U/current_J;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    };
    
    derivatives(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {derivatives_boost(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);};
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return {-1.0};};
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if (new_U<0){
      upper_m=new_m;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_m=new_m;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_m=new_m;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  return std::vector<double>{current_m, current_U, current_J, current_loglik};
}

// [[Rcpp::export]]
double root_m(double current_m, double lower_m, double upper_m, std::vector<double> &seq, const int &len, std::vector<int> &data, const double &k, const double &poisson, const double lalpha, bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_loglik,upper_loglik;
  
  loglik(upper_loglik, upper_m, seq, len, data, k, poisson);
  if(!std::isfinite(upper_loglik)) {loglik_boost(upper_loglik, upper_m, seq, len, data, k, poisson);};
  upper_loglik-=lalpha;
  
  loglik(lower_loglik, lower_m, seq, len, data, k, poisson);
  if(!std::isfinite(lower_loglik)) {loglik_boost(lower_loglik, lower_m, seq, len, data, k, poisson);};
  lower_loglik-=lalpha;
  
  if (std::abs(lower_loglik)<1e-6){return lower_m;};
  if (std::abs(upper_loglik)<1e-6){return upper_m;};
  
  if (upper_loglik*lower_loglik>0){
    if (upper_loglik<lower_loglik){
      int n=0;
      while(upper_loglik>0){
        n++;
        if (n>100) {return -1;}
        lower_m=upper_m;
        lower_loglik=upper_loglik;
        
        upper_m*=5;
        loglik(upper_loglik, upper_m, seq, len, data, k, poisson);
        if(!std::isfinite(upper_loglik)) {loglik_boost(upper_loglik, upper_m, seq, len, data, k, poisson);};
        upper_loglik-=lalpha;
      };
    } else {
      int n=0;
      while(lower_loglik>0){
        n++;
        if (n>100) {return -1;}
        upper_m=lower_m;
        upper_loglik=lower_loglik;
        
        lower_m/=10;
        if (lower_m<1e-24) {return 0;};
        loglik(lower_loglik, lower_m, seq, len, data, k, poisson);
        if(!std::isfinite(lower_loglik)) {loglik_boost(lower_loglik, lower_m, seq, len, data, k, poisson);};
        lower_loglik-=lalpha;
      };
    };
  };
  
  if (std::abs(lower_loglik)<1e-6){return lower_m;};
  if (std::abs(upper_loglik)<1e-6){return upper_m;};
  
  if(verbose) {Rcout << "Boundaries for m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "Boundary log-likelihood: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m>upper_m) || (current_m<lower_m)) {current_m=(lower_m+upper_m)/2;};
  
  derivatives(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {
    derivatives_boost(current_U, current_J, current_loglik, current_m, seq, len, data, k, poisson);
  };
  current_loglik-=lalpha;
  if(verbose) {Rcout << "Starting m: " << current_m << " Starting U: " << current_U << " Starting J: " << current_J << " Starting loglik: " << current_loglik << "\n";};
  
  while(std::abs(current_loglik)>1e-6){
    iter++;
    if (iter>50) {return -1.0;};
    
    new_m=current_m-current_loglik/current_U;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
      
      derivatives(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);
      if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {derivatives_boost(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);};
      new_loglik-=lalpha;
      
      if (((lower_loglik<upper_loglik) && (new_loglik<0)) || ((lower_loglik>upper_loglik) && (new_loglik>0))) {
        lower_m=new_m;
        lower_loglik=new_loglik;
      } else if (((lower_loglik<upper_loglik) && (new_loglik>0)) || ((lower_loglik>upper_loglik) && (new_loglik<0))) {
        upper_m=new_m;
        upper_loglik=new_loglik;
      };
    } else {
      derivatives(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);
      if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {derivatives_boost(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);};
      new_loglik-=lalpha;
    };
    
    derivatives(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {derivatives_boost(new_U, new_J, new_loglik, new_m, seq, len, data, k, poisson);};
    new_loglik-=lalpha;
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return -1.0;};
    
    current_m=new_m;
    current_U=new_U;
    current_loglik=new_loglik;
    
    // if ((new_loglik<std::min(lower_loglik, upper_loglik)) || (new_loglik>std::max(lower_loglik, upper_loglik))){
    //   
    // };
  };
  return current_m;
}

// [[Rcpp::export]]
std::vector<double> combo_optim_m(double current_m, double lower_m, double upper_m, const double &R, std::vector<double> &seq1, std::vector<double> &seq2,
                                  const int &len1, const int &len2, std::vector<int> &data1, std::vector<int> &data2,
                                  const double &k1, const double &k2, const double &poisson1, const double &poisson2, bool verbose=false){
  int iter=0;
  double new_m,new_U,new_J,new_loglik,current_U,current_J,current_loglik,lower_U,lower_J,lower_loglik,upper_U,upper_J,upper_loglik,
  U1,J1,loglik1,U2,J2,loglik2;
  
  derivatives(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);};
  derivatives(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);};
  upper_U=U1+R*U2;
  upper_J=J1+R*R*J2;
  upper_loglik=loglik1+loglik2;
  if (std::abs(upper_U)<1e-9) {return std::vector<double>{upper_m,upper_loglik};};
  
  derivatives(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);};
  derivatives(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);};
  lower_U=U1+R*U2;
  lower_J=J1+R*R*J2;
  lower_loglik=loglik1+loglik2;
  if (std::abs(lower_U)<1e-9) {return std::vector<double>{lower_m,lower_loglik};};
  
  while(upper_U>0){
    lower_m=upper_m;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_m*=5;
    derivatives(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, upper_m, seq1, len1, data1, k1, poisson1, true);};
    derivatives(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*upper_m, seq2, len2, data2, k2, poisson2, true);};
    upper_U=U1+R*U2;
    upper_J=J1+R*R*J2;
    upper_loglik=loglik1+loglik2;
    if (std::abs(upper_U)<1e-9) {return std::vector<double>{upper_m,upper_loglik};};
  };
  while(lower_U<0){
    upper_m=lower_m;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_m/=10;
    if (lower_m<1e-24) {return {0};};
    derivatives(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, lower_m, seq1, len1, data1, k1, poisson1, true);};
    derivatives(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*lower_m, seq2, len2, data2, k2, poisson2, true);};
    lower_U=U1+R*U2;
    lower_J=J1+R*R*J2;
    lower_loglik=loglik1+loglik2;
    if (std::abs(lower_U)<1e-9) {return std::vector<double>{lower_m,lower_loglik};};
  };
  
  if(verbose) {Rcout << "boundaries:: m: " << lower_m << " " << upper_m << "\n";};
  if(verbose) {Rcout << "boundaries:: U: " << lower_U << " " << upper_U << "\n";};
  if(verbose) {Rcout << "boundaries:: loglik: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if ((current_m<lower_m) || (current_m>upper_m)) {current_m=(lower_m+upper_m)/2;};
  
  derivatives(U1, J1, loglik1, current_m, seq1, len1, data1, k1, poisson1, true);
  if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, current_m, seq1, len1, data1, k1, poisson1, true);};
  derivatives(U2, J2, loglik2, R*current_m, seq2, len2, data2, k2, poisson2, true);
  if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*current_m, seq2, len2, data2, k2, poisson2, true);};
  current_U=U1+R*U2;
  current_J=J1+R*R*J2;
  current_loglik=loglik1+loglik2;
  
  if(verbose) {Rcout<< "U: " << current_U << " J: " << current_J << " loglik: " << current_loglik << " m: " << current_m << "\n\n";};
  if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {return {-1.0};};
  
  while(std::abs(current_U)>1.0e-9){
    iter++;
    if (iter>50) {return {-1.0};};
    
    new_m=current_m+current_U/current_J;
    if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
    
    if((new_m<lower_m) || (new_m>upper_m)){
      new_m=(lower_m+upper_m)/2;
      if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
    };
    
    derivatives(U1, J1, loglik1, new_m, seq1, len1, data1, k1, poisson1, true);
    if(!std::isfinite(U1) || !std::isfinite(J1) || !std::isfinite(loglik1)) {derivatives_boost(U1, J1, loglik1, new_m, seq1, len1, data1, k1, poisson1, true);};
    derivatives(U2, J2, loglik2, R*new_m, seq2, len2, data2, k2, poisson2, true);
    if(!std::isfinite(U2) || !std::isfinite(J2) || !std::isfinite(loglik2)) {derivatives_boost(U2, J2, loglik2, R*new_m, seq2, len2, data2, k2, poisson2, true);};
    new_U=U1+R*U2;
    new_J=J1+R*R*J2;
    new_loglik=loglik1+loglik2;
    if ((!std::isfinite(new_U)) || (!std::isfinite(new_J)) || (!std::isfinite(new_loglik))) {return {-1.0};};
    
    if(verbose) {Rcout << "U: " << new_U << " J: " << new_J << " loglik: " << new_loglik << "\n";};
    
    if (new_U<0){
      upper_m=new_m;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_m=new_m;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_m=new_m;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  return std::vector<double>{current_m, current_loglik};
}

// not exported
void derivatives_joint(double &U_1, double &U_2, double &J_11, double &J_12, double &J_22, double &logprob,
                       const double &sqrtm, const double &e, const double &param, const std::vector<int> &data, const int &len,
                       const std::string &option){
  int i,number;
  int samplesize=data.size();
  std::vector<double> prob(len),prob_10(len),prob_20(len),prob_01(len),prob_02(len),prob_11(len);
  
  if ((option=="lag") || (option=="growth")) {
    std::vector<double> seq(len),seq1(len),seq2(len);
    
    if (option=="lag") {
      double L=pow(2,param*param);
      double log2=0.693147180559945;
      double dLdparam=L*log2*2*param;
      double d2Ldparam2=dLdparam*(log2*2*param+1/param);
      double dLdparam2=dLdparam*dLdparam;
      seq=aux_seq_lag(L,e,len-1);
      seq1=aux_seq_lag_deriv1(L,e,len-1);
      seq2=aux_seq_lag_deriv2(L,e,len-1);
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*dLdparam2+seq1[i]*d2Ldparam2;
        seq1[i]=seq1[i]*dLdparam;
      }; // derivatives with respect to a (number of generations of lag)
    } else {
      double h=1e-4;
      std::vector<double> seqplus(len), seqminus(len);
      volatile double whp=param*param*(1.0+h);
      volatile double whm=param*param*(1.0-h);
      double dwp=whp-param*param;
      double dwm=param*param-whm;
      double dw=whp-whm;
      seq=aux_seq_boost(e, param*param, len-1);
      seqplus=aux_seq_boost(e, whp, len-1);
      seqminus=aux_seq_boost(e, whm, len-1);
      for (i=0;i<len;++i) {
        seq1[i]=(seqplus[i]-seqminus[i])/dw;
        seq2[i]=(seqplus[i]+seqminus[i]-2*seq[i])/dwp/dwm;
      };
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*2*param*2*param+seq1[i]*2;
        seq1[i]=seq1[i]*2*param;
      };
    };
    
    prob=prob_ld(sqrtm*sqrtm, seq, len);
    prob_10=prob_ld_deriv(seq, prob, len);
    prob_20=prob_ld_deriv(seq, prob_10, len);
    prob_01=prob_ld_deriv(seq1, prob, len);
    for (i=0;i<=len-1;++i){
      prob_20[i]=prob_20[i]*4*sqrtm*sqrtm+prob_10[i]*2;
      prob_10[i]*=2*sqrtm;
      prob_01[i]*=sqrtm*sqrtm;
    };
    std::vector<double> seq1_prob_01(len);
    seq1_prob_01=prob_ld_deriv(seq1, prob_01, len);
    std::vector<double> seq2_prob(len);
    seq2_prob=prob_ld_deriv(seq2, prob, len);
    for (i=0;i<=len-1;++i){
      prob_02[i]=sqrtm*sqrtm*(seq1_prob_01[i]+seq2_prob[i]);
    };
    std::vector<double> seq1_prob_10(len);
    seq1_prob_10=prob_ld_deriv(seq1, prob_10, len);
    for (i=0;i<=len-1;++i){
      prob_11[i]=sqrtm*sqrtm*seq1_prob_10[i]+prob_01[i]/sqrtm/sqrtm;
    };
  } else if (option=="mixed") {
    std::vector<double> auxseq(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len),prob_p(len),prob_p_1(len),prob_p_2(len);
    auxseq=aux_seq(e, 1, len-1);
    prob_lc=prob_ld(sqrtm*sqrtm, auxseq, len);
    prob_lc_1=prob_ld_deriv(prob_lc, auxseq, len);
    prob_lc_2=prob_ld_deriv(prob_lc_1, auxseq, len);
    for (i=0;i<=len-1;++i){
      prob_lc_2[i]=prob_lc_2[i]*4*sqrtm*sqrtm+prob_lc_1[i]*2;
      prob_lc_1[i]*=2*sqrtm;
    };
    prob_p=prob_pois_2(param, len);
    prob_p_1=prob_pois_deriv1_2(param, prob_p, len);
    prob_p_2=prob_pois_deriv2_2(param, prob_p, len);
    
    prob=prob_ld_deriv(prob_lc, prob_p, len);
    prob_10=prob_ld_deriv(prob_lc_1, prob_p, len);
    prob_20=prob_ld_deriv(prob_lc_2, prob_p, len);
    prob_01=prob_ld_deriv(prob_lc, prob_p_1, len);
    prob_02=prob_ld_deriv(prob_lc, prob_p_2, len);
    prob_11=prob_ld_deriv(prob_lc_1, prob_p_1, len);
  };
  
  U_1=0; U_2=0; J_11=0; J_12=0; J_22=0; logprob=0;
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    U_1+=prob_10[number]/prob[number];
    U_2+=prob_01[number]/prob[number];
    J_11+=prob_10[number]/prob[number]*prob_10[number]/prob[number]-prob_20[number]/prob[number];
    J_12+=prob_10[number]*prob_01[number]/prob[number]/prob[number]-prob_11[number]/prob[number];
    J_22+=prob_01[number]/prob[number]*prob_01[number]/prob[number]-prob_02[number]/prob[number];
    logprob+=log(prob[number]);
  };
}

// not exported
void derivatives_joint_boost(double &U_1, double &U_2, double &J_11, double &J_12, double &J_22, double &loglik,
                             const double &sqrtm, const double &e, const double &param, const std::vector<int> &data, const int &len,
                             const std::string &option){
  int i,number;
  int samplesize=data.size();
  std::vector<mpfr_20> prob(len),prob_10(len),prob_20(len),prob_01(len),prob_02(len),prob_11(len);
  mpfr_20 xsqrtm=sqrtm, xe=e, xparam=param;
  mpfr_20 xU_1=0, xU_2=0, xJ_11=0, xJ_12=0, xJ_22=0, xloglik=0;
  
  if ((option=="lag") || (option=="growth")) {
    std::vector<mpfr_20> seq(len),seq1(len),seq2(len);
    
    if (option=="lag") {
      mpfr_20 L=pow(2,xparam*xparam);
      mpfr_20 log2=0.693147180559945;
      mpfr_20 dLdparam=L*log2*2*xparam;
      mpfr_20 d2Ldparam2=dLdparam*(log2*2*xparam+1/xparam);
      mpfr_20 dLdparam2=dLdparam*dLdparam;
      seq=aux_seq_lag(L,xe,len-1);
      seq1=aux_seq_lag_deriv1(L,xe,len-1);
      seq2=aux_seq_lag_deriv2(L,xe,len-1);
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*dLdparam2+seq1[i]*d2Ldparam2;
        seq1[i]=seq1[i]*dLdparam;
      }; // derivatives with respect to a (number of generations of lag)
    } else {
      mpfr_20 h=1e-4;
      std::vector<mpfr_20> seqplus(len), seqminus(len);
      mpfr_20 whp=xparam*xparam*(1.0+h);
      mpfr_20 whm=xparam*xparam*(1.0-h);
      mpfr_20 dwp=whp-xparam*xparam;
      mpfr_20 dwm=xparam*xparam-whm;
      mpfr_20 dw=whp-whm;
      seq=aux_seq_boost(xe, xparam*xparam, len-1);
      seqplus=aux_seq_boost(xe, whp, len-1);
      seqminus=aux_seq_boost(xe, whm, len-1);
      for (i=0;i<len;++i) {
        seq1[i]=(seqplus[i]-seqminus[i])/dw;
        seq2[i]=(seqplus[i]+seqminus[i]-2*seq[i])/dwp/dwm;
      };
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*2*xparam*2*xparam+seq1[i]*2;
        seq1[i]=seq1[i]*2*xparam;
      };
    };
    
    prob=prob_ld(xsqrtm*xsqrtm, seq, len);
    prob_10=prob_ld_deriv(seq, prob, len);
    prob_20=prob_ld_deriv(seq, prob_10, len);
    prob_01=prob_ld_deriv(seq1, prob, len);
    for (i=0;i<=len-1;++i){
      prob_20[i]=prob_20[i]*4*xsqrtm*xsqrtm+prob_10[i]*2;
      prob_10[i]*=2*xsqrtm;
      prob_01[i]*=xsqrtm*xsqrtm;
    };
    std::vector<mpfr_20> seq1_prob_01(len);
    seq1_prob_01=prob_ld_deriv(seq1, prob_01, len);
    std::vector<mpfr_20> seq2_prob(len);
    seq2_prob=prob_ld_deriv(seq2, prob, len);
    for (i=0;i<=len-1;++i){
      prob_02[i]=xsqrtm*xsqrtm*(seq1_prob_01[i]+seq2_prob[i]);
    };
    std::vector<mpfr_20> seq1_prob_10(len);
    seq1_prob_10=prob_ld_deriv(seq1, prob_10, len);
    for (i=0;i<=len-1;++i){
      prob_11[i]=xsqrtm*xsqrtm*seq1_prob_10[i]+prob_01[i]/xsqrtm/xsqrtm;
    };
  } else if (option=="mixed") {
    std::vector<mpfr_20> auxseq(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len),prob_p(len),prob_p_1(len),prob_p_2(len);
    auxseq=aux_seq_boost(xe, 1, len-1);
    prob_lc=prob_ld(xsqrtm*xsqrtm, auxseq, len);
    prob_lc_1=prob_ld_deriv(prob_lc, auxseq, len);
    prob_lc_2=prob_ld_deriv(prob_lc_1, auxseq, len);
    for (i=0;i<=len-1;++i){
      prob_lc_2[i]=prob_lc_2[i]*4*xsqrtm*xsqrtm+prob_lc_1[i]*2;
      prob_lc_1[i]*=2*xsqrtm;
    };
    prob_p=prob_pois_2(xparam, len);
    prob_p_1=prob_pois_deriv1_2(xparam, prob_p, len);
    prob_p_2=prob_pois_deriv2_2(xparam, prob_p, len);
    
    prob=prob_ld_deriv(prob_lc, prob_p, len);
    prob_10=prob_ld_deriv(prob_lc_1, prob_p, len);
    prob_20=prob_ld_deriv(prob_lc_2, prob_p, len);
    prob_01=prob_ld_deriv(prob_lc, prob_p_1, len);
    prob_02=prob_ld_deriv(prob_lc, prob_p_2, len);
    prob_11=prob_ld_deriv(prob_lc_1, prob_p_1, len);
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    xU_1+=prob_10[number]/prob[number];
    xU_2+=prob_01[number]/prob[number];
    xJ_11+=prob_10[number]/prob[number]*prob_10[number]/prob[number]-prob_20[number]/prob[number];
    xJ_12+=prob_10[number]*prob_01[number]/prob[number]/prob[number]-prob_11[number]/prob[number];
    xJ_22+=prob_01[number]/prob[number]*prob_01[number]/prob[number]-prob_02[number]/prob[number];
    xloglik+=log(prob[number]);
  };
  
  U_1=static_cast<double>(xU_1); U_2=static_cast<double>(xU_2); J_11=static_cast<double>(xJ_11);
  J_12=static_cast<double>(xJ_12); J_22=static_cast<double>(xJ_22); loglik=static_cast<double>(xloglik);
}

// not exported
std::vector<double> optim_m_const_param(double current_m, double lower_m, double upper_m, const double &e, double &param,
                                        const int &len, std::vector<int> &data, const std::string &option) {
  std::vector<double> seq(len);
  std::vector<double> res(4);
  double poisson=0, U1=0, U2=0, J11=0, J12=0, J22=0, loglik=0;
  
  if (option == "lag") {
    seq=aux_seq_lag(pow(2,param*param), e, len-1);
  } else if (option == "growth") {
    seq=aux_seq_boost(e, param*param, len-1);
  } else if (option == "mixed") {
    seq=aux_seq(e, 1, len-1);
    poisson=param*param;
  };
  
  res=optim_m(current_m*current_m, lower_m*lower_m, upper_m*upper_m, seq, len, data, 0.0, poisson);
  current_m=sqrt(res[0]);
  
  derivatives_joint(U1, U2, J11, J12, J22, loglik, current_m, e, param, data, len, option);
  if ((!std::isfinite(U1)) || (!std::isfinite(U2)) || (!std::isfinite(J11)) || (!std::isfinite(J12)) || (!std::isfinite(J22)) || (!std::isfinite(loglik))) {
    derivatives_joint_boost(U1, U2, J11, J12, J22, loglik, current_m, e, param, data, len, option);
  };
  
  return std::vector<double> {current_m, U1, U2, J11, J12, J22, loglik};
}

// not exported
void derivatives_param(double &U_2, double &J_22, double &logprob, const double &sqrtm, const double &e, const double &param,
                       const std::vector<int> &data, const int &len, const std::string &option){
  int i,number;
  int samplesize=data.size();
  std::vector<double> prob(len),prob_01(len),prob_02(len);
  
  if ((option=="lag") || (option=="growth")) {
    std::vector<double> seq(len),seq1(len),seq2(len);
    
    if (option=="lag") {
      double L=pow(2,param*param);
      double log2=0.693147180559945;
      double dLdparam=L*log2*2*param;
      double d2Ldparam2=dLdparam*(log2*2*param+1/param);
      double dLdparam2=dLdparam*dLdparam;
      seq=aux_seq_lag(L,e,len-1);
      seq1=aux_seq_lag_deriv1(L,e,len-1);
      seq2=aux_seq_lag_deriv2(L,e,len-1);
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*dLdparam2+seq1[i]*d2Ldparam2;
        seq1[i]=seq1[i]*dLdparam;
      };
    } else {
      double h=1e-4;
      std::vector<double> seqplus(len), seqminus(len);
      volatile double whp=param*param*(1.0+h);
      volatile double whm=param*param*(1.0-h);
      double dwp=whp-param*param;
      double dwm=param*param-whm;
      double dw=whp-whm;
      seq=aux_seq_boost(e, param*param, len-1);
      seqplus=aux_seq_boost(e, whp, len-1);
      seqminus=aux_seq_boost(e, whm, len-1);
      for (i=0;i<len;++i) {
        seq1[i]=(seqplus[i]-seqminus[i])/dw;
        seq2[i]=(seqplus[i]+seqminus[i]-2*seq[i])/dwp/dwm;
      };
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*2*param*2*param+seq1[i]*2;
        seq1[i]=seq1[i]*2*param;
      };
    };
    
    prob=prob_ld(sqrtm*sqrtm, seq, len);
    prob_01=prob_ld_deriv(seq1, prob, len);
    for (i=0;i<=len-1;++i){
      prob_01[i]*=sqrtm*sqrtm;
    };
    std::vector<double> seq1_prob_01(len);
    seq1_prob_01=prob_ld_deriv(seq1, prob_01, len);
    std::vector<double> seq2_prob(len);
    seq2_prob=prob_ld_deriv(seq2, prob, len);
    for (i=0;i<=len-1;++i){
      prob_02[i]=sqrtm*sqrtm*(seq1_prob_01[i]+seq2_prob[i]);
    };
  } else if (option=="mixed") {
    std::vector<double> auxseq(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len),prob_p(len),prob_p_1(len),prob_p_2(len);
    auxseq=aux_seq(e, 1, len-1);
    prob_lc=prob_ld(sqrtm*sqrtm, auxseq, len);
    prob_p=prob_pois_2(param, len);
    prob_p_1=prob_pois_deriv1_2(param, prob_p, len);
    prob_p_2=prob_pois_deriv2_2(param, prob_p, len);
    
    prob=prob_ld_deriv(prob_lc, prob_p, len);
    prob_01=prob_ld_deriv(prob_lc, prob_p_1, len);
    prob_02=prob_ld_deriv(prob_lc, prob_p_2, len);
  };
  
  U_2=0; J_22=0; logprob=0;
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    U_2+=prob_01[number]/prob[number];
    J_22+=prob_01[number]/prob[number]*prob_01[number]/prob[number]-prob_02[number]/prob[number];
    logprob+=log(prob[number]);
  };
}

// not exported
void derivatives_param_boost(double &U_2, double &J_22, double &loglik, const double &sqrtm, const double &e, const double &param,
                             const std::vector<int> &data, const int &len, const std::string &option){
  int i,number;
  int samplesize=data.size();
  std::vector<mpfr_20> prob(len),prob_01(len),prob_02(len);
  mpfr_20 xsqrtm=sqrtm, xe=e, xparam=param;
  mpfr_20 xU_2=0, xJ_22=0, xloglik=0;
  
  if ((option=="lag") || (option=="growth")) {
    std::vector<mpfr_20> seq(len),seq1(len),seq2(len);
    
    if (option=="lag") {
      mpfr_20 L=pow(2,xparam*xparam);
      mpfr_20 log2=0.693147180559945;
      mpfr_20 dLdparam=L*log2*2*xparam;
      mpfr_20 d2Ldparam2=dLdparam*(log2*2*xparam+1/xparam);
      mpfr_20 dLdparam2=dLdparam*dLdparam;
      seq=aux_seq_lag(L,xe,len-1);
      seq1=aux_seq_lag_deriv1(L,xe,len-1);
      seq2=aux_seq_lag_deriv2(L,xe,len-1);
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*dLdparam2+seq1[i]*d2Ldparam2;
        seq1[i]=seq1[i]*dLdparam;
      };
    } else {
      mpfr_20 h=1e-4;
      std::vector<mpfr_20> seqplus(len), seqminus(len);
      mpfr_20 whp=xparam*xparam*(1.0+h);
      mpfr_20 whm=xparam*xparam*(1.0-h);
      mpfr_20 dwp=whp-xparam*xparam;
      mpfr_20 dwm=xparam*xparam-whm;
      mpfr_20 dw=whp-whm;
      seq=aux_seq_boost(xe, xparam*xparam, len-1);
      seqplus=aux_seq_boost(xe, whp, len-1);
      seqminus=aux_seq_boost(xe, whm, len-1);
      for (i=0;i<len;++i) {
        seq1[i]=(seqplus[i]-seqminus[i])/dw;
        seq2[i]=(seqplus[i]+seqminus[i]-2*seq[i])/dwp/dwm;
      };
      for (i=0;i<len;++i) {
        seq2[i]=seq2[i]*2*xparam*2*xparam+seq1[i]*2;
        seq1[i]=seq1[i]*2*xparam;
      };
    };
    
    prob=prob_ld(xsqrtm*xsqrtm, seq, len);
    prob_01=prob_ld_deriv(seq1, prob, len);
    for (i=0;i<=len-1;++i){
      prob_01[i]*=xsqrtm*xsqrtm;
    };
    std::vector<mpfr_20> seq1_prob_01(len);
    seq1_prob_01=prob_ld_deriv(seq1, prob_01, len);
    std::vector<mpfr_20> seq2_prob(len);
    seq2_prob=prob_ld_deriv(seq2, prob, len);
    for (i=0;i<=len-1;++i){
      prob_02[i]=xsqrtm*xsqrtm*(seq1_prob_01[i]+seq2_prob[i]);
    };
  } else if (option=="mixed") {
    std::vector<mpfr_20> auxseq(len),prob_lc(len),prob_lc_1(len),prob_lc_2(len),prob_p(len),prob_p_1(len),prob_p_2(len);
    auxseq=aux_seq_boost(xe, 1, len-1);
    prob_lc=prob_ld(xsqrtm*xsqrtm, auxseq, len);
    prob_p=prob_pois_2(xparam, len);
    prob_p_1=prob_pois_deriv1_2(xparam, prob_p, len);
    prob_p_2=prob_pois_deriv2_2(xparam, prob_p, len);
    
    prob=prob_ld_deriv(prob_lc, prob_p, len);
    prob_01=prob_ld_deriv(prob_lc, prob_p_1, len);
    prob_02=prob_ld_deriv(prob_lc, prob_p_2, len);
  };
  
  for(i=0;i<=samplesize-1;++i){
    number=data[i];
    xU_2+=prob_01[number]/prob[number];
    xJ_22+=prob_01[number]/prob[number]*prob_01[number]/prob[number]-prob_02[number]/prob[number];
    xloglik+=log(prob[number]);
  };
  
  U_2=static_cast<double>(xU_2); J_22=static_cast<double>(xJ_22); loglik=static_cast<double>(xloglik);
}

// not exported
std::vector<double> optim_param_const_m(double current_param, double lower_param, double upper_param, const double &e, double &m,
                                        const int &len, std::vector<int> &data, const std::string &option) {
  double U1=0, U2=0, J11=0, J12=0, J22=0, loglik=0;
  bool exit=false;
  int iter=0;
  double new_param,new_U,new_J,new_loglik,current_U=0,current_J,current_loglik,lower_U,lower_J,lower_loglik,upper_U,upper_J,upper_loglik;
  
  derivatives_param(upper_U, upper_J, upper_loglik, m, e, upper_param, data, len, option);
  if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_param_boost(upper_U, upper_J, upper_loglik, m, e, upper_param, data, len, option);};
  
  derivatives_param(lower_U, lower_J, lower_loglik, m, e, lower_param, data, len, option);
  if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_param_boost(lower_U, lower_J, lower_loglik, m, e, lower_param, data, len, option);};
  
  if(std::abs(upper_U)<1e-9){
    current_param = upper_param;
    exit = true;
  };
  
  if(std::abs(lower_U)<1e-9 && exit==false){
    current_param = lower_param;
    exit = true;
  };
  
  while(upper_U>0 && exit==false){
    lower_param=upper_param;lower_U=upper_U;lower_J=upper_J;lower_loglik=upper_loglik;
    upper_param *= 2.25;
    derivatives_param(upper_U, upper_J, upper_loglik, m, e, upper_param, data, len, option);
    if(!std::isfinite(upper_U) || !std::isfinite(upper_J) || !std::isfinite(upper_loglik)) {derivatives_param_boost(upper_U, upper_J, upper_loglik, m, e, upper_param, data, len, option);};
    if(std::abs(upper_U) < 1e-9 && exit==false){
      current_param = upper_param;
      exit = true;
    };
  };
  while(lower_U<0 && exit==false){
    upper_param=lower_param;upper_U=lower_U;upper_J=lower_J;upper_loglik=lower_loglik;
    lower_param /= 3;
    derivatives_param(lower_U, lower_J, lower_loglik, m, e, lower_param, data, len, option);
    if(!std::isfinite(lower_U) || !std::isfinite(lower_J) || !std::isfinite(lower_loglik)) {derivatives_param_boost(lower_U, lower_J, lower_loglik, m, e, lower_param, data, len, option);};
    if((std::abs(lower_U) < 1e-9 || lower_param < 1e-6) && exit==false){
      current_param = lower_param;
      exit = true;
    };
  };
  
  
  if (exit==false) {
    if ((current_param<lower_param) || (current_param>upper_param)) {current_param=(lower_param+upper_param)/2;};
    
    derivatives_param(current_U, current_J, current_loglik, m, e, current_param, data, len, option);
    if(!std::isfinite(current_U) || !std::isfinite(current_J) || !std::isfinite(current_loglik)) {derivatives_param_boost(current_U, current_J, current_loglik, m, e, current_param, data, len, option);};
    
    if ((!std::isfinite(current_U)) || (!std::isfinite(current_J)) || (!std::isfinite(current_loglik))) {
      current_param = -1.0;
      exit = true;
    };
  }
  
  while(std::abs(current_U)>1.0e-9 && exit==false){
    iter++;
    if (iter>50) {
      current_param = -1.0;
      exit = true;  
    };
    
    new_param=current_param+current_U/current_J;
    
    if((new_param<lower_param) || (new_param>upper_param)){
      new_param=(lower_param+upper_param)/2;
    };
    
    derivatives_param(new_U, new_J, new_loglik, m, e, new_param, data, len, option);
    if(!std::isfinite(new_U) || !std::isfinite(new_J) || !std::isfinite(new_loglik)) {derivatives_param_boost(new_U, new_J, new_loglik, m, e, new_param, data, len, option);};
    if(!std::isfinite(new_U) || !std::isfinite(new_J) || !std::isfinite(new_loglik)) {
      current_param = -1.0;
      exit = true;  
    };
    
    if (new_U<0){
      upper_param=new_param;upper_U=new_U;upper_J=new_J;upper_loglik=new_loglik;
    } else if (new_U>0){
      lower_param=new_param;lower_U=new_U;lower_J=new_J;lower_loglik=new_loglik;
    };
    current_param=new_param;current_U=new_U;current_J=new_J;current_loglik=new_loglik;
  };
  
  derivatives_joint(U1, U2, J11, J12, J22, loglik, m, e, current_param, data, len, option);
  if ((!std::isfinite(U1)) || (!std::isfinite(U2)) || (!std::isfinite(J11)) || (!std::isfinite(J12)) || (!std::isfinite(J22)) || (!std::isfinite(loglik))) {
    derivatives_joint_boost(U1, U2, J11, J12, J22, loglik, m, e, current_param, data, len, option);
  };
  
  return std::vector<double> {current_param, U1, U2, J11, J12, J22, loglik};
}

// not exported
void calc_deltas(double &U_1, double &U_2, double &J_11, double &J_12, double &J_22, double &delta_m, double &delta_param){
  double invdet=1.0/(J_11*J_22-J_12*J_12);
  delta_m=J_22*invdet*U_1-J_12*invdet*U_2;
  delta_param=-J_12*invdet*U_1+J_11*invdet*U_2;
};

int sgn(double x){
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

double l2norm(double U1, double U2){
  return sqrt(U1*U1 + U2*U2);
}

// not exported
void eigenvalue(double &eigen1, double &eigen2, double &J_11, double &J_12, double &J_22){
  double trace=J_11+J_22;
  double det=J_11*J_22-J_12*J_12;
  double sq=sqrt(trace*trace-4*det);
  eigen1=(trace+sq)/2;
  eigen2=(trace-sq)/2;
}

// not exported
void eigenvector(double &eigen11, double &eigen12, double &eigen21, double &eigen22, double &J_11, double &J_12, double &J_22){
  double tau=(J_22-J_11)/2/J_12;
  double t=sgn(tau)/(std::abs(tau)+sqrt(tau*tau+1));
  eigen11=1/sqrt(t*t+1);
  eigen22=eigen11;
  eigen12=eigen11*t;
  eigen21=-eigen12;
}

// not exported
void decompose(double &J_11, double &J_12, double &J_22, double &delta) {
  double eigen1, eigen2, eigen11, eigen12, eigen21, eigen22;
  double first_J_11, first_J_12, first_J_21, first_J_22;
  double second_J_11, second_J_12, second_J_21, second_J_22;
  eigenvalue(eigen1, eigen2, J_11, J_12, J_22);
  eigenvector(eigen11, eigen12, eigen21, eigen22, J_11, J_12, J_22);
  if (std::abs(eigen1) < delta) {
    eigen1=delta;
  } else {
    eigen1=std::abs(eigen1);
  };
  if (std::abs(eigen2) < delta) {
    eigen2=delta;
  } else {
    eigen2=std::abs(eigen2);
  };
  first_J_11=eigen11*eigen1;
  first_J_12=eigen12*eigen2;
  first_J_21=eigen21*eigen1;
  first_J_22=eigen22*eigen2;
  double invdet=1/(eigen11*eigen22-eigen12*eigen21);
  second_J_11=first_J_11*invdet*eigen22-first_J_12*eigen21*invdet;
  second_J_12=-first_J_11*invdet*eigen12+first_J_12*eigen11*invdet;
  second_J_21=first_J_21*invdet*eigen11-first_J_22*eigen21*invdet;
  second_J_22=-first_J_21*invdet*eigen12+first_J_22*eigen11*invdet;

  J_11=second_J_11;
  J_12=second_J_12;
  J_22=second_J_22;
}

// [[Rcpp::export]]
std::vector<double> optim_m_param_joint(double current_m, double lower_m, double upper_m, double current_param, double lower_param, double upper_param,
                                        const double e, std::vector<int> &data, const int len, const std::string &option, bool verbose=false){
  std::vector<double> upper(6), lower(6), current(6);
  double upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik;
  double lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik;
  double current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik;
  double new_m, new_param, new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik;
  double eigen1, eigen2;
  double delta_m, delta_param;
  int iter=0;
  
  if (verbose) {Rcout << "***HYBRID BISECTION ALGORITHM***\n";};
  
  upper=optim_m_const_param(upper_m, upper_m/2, upper_m*2, e, upper_param, len, data, option);
  upper_m=upper[0]; upper_U1=upper[1]; upper_U2=upper[2]; upper_J11=upper[3]; upper_J12=upper[4]; upper_J22=upper[5]; upper_loglik=upper[6];
  
  lower=optim_m_const_param(lower_m, lower_m/2, lower_m*2, e, lower_param, len, data, option);
  lower_m=lower[0]; lower_U1=lower[1]; lower_U2=lower[2]; lower_J11=lower[3]; lower_J12=lower[4]; lower_J22=lower[5]; lower_loglik=lower[6];
  
  if (verbose) {Rcout << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
  if (verbose) {Rcout << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
  
  if(std::abs(sqrt(upper_U1 * upper_U1 + upper_U2 * upper_U2)) < 1e-6){
    eigenvalue(eigen1, eigen2, upper_J11, upper_J12, upper_J22);
    if((eigen1 >= 0) && (eigen2 >= 0)){
      return std::vector<double> {upper_m*upper_m, upper_param*upper_param, upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik};
    };
  };
  
  if(std::abs(sqrt(lower_U1 * lower_U1 + lower_U2 * lower_U2)) < 1e-6){
    eigenvalue(eigen1, eigen2, lower_J11, lower_J12, lower_J22);
    if((eigen1 >= 0) && (eigen2 >= 0)){
      return std::vector<double> {lower_m*lower_m, lower_param*lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik};
    };
  };
  
  while(upper_U2 > 0){
    lower_param=upper_param; lower_m=upper_m; lower_U1=upper_U1; lower_U2=upper_U2; lower_J11=upper_J11; lower_J12=upper_J12; lower_J22=upper_J22; lower_loglik=upper_loglik;
    upper_param*=1.7;
    upper=optim_m_const_param(upper_m, upper_m/2, upper_m*2, e, upper_param, len, data, option);
    upper_m=upper[0]; upper_U1=upper[1]; upper_U2=upper[2]; upper_J11=upper[3]; upper_J12=upper[4]; upper_J22=upper[5]; upper_loglik=upper[6];
    if (verbose) {Rcout << "boundaries were readjusted\n" << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
    if(std::abs(sqrt(upper_U1 * upper_U1 + upper_U2 * upper_U2)) < 1e-6){
      eigenvalue(eigen1, eigen2, upper_J11, upper_J12, upper_J22);
      if((eigen1 >= 0) && (eigen2 >= 0)){
        return std::vector<double> {upper_m*upper_m, upper_param*upper_param, upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik};
      };
    };
  };
  
  while(lower_U2 < 0){
    upper_param=lower_param; upper_m=lower_m; upper_U1=lower_U1; upper_U2=lower_U2; upper_J11=lower_J11; upper_J12=lower_J12; upper_J22=lower_J22; upper_loglik=lower_loglik;
    lower_param/=1.7;
    lower=optim_m_const_param(lower_m, lower_m/2, lower_m*2, e, lower_param, len, data, option);
    lower_m=lower[0]; lower_U1=lower[1]; lower_U2=lower[2]; lower_J11=lower[3]; lower_J12=lower[4]; lower_J22=lower[5]; lower_loglik=lower[6];
    if (verbose) {Rcout << "boundaries were readjusted\n" << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
    if (lower_param < 1e-6){
      return std::vector<double> {lower_m*lower_m, lower_param*lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik};
    };
    if(std::abs(sqrt(lower_U1 * lower_U1 + lower_U2 * lower_U2)) < 1e-6){
      eigenvalue(eigen1, eigen2, lower_J11, lower_J12, lower_J22);
      if((eigen1 >= 0) && (eigen2 >= 0)){
        return std::vector<double> {lower_m*lower_m, lower_param*lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik};
      };
    };
  };
  
  if (verbose) {Rcout << "Set parameter bounds to " << lower_param << " and " << upper_param << "\n";};
  
  if (verbose) {Rcout << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
  if (verbose) {Rcout << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
  
  if (lower_U1 < 1e-3 && upper_U1 < 1e-3 && lower_U2 < 1e-3 && upper_U2 < 1e-3 && std::abs(upper_loglik-lower_loglik) < 1e-3) {
    if (verbose) {Rcout << "param value might be 0\n";};
    current_param = 1e-10;
    current=optim_m_const_param(current_m, current_m/2, current_m*2, e, current_param, len, data, option);
    current_m=current[0]; current_U1=current[1]; current_U2=current[2]; current_J11=current[3]; current_J12=current[4]; current_J22=current[5]; current_loglik=current[6];
    if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
    if (std::abs(current_loglik-lower_loglik) < 1e-3) {
      return std::vector<double> {current_m*current_m, current_param*current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik};
    };
  };
  
  if (current_param<lower_param || current_param>upper_param){
    current_param=(lower_param+upper_param)/2;
    current_m=(lower_m+upper_m)/2;
  };
  derivatives_joint(current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik, current_m, e, current_param, data, len, option);
  if ((!std::isfinite(current_U1)) || (!std::isfinite(current_U2)) || (!std::isfinite(current_J11)) || (!std::isfinite(current_J12)) || (!std::isfinite(current_J22)) || (!std::isfinite(current_loglik))) {
    derivatives_joint_boost(current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik, current_m, e, current_param, data, len, option);
  };
  if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
  
  while(std::abs(sqrt(current_U1 * current_U1 + current_U2 * current_U2)) > 1e-6) {
    iter++;
    if (iter > 50) {return {-1.0};};
    
    if (verbose) {Rcout << std::endl << "Iteration " << iter << " with guesses m=" << current_m << " and parameter=" << current_param << "\n";};
    
    calc_deltas(current_U1, current_U2, current_J11, current_J12, current_J22, delta_m, delta_param);
    new_m=current_m+delta_m;
    new_param=current_param+delta_param;
    
    if ((new_param < lower_param) || (new_param > upper_param) || (new_m < 0)) {
      if (verbose) {Rcout << "Joint m and parameter: Bisection step" << "\n";}
      for (int i=1; i<=2; ++i){
        new_param=(lower_param+upper_param)/2;
        current=optim_m_const_param(current_m, current_m/2, current_m*2, e, new_param, len, data, option);
        new_m=current[0]; new_U1=current[1]; new_U2=current[2]; new_J11=current[3]; new_J12=current[4]; new_J22=current[5]; new_loglik=current[6];
        if (new_U2 > 0) {
          lower_m=new_m; lower_param=new_param; lower_U1=new_U1; lower_U2=new_U2; lower_J11=new_J11; lower_J12=new_J12; lower_J22=new_J22; lower_loglik=new_loglik;
        } else if (new_U2 < 0) {
          upper_m=new_m; upper_param=new_param; upper_U1=new_U1; upper_U2=new_U2; upper_J11=new_J11; upper_J12=new_J12; upper_J22=new_J22; upper_loglik=new_loglik;
        };
      }
    } else {
      if (verbose) {Rcout << "Joint m and parameter: Newton step" << "\n";}
      derivatives_joint(new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik, new_m, e, new_param, data, len, option);
      if ((!std::isfinite(new_U1)) || (!std::isfinite(new_U2)) || (!std::isfinite(new_J11)) || (!std::isfinite(new_J12)) || (!std::isfinite(new_J22)) || (!std::isfinite(new_loglik))) {
        derivatives_joint_boost(new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik, new_m, e, new_param, data, len, option);
      };
    };
    
    eigenvalue(eigen1, eigen2, new_J11, new_J12, new_J22);
    if (verbose) {Rcout << "Eigenvalues are " << eigen1 << " and " << eigen2 << "\n";};
    
    if (verbose) {Rcout << "Difference is " << std::abs(new_m-current_m)/current_m + std::abs(new_param-current_param)/current_param << "\n";};
    if (std::abs(new_m-current_m)/current_m + std::abs(new_param-current_param)/current_param < 1e-6) {
      return std::vector<double> {current_m*current_m, current_param*current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik};
    };
    
    current_m=new_m; current_param=new_param; current_U1=new_U1; current_U2=new_U2; current_J11=new_J11; current_J12=new_J12; current_J22=new_J22; current_loglik=new_loglik;
    if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
    
  };
  
  return std::vector<double> {current_m*current_m, current_param*current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik};
  
}

// [[Rcpp::export]]
std::vector<double> optim_m_param_joint_2(double current_m, double current_param, const double e, std::vector<int> &data,
                                          const int len, const std::string &option, bool verbose=false){
  double current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik;
  double new_m, new_param, new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik;
  double eigen1=-1, eigen2=-1;
  double delta_m, delta_param;
  int iter=0;
  double delta=0.0005;
  
  if (verbose) {Rcout << "***SUFFICIENT DESCENT ALGORITHM***\n";};
  
  derivatives_joint(current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik, current_m, e, current_param, data, len, option);
  if ((!std::isfinite(current_U1)) || (!std::isfinite(current_U2)) || (!std::isfinite(current_J11)) || (!std::isfinite(current_J12)) || (!std::isfinite(current_J22)) || (!std::isfinite(current_loglik))) {
    derivatives_joint_boost(current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik, current_m, e, current_param, data, len, option);
  };
  eigenvalue(eigen1, eigen2, current_J11, current_J12, current_J22);
  if (verbose) {Rcout << "Eigenvalues are " << eigen1 << " and " << eigen2 << "\n";};
  
  while((std::abs(l2norm(current_U1, current_U2)) > 1e-6) || (eigen1 < 0) || (eigen2 < 0)) {
    iter++;
    if (iter > 50) {return {-1.0};};
    
    if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
    
    if (verbose) {Rcout << std::endl << "Iteration " << iter << " with guesses m=" << current_m << " and parameter=" << current_param << "\n";};
    
    if ((eigen1 < 0) || (eigen2 < 0)) {
      decompose(current_J11, current_J12, current_J22, delta);
      if (verbose) {Rcout << "Setting new Hessian: " << "current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << "\n";};
      eigenvalue(eigen1, eigen2, current_J11, current_J12, current_J22);
      if (verbose) {Rcout << "New eigenvalues are " << eigen1 << " and " << eigen2 << "\n";};
    };
    calc_deltas(current_U1, current_U2, current_J11, current_J12, current_J22, delta_m, delta_param);
    new_m=current_m+delta_m;
    new_param=current_param+delta_param;
    double iter2=0;
    while((new_m < 0) || (new_param < 0)) {
      iter2++;
      if (iter > 50) {return {-1.0};};
      delta_m*=0.5;
      delta_param*=0.5;
      new_m=current_m+delta_m;
      new_param=current_param+delta_param;
    };
    derivatives_joint(new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik, new_m, e, new_param, data, len, option);
    if ((!std::isfinite(new_U1)) || (!std::isfinite(new_U2)) || (!std::isfinite(new_J11)) || (!std::isfinite(new_J12)) || (!std::isfinite(new_J22)) || (!std::isfinite(new_loglik))) {
      if (verbose) {Rcout << "boost" << std::endl;};
      derivatives_joint_boost(new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik, new_m, e, new_param, data, len, option);
    };
    if (verbose) {Rcout << "new_m: " << new_m << " new_param: " << new_param << " new_U1: " << new_U1 << " new_U2: " << new_U2 << " new_J11: " << new_J11 << " new_J12: " << new_J12 << " new_J22: " << new_J22 << " new_loglik: " << new_loglik << "\n";};
    
    eigenvalue(eigen1, eigen2, new_J11, new_J12, new_J22);
    if (verbose) {Rcout << "Eigenvalues are " << eigen1 << " and " << eigen2 << "\n";};
    
    current_m=new_m; current_param=new_param; current_U1=new_U1; current_U2=new_U2; current_J11=new_J11; current_J12=new_J12; current_J22=new_J22; current_loglik=new_loglik;
    
  };
  
  return std::vector<double> {current_m*current_m, current_param*current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik};
  
}

// [[Rcpp::export]]
std::vector<double> root_m_param_joint(double current_m, double lower_m, double upper_m, double current_param, double lower_param, double upper_param,
                                       const double e, std::vector<int> &data, const int len, const double lalpha, const int &rootparam,
                                       const int &lowerroot, const std::string &option, bool verbose=false) {
  
  std::vector<double> upper(6), lower(6), current(6), newx(6);
  double upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik;
  double lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik;
  double current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik;
  double new_m, new_param, new_U1, new_U2, new_J11, new_J12, new_J22, new_loglik;
  int iter=0;
  int n=0;
  
  switch(rootparam){
  case 1:
    if (verbose) {Rcout << "rootparam: " << 1 << "\n";}
    upper=optim_m_const_param(upper_m, upper_m/2, upper_m*2, e, upper_param, len, data, option);
    lower=optim_m_const_param(lower_m, lower_m/2, lower_m*2, e, lower_param, len, data, option);
    upper_m=upper[0]; upper_U1=upper[1]; upper_U2=upper[2]; upper_J11=upper[3]; upper_J12=upper[4]; upper_J22=upper[5]; upper_loglik=upper[6]-lalpha;
    lower_m=lower[0]; lower_U1=lower[1]; lower_U2=lower[2]; lower_J11=lower[3]; lower_J12=lower[4]; lower_J22=lower[5]; lower_loglik=lower[6]-lalpha;
    break;
  case 0:
    if (verbose) {Rcout << "rootparam: " << 0 << "\n";}
    upper=optim_param_const_m(upper_param, upper_param/2, upper_param*2, e, upper_m, len, data, option);
    lower=optim_param_const_m(lower_param, lower_param/2, lower_param*2, e, lower_m, len, data, option);
    upper_param=upper[0]; upper_U1=upper[1]; upper_U2=upper[2]; upper_J11=upper[3]; upper_J12=upper[4]; upper_J22=upper[5]; upper_loglik=upper[6]-lalpha;
    lower_param=lower[0]; lower_U1=lower[1]; lower_U2=lower[2]; lower_J11=lower[3]; lower_J12=lower[4]; lower_J22=lower[5]; lower_loglik=lower[6]-lalpha;
    break;
  };
  
  if (verbose) {Rcout << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
  if (verbose) {Rcout << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
  
  if (std::abs(lower_loglik) < 1e-6) {return std::vector<double>{lower_m*lower_m, lower_param*lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik};};
  if (std::abs(upper_loglik) < 1e-6) {return std::vector<double>{upper_m*upper_m, upper_param*upper_param, upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik};};
  
  if (upper_loglik*lower_loglik>0){
    switch(lowerroot){
    case 1:
      while(lower_loglik > 0){
        n++;
        if (n > 100) {return std::vector<double>{-1,-1,-1,-1,-1,-1,-1,-1};}
        upper_m=lower_m; upper_param=lower_param; upper_U1=lower_U1; upper_U2=lower_U2; upper_J11=lower_J11; upper_J12=lower_J12; upper_J22=lower_J22; upper_loglik=lower_loglik;
        
        switch(rootparam){
        case 1:
          lower_param /= 3;
          lower=optim_m_const_param(lower_m, lower_m/2, lower_m*2, e, lower_param, len, data, option);
          lower_m=lower[0]; lower_U1=lower[1]; lower_U2=lower[2]; lower_J11=lower[3]; lower_J12=lower[4]; lower_J22=lower[5]; lower_loglik=lower[6]-lalpha;
          break;
        case 0:
          lower_m /= 3;
          lower=optim_param_const_m(lower_param, lower_param/2, lower_param*2, e, lower_m, len, data, option);
          lower_param=lower[0]; lower_U1=lower[1]; lower_U2=lower[2]; lower_J11=lower[3]; lower_J12=lower[4]; lower_J22=lower[5]; lower_loglik=lower[6]-lalpha;
          break;
        };
        if (verbose) {Rcout << "boundaries were readjusted\n" << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
        if (lower_m < 1e-6) {return std::vector<double>{lower_m*lower_m, lower_param*lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik};};
        
      };
      break;
    case 0:
      while(upper_loglik > 0){
        n++;
        if (verbose) {Rcout << "n: " << n << "\n";}
        if (n > 100) {return std::vector<double>{-1,-1,-1,-1,-1,-1,-1,-1};}
        lower_m=upper_m; lower_param=upper_param; lower_U1=upper_U1; lower_U2=upper_U2; lower_J11=upper_J11; lower_J12=upper_J12; lower_J22=upper_J22; lower_loglik=upper_loglik;
        
        switch(rootparam){
        case 1:
          upper_param *= 2.25;
          upper=optim_m_const_param(upper_m, upper_m/2, upper_m*2, e, upper_param, len, data, option);
          upper_m=upper[0]; upper_U1=upper[1]; upper_U2=upper[2]; upper_J11=upper[3]; upper_J12=upper[4]; upper_J22=upper[5]; upper_loglik=upper[6]-lalpha;
          break;
        case 0:
          upper_m *= 2.25;
          upper=optim_param_const_m(upper_param, upper_param/2, upper_param*2, e, upper_m, len, data, option);
          upper_param=upper[0]; upper_U1=upper[1]; upper_U2=upper[2]; upper_J11=upper[3]; upper_J12=upper[4]; upper_J22=upper[5]; upper_loglik=upper[6]-lalpha;
          break;
        };
        if (verbose) {Rcout << "boundaries were readjusted\n" << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
        
      };
      break;
    };
  };
  
  if (std::abs(lower_loglik) < 1e-6) {return std::vector<double>{lower_m, lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik};};
  if (std::abs(upper_loglik) < 1e-6) {return std::vector<double>{upper_m, upper_param, upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik};};
  
  if (verbose) {Rcout << "Set parameter bounds to " << lower_param << " and " << upper_param << "\n";};
  if (verbose) {Rcout << "Boundary log-likelihood: " << lower_loglik << " " << upper_loglik << "\n";};
  
  if (verbose) {Rcout << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
  if (verbose) {Rcout << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
  
  if ((rootparam == 1 && (current_param<lower_param || current_param>upper_param)) || (rootparam == 0 && (current_m<lower_m || current_m>upper_m))){
    current_param=(lower_param+upper_param)/2;
    current_m=(lower_m+upper_m)/2;
  };
  switch(rootparam){
  case 1:
    current=optim_m_const_param(current_m, current_m/2, current_m*2, e, current_param, len, data, option);
    current_m=current[0]; current_U1=current[1]; current_U2=current[2]; current_J11=current[3]; current_J12=current[4]; current_J22=current[5]; current_loglik=current[6]-lalpha;
    break;
  case 0:
    current=optim_param_const_m(current_param, current_param/2, current_param*2, e, current_m, len, data, option);
    current_param=current[0]; current_U1=current[1]; current_U2=current[2]; current_J11=current[3]; current_J12=current[4]; current_J22=current[5]; current_loglik=current[6]-lalpha;
    break;
  };
  
  if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
  
  while (std::abs(current_loglik) > 1e-6) {
    iter++;
    if (iter>50) {return std::vector<double>{-1,-1,-1,-1,-1,-1,-1,-1};};
    
    switch(rootparam){
    case 1:
      new_param=current_param-current_loglik/current_U2;
      if(verbose) {Rcout << "newton:: iter: " << iter << " param: " << new_param << "\n";};
      break;
    case 0:
      new_m=current_m-current_loglik/current_U1;
      if(verbose) {Rcout << "newton:: iter: " << iter << " m: " << new_m << "\n";};
      break;
    };
    
    if((((new_m<lower_m) || (new_m>upper_m)) && rootparam==0) || (((new_param<lower_param) || (new_param>upper_param)) && rootparam==1)){
      
      switch(rootparam){
      case 1:
        new_param=(lower_param+upper_param)/2;
        if(verbose) {Rcout << "bisection:: param: " << new_param << "\n";};
        newx=optim_m_const_param(new_m, new_m/2, new_m*2, e, new_param, len, data, option);
        new_m=newx[0]; new_U1=newx[1]; new_U2=newx[2]; new_J11=newx[3]; new_J12=newx[4]; new_J22=newx[5]; new_loglik=newx[6]-lalpha;
        break;
      case 0:
        new_m=(lower_m+upper_m)/2;
        if(verbose) {Rcout << "bisection:: m: " << new_m << "\n";};
        newx=optim_param_const_m(new_param, new_param/2, new_param*2, e, new_m, len, data, option);
        new_param=newx[0]; new_U1=newx[1]; new_U2=newx[2]; new_J11=newx[3]; new_J12=newx[4]; new_J22=newx[5]; new_loglik=newx[6]-lalpha;
        break;
      };
      
      if (((lower_loglik<upper_loglik) && (new_loglik<0)) || ((lower_loglik>upper_loglik) && (new_loglik>0))) {
        switch(rootparam){
        case 1:
          lower_param=new_param;
          break;
        case 0:
          lower_m=new_m;
          break;
        };
        lower_loglik=new_loglik;
      } else if (((lower_loglik<upper_loglik) && (new_loglik>0)) || ((lower_loglik>upper_loglik) && (new_loglik<0))) {
        switch(rootparam){
        case 1:
          upper_param=new_param;
          break;
        case 0:
          upper_m=new_m;
          break;
        };       
        upper_loglik=new_loglik;
      };
      
    } else {
      
      switch(rootparam){
      case 1:
        newx=optim_m_const_param(new_m, new_m/2, new_m*2, e, new_param, len, data, option);
        new_m=newx[0]; new_U1=newx[1]; new_U2=newx[2]; new_J11=newx[3]; new_J12=newx[4]; new_J22=newx[5]; new_loglik=newx[6]-lalpha;
        break;
      case 0:
        newx=optim_param_const_m(new_param, new_param/2, new_param*2, e, new_m, len, data, option);
        new_param=newx[0]; new_U1=newx[1]; new_U2=newx[2]; new_J11=newx[3]; new_J12=newx[4]; new_J22=newx[5]; new_loglik=newx[6]-lalpha;
        break;
      };
      
    };
    
    current_m=new_m; current_param=new_param; current_U1=new_U1; current_U2=new_U2; current_J11=new_J11; current_J12=new_J12; current_J22=new_J22; current_loglik=new_loglik;
    if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
    
  }
  
  return std::vector<double>{current_m*current_m, current_param*current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik};
  
}

// [[Rcpp::export]]
std::vector<double> lower_root_star_m_param_joint(std::vector<double> max, const double e, std::vector<int> &data, const int len,
                                                  double alpha, const int &rootparam, const std::string &option, bool verbose=false) {
  
  double max_m=sqrt(max[0]), max_param=sqrt(max[1]), max_U1=max[2], max_U2=max[3], max_J11=max[4], max_J12=max[5], max_J22=max[6], max_loglik=max[7];
  double boundary_m=max_m, boundary_param=max_param;
  std::vector<double> boundary(6), lower(6), upper(6), current(6);
  double boundary_U1, boundary_U2, boundary_J11, boundary_J12, boundary_J22, boundary_loglik;
  double lower_m, lower_param, lower_U1, lower_U2, lower_J11, lower_J12, lower_J22, lower_loglik;
  double upper_m, upper_param, upper_U1, upper_U2, upper_J11, upper_J12, upper_J22, upper_loglik;
  double current_m, current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik;
  double rb, r, p;
  int iter=0;
  
  switch(rootparam){
  case 1:
    boundary_param = 1e-10;
    boundary=optim_m_const_param(boundary_m, boundary_m/2, boundary_m*2, e, boundary_param, len, data, option);
    boundary_m=boundary[0]; boundary_U1=boundary[1]; boundary_U2=boundary[2]; boundary_J11=boundary[3]; boundary_J12=boundary[4]; boundary_J22=boundary[5]; boundary_loglik=boundary[6];
    break;
  case 0:
    boundary_m = 1e-10;
    boundary=optim_param_const_m(boundary_param, boundary_param/2, boundary_param*2, e, boundary_m, len, data, option);
    boundary_param=boundary[0]; boundary_U1=boundary[1]; boundary_U2=boundary[2]; boundary_J11=boundary[3]; boundary_J12=boundary[4]; boundary_J22=boundary[5]; boundary_loglik=boundary[6];
    break;
  };
  
  rb=sqrt(-2*boundary_loglik + 2*max_loglik);
  if (verbose) {Rcout << "rb: " << rb << " pnorm: " << Rf_pnorm5(-rb, 0, 1, 1, 0) << "\n";};
  
  lower_m=boundary_m; lower_param=boundary_param; lower_U1=boundary_U1; lower_U2=boundary_U2; lower_J11=boundary_J11; lower_J12=boundary_J12; lower_J22=boundary_J22; lower_loglik=boundary_loglik;
  upper_m=max_m; upper_param=max_param; upper_U1=max_U1; upper_U2=max_U2; upper_J11=max_J11; upper_J12=max_J12; upper_J22=max_J22; upper_loglik=max_loglik;
  if (verbose) {Rcout << "upper_m: " << upper_m << " upper_param: " << upper_param << " upper_U1: " << upper_U1 << " upper_U2: " << upper_U2 << " upper_J11: " << upper_J11 << " upper_J12: " << upper_J12 << " upper_J22: " << upper_J22 << " upper_loglik: " << upper_loglik << "\n";};
  if (verbose) {Rcout << "lower_m: " << lower_m << " lower_param: " << lower_param << " lower_U1: " << lower_U1 << " lower_U2: " << lower_U2 << " lower_J11: " << lower_J11 << " lower_J12: " << lower_J12 << " lower_J22: " << lower_J22 << " lower_loglik: " << lower_loglik << "\n";};
  current_m = max_m*2;
  current_param = max_param*2;
  
  do {
    iter++;
    if (verbose) {Rcout << "iter: " << iter << "\n";};
    if (iter>50) {return std::vector<double>{-1,-1,-1,-1,-1,-1,-1,-1};};
    switch(rootparam){
    case 1:
      current_param = (lower_param+upper_param)/2;
      current=optim_m_const_param(current_m, current_m/2, current_m*2, e, current_param, len, data, option);
      current_m=current[0]; current_U1=current[1]; current_U2=current[2]; current_J11=current[3]; current_J12=current[4]; current_J22=current[5]; current_loglik=current[6];
      break;
    case 0:
      current_m = (lower_m+upper_m)/2;
      current=optim_param_const_m(current_param, current_param/2, current_param*2, e, current_m, len, data, option);
      current_param=current[0]; current_U1=current[1]; current_U2=current[2]; current_J11=current[3]; current_J12=current[4]; current_J22=current[5]; current_loglik=current[6];
      break;
    };
    if (verbose) {Rcout << "current_m: " << current_m << " current_param: " << current_param << " current_U1: " << current_U1 << " current_U2: " << current_U2 << " current_J11: " << current_J11 << " current_J12: " << current_J12 << " current_J22: " << current_J22 << " current_loglik: " << current_loglik << "\n";};
    r = sqrt(-2*current_loglik + 2*max_loglik);
    p = Rf_pnorm5(-r, 0, 1, 1, 0) + Rf_pnorm5( -0.5 * (rb - r + r*r/(rb-r) ), 0, 1, 1, 0) - alpha;
    if (p < 0) {
      lower_m=current_m; lower_param=current_param; lower_U1=current_U1; lower_U2=current_U2; lower_J11=current_J11; lower_J12=current_J12; lower_J22=current_J22; lower_loglik=current_loglik;
    } else {
      upper_m=current_m; upper_param=current_param; upper_U1=current_U1; upper_U2=current_U2; upper_J11=current_J11; upper_J12=current_J12; upper_J22=current_J22; upper_loglik=current_loglik;
    }
    if (verbose) {Rcout << "p: " << p << "\n";};
  } while (std::abs(p)>1e-6);
  
  return std::vector<double>{current_m*current_m, current_param*current_param, current_U1, current_U2, current_J11, current_J12, current_J22, current_loglik};
  
}

// [[Rcpp::export]]
std::vector<double> fisher_info(const double sqrtm, const double sqrtparam, const double e, const int len, const std::string &option) {
  double I11=0, I12=0, I22=0;
  int i;
  std::vector<double> prob(len),prob_10(len),prob_01(len);
  
  if ((option=="lag") || (option=="growth")) {
    std::vector<double> seq(len),seq1(len);
    
    if (option=="lag") {
      double L=pow(2,sqrtparam*sqrtparam);
      double log2=0.693147180559945;
      double dLdparam=L*log2*2*sqrtparam;
      seq=aux_seq_lag(L,e,len-1);
      seq1=aux_seq_lag_deriv1(L,e,len-1);
      for (i=0;i<len;++i) {
        seq1[i]=seq1[i]*dLdparam;
      }; // derivatives with respect to a (number of generations of lag)
    } else {
      double h=1e-4;
      std::vector<double> seqplus(len), seqminus(len);
      volatile double whp=sqrtparam*sqrtparam*(1.0+h);
      volatile double whm=sqrtparam*sqrtparam*(1.0-h);
      double dw=whp-whm;
      seq=aux_seq_boost(e, sqrtparam*sqrtparam, len-1);
      seqplus=aux_seq_boost(e, whp, len-1);
      seqminus=aux_seq_boost(e, whm, len-1);
      for (i=0;i<len;++i) {
        seq1[i]=(seqplus[i]-seqminus[i])/dw;
      };
      for (i=0;i<len;++i) {
        seq1[i]=seq1[i]*2*sqrtparam;
      };
    };
    
    prob=prob_ld(sqrtm*sqrtm, seq, len);
    prob_10=prob_ld_deriv(seq, prob, len);
    prob_01=prob_ld_deriv(seq1, prob, len);
    for (i=0;i<=len-1;++i){
      prob_10[i]*=2*sqrtm;
      prob_01[i]*=sqrtm*sqrtm;
    };
  } else if (option=="mixed") {
    std::vector<double> auxseq(len),prob_lc(len),prob_lc_1(len),prob_p(len),prob_p_1(len);
    auxseq=aux_seq(e, 1, len-1);
    prob_lc=prob_ld(sqrtm*sqrtm, auxseq, len);
    prob_lc_1=prob_ld_deriv(prob_lc, auxseq, len);
    for (i=0;i<=len-1;++i){
      prob_lc_1[i]*=2*sqrtm;
    };
    prob_p=prob_pois_2(sqrtparam, len);
    prob_p_1=prob_pois_deriv1_2(sqrtparam, prob_p, len);
    
    prob=prob_ld_deriv(prob_lc, prob_p, len);
    prob_10=prob_ld_deriv(prob_lc_1, prob_p, len);
    prob_01=prob_ld_deriv(prob_lc, prob_p_1, len);
  };
  
  for(i=0;i<=len-1;++i){
    I11+=prob_10[i]*prob_10[i]/prob[i];
    I22+=prob_01[i]*prob_01[i]/prob[i];
    I12+=prob_10[i]*prob_01[i]/prob[i];
  };
  
  return std::vector<double> {I11, I12, I22};
}


