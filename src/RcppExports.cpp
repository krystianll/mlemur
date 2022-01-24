// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rluria_list
Rcpp::List rluria_list(const int n, const double rate, const double N0, const double Nt, const double mut_fit, const int type, const double wt_dp, const double mut_dp, const double lag, const double e, const double cv, const int trim);
RcppExport SEXP _mlemur_rluria_list(SEXP nSEXP, SEXP rateSEXP, SEXP N0SEXP, SEXP NtSEXP, SEXP mut_fitSEXP, SEXP typeSEXP, SEXP wt_dpSEXP, SEXP mut_dpSEXP, SEXP lagSEXP, SEXP eSEXP, SEXP cvSEXP, SEXP trimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< const double >::type N0(N0SEXP);
    Rcpp::traits::input_parameter< const double >::type Nt(NtSEXP);
    Rcpp::traits::input_parameter< const double >::type mut_fit(mut_fitSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const double >::type wt_dp(wt_dpSEXP);
    Rcpp::traits::input_parameter< const double >::type mut_dp(mut_dpSEXP);
    Rcpp::traits::input_parameter< const double >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const int >::type trim(trimSEXP);
    rcpp_result_gen = Rcpp::wrap(rluria_list(n, rate, N0, Nt, mut_fit, type, wt_dp, mut_dp, lag, e, cv, trim));
    return rcpp_result_gen;
END_RCPP
}
// rluria_vec
NumericVector rluria_vec(const int n, const double rate, const double N0, const double Nt, const double mut_fit, const int type, const double wt_dp, const double mut_dp, const double lag, const double e, const double cv, const int trim);
RcppExport SEXP _mlemur_rluria_vec(SEXP nSEXP, SEXP rateSEXP, SEXP N0SEXP, SEXP NtSEXP, SEXP mut_fitSEXP, SEXP typeSEXP, SEXP wt_dpSEXP, SEXP mut_dpSEXP, SEXP lagSEXP, SEXP eSEXP, SEXP cvSEXP, SEXP trimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< const double >::type N0(N0SEXP);
    Rcpp::traits::input_parameter< const double >::type Nt(NtSEXP);
    Rcpp::traits::input_parameter< const double >::type mut_fit(mut_fitSEXP);
    Rcpp::traits::input_parameter< const int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const double >::type wt_dp(wt_dpSEXP);
    Rcpp::traits::input_parameter< const double >::type mut_dp(mut_dpSEXP);
    Rcpp::traits::input_parameter< const double >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< const int >::type trim(trimSEXP);
    rcpp_result_gen = Rcpp::wrap(rluria_vec(n, rate, N0, Nt, mut_fit, type, wt_dp, mut_dp, lag, e, cv, trim));
    return rcpp_result_gen;
END_RCPP
}
// aux_seq_integrate
std::vector<double> aux_seq_integrate(double e, double w, double d, double lag, double phi, int n);
RcppExport SEXP _mlemur_aux_seq_integrate(SEXP eSEXP, SEXP wSEXP, SEXP dSEXP, SEXP lagSEXP, SEXP phiSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_seq_integrate(e, w, d, lag, phi, n));
    return rcpp_result_gen;
END_RCPP
}
// aux_seq
std::vector<double> aux_seq(double e, double w, int n, int option);
RcppExport SEXP _mlemur_aux_seq(SEXP eSEXP, SEXP wSEXP, SEXP nSEXP, SEXP optionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type option(optionSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_seq(e, w, n, option));
    return rcpp_result_gen;
END_RCPP
}
// aux_seq_deriv1_deriv2
Rcpp::List aux_seq_deriv1_deriv2(double e, double w, int n, double h, bool boost);
RcppExport SEXP _mlemur_aux_seq_deriv1_deriv2(SEXP eSEXP, SEXP wSEXP, SEXP nSEXP, SEXP hSEXP, SEXP boostSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< bool >::type boost(boostSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_seq_deriv1_deriv2(e, w, n, h, boost));
    return rcpp_result_gen;
END_RCPP
}
// aux_seq_lag_ext
std::vector<double> aux_seq_lag_ext(double L, double e, int n, bool boost, int deriv);
RcppExport SEXP _mlemur_aux_seq_lag_ext(SEXP LSEXP, SEXP eSEXP, SEXP nSEXP, SEXP boostSEXP, SEXP derivSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type L(LSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type boost(boostSEXP);
    Rcpp::traits::input_parameter< int >::type deriv(derivSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_seq_lag_ext(L, e, n, boost, deriv));
    return rcpp_result_gen;
END_RCPP
}
// aux_seq_death_ext
std::vector<double> aux_seq_death_ext(double e, double w, double d, int n, bool boost);
RcppExport SEXP _mlemur_aux_seq_death_ext(SEXP eSEXP, SEXP wSEXP, SEXP dSEXP, SEXP nSEXP, SEXP boostSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< bool >::type boost(boostSEXP);
    rcpp_result_gen = Rcpp::wrap(aux_seq_death_ext(e, w, d, n, boost));
    return rcpp_result_gen;
END_RCPP
}
// prob_ld
std::vector<double> prob_ld(const double m, std::vector<double>& seq, const int len);
RcppExport SEXP _mlemur_prob_ld(SEXP mSEXP, SEXP seqSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(prob_ld(m, seq, len));
    return rcpp_result_gen;
END_RCPP
}
// xi_seq
std::vector<double> xi_seq(const double A, std::vector<double>& seq, const int len);
RcppExport SEXP _mlemur_xi_seq(SEXP ASEXP, SEXP seqSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(xi_seq(A, seq, len));
    return rcpp_result_gen;
END_RCPP
}
// prob_b0
std::vector<double> prob_b0(const double A, const double k, const double seq0, std::vector<double>& xi, const int len);
RcppExport SEXP _mlemur_prob_b0(SEXP ASEXP, SEXP kSEXP, SEXP seq0SEXP, SEXP xiSEXP, SEXP lenSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type A(ASEXP);
    Rcpp::traits::input_parameter< const double >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double >::type seq0(seq0SEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    rcpp_result_gen = Rcpp::wrap(prob_b0(A, k, seq0, xi, len));
    return rcpp_result_gen;
END_RCPP
}
// prob_mutations
std::vector<double> prob_mutations(double m, int n, double e, double w, double cv, double death, double lag, double phi, double poisson);
RcppExport SEXP _mlemur_prob_mutations(SEXP mSEXP, SEXP nSEXP, SEXP eSEXP, SEXP wSEXP, SEXP cvSEXP, SEXP deathSEXP, SEXP lagSEXP, SEXP phiSEXP, SEXP poissonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< double >::type death(deathSEXP);
    Rcpp::traits::input_parameter< double >::type lag(lagSEXP);
    Rcpp::traits::input_parameter< double >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< double >::type poisson(poissonSEXP);
    rcpp_result_gen = Rcpp::wrap(prob_mutations(m, n, e, w, cv, death, lag, phi, poisson));
    return rcpp_result_gen;
END_RCPP
}
// calc_probs
Rcpp::List calc_probs(const double m, const double cv, std::vector<double>& seq, const int len, const double poisson, bool boost);
RcppExport SEXP _mlemur_calc_probs(SEXP mSEXP, SEXP cvSEXP, SEXP seqSEXP, SEXP lenSEXP, SEXP poissonSEXP, SEXP boostSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< const double >::type poisson(poissonSEXP);
    Rcpp::traits::input_parameter< bool >::type boost(boostSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_probs(m, cv, seq, len, poisson, boost));
    return rcpp_result_gen;
END_RCPP
}
// calc_probs_2
Rcpp::List calc_probs_2(const double m, const double cv, std::vector<double>& seq, const int len, bool boost);
RcppExport SEXP _mlemur_calc_probs_2(SEXP mSEXP, SEXP cvSEXP, SEXP seqSEXP, SEXP lenSEXP, SEXP boostSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< bool >::type boost(boostSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_probs_2(m, cv, seq, len, boost));
    return rcpp_result_gen;
END_RCPP
}
// logprob
double logprob(double m, int len, std::vector<int>& data, std::vector<double>& seq, double k, double poisson);
RcppExport SEXP _mlemur_logprob(SEXP mSEXP, SEXP lenSEXP, SEXP dataSEXP, SEXP seqSEXP, SEXP kSEXP, SEXP poissonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< double >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type poisson(poissonSEXP);
    rcpp_result_gen = Rcpp::wrap(logprob(m, len, data, seq, k, poisson));
    return rcpp_result_gen;
END_RCPP
}
// logprob_boost
double logprob_boost(double xm, int len, std::vector<int>& data, std::vector<double>& seq, double xk, double xpoisson);
RcppExport SEXP _mlemur_logprob_boost(SEXP xmSEXP, SEXP lenSEXP, SEXP dataSEXP, SEXP seqSEXP, SEXP xkSEXP, SEXP xpoissonSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type xm(xmSEXP);
    Rcpp::traits::input_parameter< int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< double >::type xk(xkSEXP);
    Rcpp::traits::input_parameter< double >::type xpoisson(xpoissonSEXP);
    rcpp_result_gen = Rcpp::wrap(logprob_boost(xm, len, data, seq, xk, xpoisson));
    return rcpp_result_gen;
END_RCPP
}
// optim_m
std::vector<double> optim_m(double current_m, double lower_m, double upper_m, std::vector<double>& seq, const int& len, std::vector<int>& data, const double& k, const double& poisson, bool verbose);
RcppExport SEXP _mlemur_optim_m(SEXP current_mSEXP, SEXP lower_mSEXP, SEXP upper_mSEXP, SEXP seqSEXP, SEXP lenSEXP, SEXP dataSEXP, SEXP kSEXP, SEXP poissonSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type current_m(current_mSEXP);
    Rcpp::traits::input_parameter< double >::type lower_m(lower_mSEXP);
    Rcpp::traits::input_parameter< double >::type upper_m(upper_mSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const int& >::type len(lenSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type poisson(poissonSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_m(current_m, lower_m, upper_m, seq, len, data, k, poisson, verbose));
    return rcpp_result_gen;
END_RCPP
}
// root_m
double root_m(double current_m, double lower_m, double upper_m, std::vector<double>& seq, const int& len, std::vector<int>& data, const double& k, const double& poisson, const double lalpha, bool verbose);
RcppExport SEXP _mlemur_root_m(SEXP current_mSEXP, SEXP lower_mSEXP, SEXP upper_mSEXP, SEXP seqSEXP, SEXP lenSEXP, SEXP dataSEXP, SEXP kSEXP, SEXP poissonSEXP, SEXP lalphaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type current_m(current_mSEXP);
    Rcpp::traits::input_parameter< double >::type lower_m(lower_mSEXP);
    Rcpp::traits::input_parameter< double >::type upper_m(upper_mSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const int& >::type len(lenSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const double& >::type poisson(poissonSEXP);
    Rcpp::traits::input_parameter< const double >::type lalpha(lalphaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(root_m(current_m, lower_m, upper_m, seq, len, data, k, poisson, lalpha, verbose));
    return rcpp_result_gen;
END_RCPP
}
// combo_optim_m
std::vector<double> combo_optim_m(double current_m, double lower_m, double upper_m, const double& R, std::vector<double>& seq1, std::vector<double>& seq2, const int& len1, const int& len2, std::vector<int>& data1, std::vector<int>& data2, const double& k1, const double& k2, const double& poisson1, const double& poisson2, bool verbose);
RcppExport SEXP _mlemur_combo_optim_m(SEXP current_mSEXP, SEXP lower_mSEXP, SEXP upper_mSEXP, SEXP RSEXP, SEXP seq1SEXP, SEXP seq2SEXP, SEXP len1SEXP, SEXP len2SEXP, SEXP data1SEXP, SEXP data2SEXP, SEXP k1SEXP, SEXP k2SEXP, SEXP poisson1SEXP, SEXP poisson2SEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type current_m(current_mSEXP);
    Rcpp::traits::input_parameter< double >::type lower_m(lower_mSEXP);
    Rcpp::traits::input_parameter< double >::type upper_m(upper_mSEXP);
    Rcpp::traits::input_parameter< const double& >::type R(RSEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq1(seq1SEXP);
    Rcpp::traits::input_parameter< std::vector<double>& >::type seq2(seq2SEXP);
    Rcpp::traits::input_parameter< const int& >::type len1(len1SEXP);
    Rcpp::traits::input_parameter< const int& >::type len2(len2SEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data1(data1SEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data2(data2SEXP);
    Rcpp::traits::input_parameter< const double& >::type k1(k1SEXP);
    Rcpp::traits::input_parameter< const double& >::type k2(k2SEXP);
    Rcpp::traits::input_parameter< const double& >::type poisson1(poisson1SEXP);
    Rcpp::traits::input_parameter< const double& >::type poisson2(poisson2SEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(combo_optim_m(current_m, lower_m, upper_m, R, seq1, seq2, len1, len2, data1, data2, k1, k2, poisson1, poisson2, verbose));
    return rcpp_result_gen;
END_RCPP
}
// optim_m_param_joint
std::vector<double> optim_m_param_joint(double current_m, double lower_m, double upper_m, double current_param, double lower_param, double upper_param, const double e, std::vector<int>& data, const int len, const std::string& option, bool verbose);
RcppExport SEXP _mlemur_optim_m_param_joint(SEXP current_mSEXP, SEXP lower_mSEXP, SEXP upper_mSEXP, SEXP current_paramSEXP, SEXP lower_paramSEXP, SEXP upper_paramSEXP, SEXP eSEXP, SEXP dataSEXP, SEXP lenSEXP, SEXP optionSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type current_m(current_mSEXP);
    Rcpp::traits::input_parameter< double >::type lower_m(lower_mSEXP);
    Rcpp::traits::input_parameter< double >::type upper_m(upper_mSEXP);
    Rcpp::traits::input_parameter< double >::type current_param(current_paramSEXP);
    Rcpp::traits::input_parameter< double >::type lower_param(lower_paramSEXP);
    Rcpp::traits::input_parameter< double >::type upper_param(upper_paramSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type option(optionSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_m_param_joint(current_m, lower_m, upper_m, current_param, lower_param, upper_param, e, data, len, option, verbose));
    return rcpp_result_gen;
END_RCPP
}
// optim_m_param_joint_2
std::vector<double> optim_m_param_joint_2(double current_m, double current_param, const double e, std::vector<int>& data, const int len, const std::string& option, bool verbose);
RcppExport SEXP _mlemur_optim_m_param_joint_2(SEXP current_mSEXP, SEXP current_paramSEXP, SEXP eSEXP, SEXP dataSEXP, SEXP lenSEXP, SEXP optionSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type current_m(current_mSEXP);
    Rcpp::traits::input_parameter< double >::type current_param(current_paramSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type option(optionSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(optim_m_param_joint_2(current_m, current_param, e, data, len, option, verbose));
    return rcpp_result_gen;
END_RCPP
}
// root_m_param_joint
std::vector<double> root_m_param_joint(double current_m, double lower_m, double upper_m, double current_param, double lower_param, double upper_param, const double e, std::vector<int>& data, const int len, const double lalpha, const int& rootparam, const int& lowerroot, const std::string& option, bool verbose);
RcppExport SEXP _mlemur_root_m_param_joint(SEXP current_mSEXP, SEXP lower_mSEXP, SEXP upper_mSEXP, SEXP current_paramSEXP, SEXP lower_paramSEXP, SEXP upper_paramSEXP, SEXP eSEXP, SEXP dataSEXP, SEXP lenSEXP, SEXP lalphaSEXP, SEXP rootparamSEXP, SEXP lowerrootSEXP, SEXP optionSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type current_m(current_mSEXP);
    Rcpp::traits::input_parameter< double >::type lower_m(lower_mSEXP);
    Rcpp::traits::input_parameter< double >::type upper_m(upper_mSEXP);
    Rcpp::traits::input_parameter< double >::type current_param(current_paramSEXP);
    Rcpp::traits::input_parameter< double >::type lower_param(lower_paramSEXP);
    Rcpp::traits::input_parameter< double >::type upper_param(upper_paramSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< const double >::type lalpha(lalphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type rootparam(rootparamSEXP);
    Rcpp::traits::input_parameter< const int& >::type lowerroot(lowerrootSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type option(optionSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(root_m_param_joint(current_m, lower_m, upper_m, current_param, lower_param, upper_param, e, data, len, lalpha, rootparam, lowerroot, option, verbose));
    return rcpp_result_gen;
END_RCPP
}
// lower_root_star_m_param_joint
std::vector<double> lower_root_star_m_param_joint(std::vector<double> max, const double e, std::vector<int>& data, const int len, double alpha, const int& rootparam, const std::string& option, bool verbose);
RcppExport SEXP _mlemur_lower_root_star_m_param_joint(SEXP maxSEXP, SEXP eSEXP, SEXP dataSEXP, SEXP lenSEXP, SEXP alphaSEXP, SEXP rootparamSEXP, SEXP optionSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type max(maxSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type rootparam(rootparamSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type option(optionSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(lower_root_star_m_param_joint(max, e, data, len, alpha, rootparam, option, verbose));
    return rcpp_result_gen;
END_RCPP
}
// fisher_info
std::vector<double> fisher_info(const double sqrtm, const double sqrtparam, const double e, const int len, const std::string& option);
RcppExport SEXP _mlemur_fisher_info(SEXP sqrtmSEXP, SEXP sqrtparamSEXP, SEXP eSEXP, SEXP lenSEXP, SEXP optionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type sqrtm(sqrtmSEXP);
    Rcpp::traits::input_parameter< const double >::type sqrtparam(sqrtparamSEXP);
    Rcpp::traits::input_parameter< const double >::type e(eSEXP);
    Rcpp::traits::input_parameter< const int >::type len(lenSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type option(optionSEXP);
    rcpp_result_gen = Rcpp::wrap(fisher_info(sqrtm, sqrtparam, e, len, option));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mlemur_rluria_list", (DL_FUNC) &_mlemur_rluria_list, 12},
    {"_mlemur_rluria_vec", (DL_FUNC) &_mlemur_rluria_vec, 12},
    {"_mlemur_aux_seq_integrate", (DL_FUNC) &_mlemur_aux_seq_integrate, 6},
    {"_mlemur_aux_seq", (DL_FUNC) &_mlemur_aux_seq, 4},
    {"_mlemur_aux_seq_deriv1_deriv2", (DL_FUNC) &_mlemur_aux_seq_deriv1_deriv2, 5},
    {"_mlemur_aux_seq_lag_ext", (DL_FUNC) &_mlemur_aux_seq_lag_ext, 5},
    {"_mlemur_aux_seq_death_ext", (DL_FUNC) &_mlemur_aux_seq_death_ext, 5},
    {"_mlemur_prob_ld", (DL_FUNC) &_mlemur_prob_ld, 3},
    {"_mlemur_xi_seq", (DL_FUNC) &_mlemur_xi_seq, 3},
    {"_mlemur_prob_b0", (DL_FUNC) &_mlemur_prob_b0, 5},
    {"_mlemur_prob_mutations", (DL_FUNC) &_mlemur_prob_mutations, 9},
    {"_mlemur_calc_probs", (DL_FUNC) &_mlemur_calc_probs, 6},
    {"_mlemur_calc_probs_2", (DL_FUNC) &_mlemur_calc_probs_2, 5},
    {"_mlemur_logprob", (DL_FUNC) &_mlemur_logprob, 6},
    {"_mlemur_logprob_boost", (DL_FUNC) &_mlemur_logprob_boost, 6},
    {"_mlemur_optim_m", (DL_FUNC) &_mlemur_optim_m, 9},
    {"_mlemur_root_m", (DL_FUNC) &_mlemur_root_m, 10},
    {"_mlemur_combo_optim_m", (DL_FUNC) &_mlemur_combo_optim_m, 15},
    {"_mlemur_optim_m_param_joint", (DL_FUNC) &_mlemur_optim_m_param_joint, 11},
    {"_mlemur_optim_m_param_joint_2", (DL_FUNC) &_mlemur_optim_m_param_joint_2, 7},
    {"_mlemur_root_m_param_joint", (DL_FUNC) &_mlemur_root_m_param_joint, 14},
    {"_mlemur_lower_root_star_m_param_joint", (DL_FUNC) &_mlemur_lower_root_star_m_param_joint, 8},
    {"_mlemur_fisher_info", (DL_FUNC) &_mlemur_fisher_info, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_mlemur(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
