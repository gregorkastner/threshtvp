// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// get_threshold
NumericVector get_threshold(NumericMatrix Achg1, NumericMatrix SIGMAS1, NumericMatrix SIGMAS2, NumericMatrix threshgrid, NumericMatrix Aapprox);
RcppExport SEXP _threshtvp_get_threshold(SEXP Achg1SEXP, SEXP SIGMAS1SEXP, SEXP SIGMAS2SEXP, SEXP threshgridSEXP, SEXP AapproxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Achg1(Achg1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SIGMAS1(SIGMAS1SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type SIGMAS2(SIGMAS2SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type threshgrid(threshgridSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Aapprox(AapproxSEXP);
    rcpp_result_gen = Rcpp::wrap(get_threshold(Achg1, SIGMAS1, SIGMAS2, threshgrid, Aapprox));
    return rcpp_result_gen;
END_RCPP
}
// dinvgamma
double dinvgamma(const double x, const double a, const double b);
RcppExport SEXP _threshtvp_dinvgamma(SEXP xSEXP, SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type x(xSEXP);
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(dinvgamma(x, a, b));
    return rcpp_result_gen;
END_RCPP
}
// KF
List KF(NumericMatrix y, NumericMatrix Z, NumericMatrix Ht, NumericMatrix Qtt, int m, int p, int t, NumericVector B0, NumericMatrix V0);
RcppExport SEXP _threshtvp_KF(SEXP ySEXP, SEXP ZSEXP, SEXP HtSEXP, SEXP QttSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP B0SEXP, SEXP V0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Qtt(QttSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V0(V0SEXP);
    rcpp_result_gen = Rcpp::wrap(KF(y, Z, Ht, Qtt, m, p, t, B0, V0));
    return rcpp_result_gen;
END_RCPP
}
// get_lik
double get_lik(NumericMatrix y, NumericMatrix Z, NumericMatrix Ht, NumericMatrix Qtt, int m, int p, int t, NumericVector B0, NumericMatrix V0);
RcppExport SEXP _threshtvp_get_lik(SEXP ySEXP, SEXP ZSEXP, SEXP HtSEXP, SEXP QttSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP B0SEXP, SEXP V0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Qtt(QttSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V0(V0SEXP);
    rcpp_result_gen = Rcpp::wrap(get_lik(y, Z, Ht, Qtt, m, p, t, B0, V0));
    return rcpp_result_gen;
END_RCPP
}
// KF_fast
List KF_fast(NumericMatrix y, NumericMatrix Z, NumericMatrix Ht, NumericMatrix Qtt, int m, int p, int t, NumericVector B0, NumericMatrix V0);
RcppExport SEXP _threshtvp_KF_fast(SEXP ySEXP, SEXP ZSEXP, SEXP HtSEXP, SEXP QttSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP B0SEXP, SEXP V0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Qtt(QttSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V0(V0SEXP);
    rcpp_result_gen = Rcpp::wrap(KF_fast(y, Z, Ht, Qtt, m, p, t, B0, V0));
    return rcpp_result_gen;
END_RCPP
}
// KF_MH
List KF_MH(NumericMatrix y, NumericMatrix Z, NumericMatrix Ht, NumericMatrix Qtt, int m, int p, int t, NumericVector B0, NumericMatrix V0);
RcppExport SEXP _threshtvp_KF_MH(SEXP ySEXP, SEXP ZSEXP, SEXP HtSEXP, SEXP QttSEXP, SEXP mSEXP, SEXP pSEXP, SEXP tSEXP, SEXP B0SEXP, SEXP V0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Ht(HtSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Qtt(QttSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type B0(B0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type V0(V0SEXP);
    rcpp_result_gen = Rcpp::wrap(KF_MH(y, Z, Ht, Qtt, m, p, t, B0, V0));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_threshtvp_get_threshold", (DL_FUNC) &_threshtvp_get_threshold, 5},
    {"_threshtvp_dinvgamma", (DL_FUNC) &_threshtvp_dinvgamma, 3},
    {"_threshtvp_KF", (DL_FUNC) &_threshtvp_KF, 9},
    {"_threshtvp_get_lik", (DL_FUNC) &_threshtvp_get_lik, 9},
    {"_threshtvp_KF_fast", (DL_FUNC) &_threshtvp_KF_fast, 9},
    {"_threshtvp_KF_MH", (DL_FUNC) &_threshtvp_KF_MH, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_threshtvp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}