// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// powThresh
NumericVector powThresh(NumericVector z, double lambda, double q);
RcppExport SEXP _powopt_powThresh(SEXP zSEXP, SEXP lambdaSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(powThresh(z, lambda, q));
    return rcpp_result_gen;
END_RCPP
}
// powObj
double powObj(NumericVector beta, NumericMatrix Q, NumericVector l, double yty, double sigmasq, double lambda, double q);
RcppExport SEXP _powopt_powObj(SEXP betaSEXP, SEXP QSEXP, SEXP lSEXP, SEXP ytySEXP, SEXP sigmasqSEXP, SEXP lambdaSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Q(QSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type l(lSEXP);
    Rcpp::traits::input_parameter< double >::type yty(ytySEXP);
    Rcpp::traits::input_parameter< double >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(powObj(beta, Q, l, yty, sigmasq, lambda, q));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_powopt_powThresh", (DL_FUNC) &_powopt_powThresh, 3},
    {"_powopt_powObj", (DL_FUNC) &_powopt_powObj, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_powopt(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
