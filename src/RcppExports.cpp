// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// vec2mat
Eigen::MatrixXd vec2mat(Eigen::VectorXd x, int type, int num);
RcppExport SEXP _eimpute_vec2mat(SEXP xSEXP, SEXP typeSEXP, SEXP numSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type num(numSEXP);
    rcpp_result_gen = Rcpp::wrap(vec2mat(x, type, num));
    return rcpp_result_gen;
END_RCPP
}
// biscale_alt
List biscale_alt(Eigen::MatrixXd x, Eigen::MatrixXi ind, Eigen::VectorXd obsrow, Eigen::VectorXd obscol, int max_it, double tol, Eigen::VectorXd alpha, Eigen::VectorXd beta, Eigen::VectorXd tau, Eigen::VectorXd gamma, bool row_mean, bool col_mean, bool row_std, bool col_std);
RcppExport SEXP _eimpute_biscale_alt(SEXP xSEXP, SEXP indSEXP, SEXP obsrowSEXP, SEXP obscolSEXP, SEXP max_itSEXP, SEXP tolSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP tauSEXP, SEXP gammaSEXP, SEXP row_meanSEXP, SEXP col_meanSEXP, SEXP row_stdSEXP, SEXP col_stdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXi >::type ind(indSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type obsrow(obsrowSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type obscol(obscolSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< bool >::type row_mean(row_meanSEXP);
    Rcpp::traits::input_parameter< bool >::type col_mean(col_meanSEXP);
    Rcpp::traits::input_parameter< bool >::type row_std(row_stdSEXP);
    Rcpp::traits::input_parameter< bool >::type col_std(col_stdSEXP);
    rcpp_result_gen = Rcpp::wrap(biscale_alt(x, ind, obsrow, obscol, max_it, tol, alpha, beta, tau, gamma, row_mean, col_mean, row_std, col_std));
    return rcpp_result_gen;
END_RCPP
}
// trun_svd
Eigen::MatrixXd trun_svd(Eigen::MatrixXd X, int k);
RcppExport SEXP _eimpute_trun_svd(SEXP XSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(trun_svd(X, k));
    return rcpp_result_gen;
END_RCPP
}
// kkt_fix
List kkt_fix(Eigen::MatrixXi& omega, Eigen::MatrixXd& noise, Eigen::VectorXd& X, int m, int n, int rank, int max_it, double tol, int type, bool init);
RcppExport SEXP _eimpute_kkt_fix(SEXP omegaSEXP, SEXP noiseSEXP, SEXP XSEXP, SEXP mSEXP, SEXP nSEXP, SEXP rankSEXP, SEXP max_itSEXP, SEXP tolSEXP, SEXP typeSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXi& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type noise(noiseSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type rank(rankSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(kkt_fix(omega, noise, X, m, n, rank, max_it, tol, type, init));
    return rcpp_result_gen;
END_RCPP
}
// cv_rank
Eigen::VectorXd cv_rank(Eigen::MatrixXi& omega, Eigen::MatrixXd& noise, Eigen::VectorXd& X, int m, int n, int r_min, int r_max, int n_fold, int max_it, double tol, int type, bool init);
RcppExport SEXP _eimpute_cv_rank(SEXP omegaSEXP, SEXP noiseSEXP, SEXP XSEXP, SEXP mSEXP, SEXP nSEXP, SEXP r_minSEXP, SEXP r_maxSEXP, SEXP n_foldSEXP, SEXP max_itSEXP, SEXP tolSEXP, SEXP typeSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXi& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type noise(noiseSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type r_min(r_minSEXP);
    Rcpp::traits::input_parameter< int >::type r_max(r_maxSEXP);
    Rcpp::traits::input_parameter< int >::type n_fold(n_foldSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(cv_rank(omega, noise, X, m, n, r_min, r_max, n_fold, max_it, tol, type, init));
    return rcpp_result_gen;
END_RCPP
}
// ic_rank
Eigen::VectorXd ic_rank(Eigen::MatrixXi& omega, Eigen::MatrixXd& noise, Eigen::VectorXd& X, int m, int n, int r_min, int r_max, int max_it, double tol, int type, bool init);
RcppExport SEXP _eimpute_ic_rank(SEXP omegaSEXP, SEXP noiseSEXP, SEXP XSEXP, SEXP mSEXP, SEXP nSEXP, SEXP r_minSEXP, SEXP r_maxSEXP, SEXP max_itSEXP, SEXP tolSEXP, SEXP typeSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXi& >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd& >::type noise(noiseSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type r_min(r_minSEXP);
    Rcpp::traits::input_parameter< int >::type r_max(r_maxSEXP);
    Rcpp::traits::input_parameter< int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(ic_rank(omega, noise, X, m, n, r_min, r_max, max_it, tol, type, init));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eimpute_vec2mat", (DL_FUNC) &_eimpute_vec2mat, 3},
    {"_eimpute_biscale_alt", (DL_FUNC) &_eimpute_biscale_alt, 14},
    {"_eimpute_trun_svd", (DL_FUNC) &_eimpute_trun_svd, 2},
    {"_eimpute_kkt_fix", (DL_FUNC) &_eimpute_kkt_fix, 10},
    {"_eimpute_cv_rank", (DL_FUNC) &_eimpute_cv_rank, 12},
    {"_eimpute_ic_rank", (DL_FUNC) &_eimpute_ic_rank, 11},
    {NULL, NULL, 0}
};

RcppExport void R_init_eimpute(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
