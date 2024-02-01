// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// setParallelCRT
void setParallelCRT(SEXP parallel_, int cores_);
RcppExport SEXP _crctStepdown_setParallelCRT(SEXP parallel_SEXP, SEXP cores_SEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type parallel_(parallel_SEXP);
    Rcpp::traits::input_parameter< int >::type cores_(cores_SEXP);
    setParallelCRT(parallel_, cores_);
    return R_NilValue;
END_RCPP
}
// qscore_impl
double qscore_impl(const Eigen::VectorXd& resids, Eigen::VectorXd tr, const Eigen::VectorXd& xb, const Eigen::MatrixXd& invS, const std::string& family2, const Eigen::ArrayXXd& Z, bool weight);
RcppExport SEXP _crctStepdown_qscore_impl(SEXP residsSEXP, SEXP trSEXP, SEXP xbSEXP, SEXP invSSEXP, SEXP family2SEXP, SEXP ZSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type tr(trSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family2(family2SEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(qscore_impl(resids, tr, xb, invS, family2, Z, weight));
    return rcpp_result_gen;
END_RCPP
}
// permutation_test_impl
Eigen::VectorXd permutation_test_impl(const Eigen::VectorXd& resids, const Eigen::MatrixXd& tr_mat, const Eigen::VectorXd& xb, const Eigen::MatrixXd& invS, const std::string& family2, const Eigen::ArrayXXd& Z, bool weight, int iter, bool verbose);
RcppExport SEXP _crctStepdown_permutation_test_impl(SEXP residsSEXP, SEXP tr_matSEXP, SEXP xbSEXP, SEXP invSSEXP, SEXP family2SEXP, SEXP ZSEXP, SEXP weightSEXP, SEXP iterSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type resids(residsSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type tr_mat(tr_matSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type xb(xbSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type family2(family2SEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(permutation_test_impl(resids, tr_mat, xb, invS, family2, Z, weight, iter, verbose));
    return rcpp_result_gen;
END_RCPP
}
// confint_search
Rcpp::List confint_search(Eigen::VectorXd start, Eigen::VectorXd b, int n, int nmodel, const Rcpp::List& Xnull_, const Rcpp::List& y, const Rcpp::NumericVector& tr_, const Eigen::MatrixXd& new_tr_mat, const Rcpp::List& invS, Rcpp::List family, Rcpp::List family2, const Eigen::ArrayXXd& Z, const Rcpp::String& type, int nsteps, bool weight, double alpha, bool verbose);
RcppExport SEXP _crctStepdown_confint_search(SEXP startSEXP, SEXP bSEXP, SEXP nSEXP, SEXP nmodelSEXP, SEXP Xnull_SEXP, SEXP ySEXP, SEXP tr_SEXP, SEXP new_tr_matSEXP, SEXP invSSEXP, SEXP familySEXP, SEXP family2SEXP, SEXP ZSEXP, SEXP typeSEXP, SEXP nstepsSEXP, SEXP weightSEXP, SEXP alphaSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type start(startSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type nmodel(nmodelSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Xnull_(Xnull_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type tr_(tr_SEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type new_tr_mat(new_tr_matSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type invS(invSSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type family(familySEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type family2(family2SEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXXd& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< const Rcpp::String& >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type nsteps(nstepsSEXP);
    Rcpp::traits::input_parameter< bool >::type weight(weightSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(confint_search(start, b, n, nmodel, Xnull_, y, tr_, new_tr_mat, invS, family, family2, Z, type, nsteps, weight, alpha, verbose));
    return rcpp_result_gen;
END_RCPP
}
// simpleLM
Rcpp::List simpleLM(SEXP y_, SEXP X_);
RcppExport SEXP _crctStepdown_simpleLM(SEXP y_SEXP, SEXP X_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type y_(y_SEXP);
    Rcpp::traits::input_parameter< SEXP >::type X_(X_SEXP);
    rcpp_result_gen = Rcpp::wrap(simpleLM(y_, X_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_crctStepdown_setParallelCRT", (DL_FUNC) &_crctStepdown_setParallelCRT, 2},
    {"_crctStepdown_qscore_impl", (DL_FUNC) &_crctStepdown_qscore_impl, 7},
    {"_crctStepdown_permutation_test_impl", (DL_FUNC) &_crctStepdown_permutation_test_impl, 9},
    {"_crctStepdown_confint_search", (DL_FUNC) &_crctStepdown_confint_search, 17},
    {"_crctStepdown_simpleLM", (DL_FUNC) &_crctStepdown_simpleLM, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_crctStepdown(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
