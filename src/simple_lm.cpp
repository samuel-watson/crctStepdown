#include "eigenext.h"
#include "openmpheader.h"

using namespace Rcpp;

//' A very basic linear model solver
//'
//' Returns the OLS paramter estimates and fitted values. Used internally for quick
//' fitting of null models.
//'
//' @param y_ A vector of outcome values
//' @param X_ The design matrix of fixed effects
//' @return A list with the parameter values and fitted values
// [[Rcpp::export]]
Rcpp::List simpleLM(SEXP y_,
                    SEXP X_){
  Eigen::MatrixXd X = as<Eigen::MatrixXd>(X_);
  Eigen::VectorXd y = as<Eigen::VectorXd>(y_);
  Eigen::VectorXd XtX = X.transpose() * X;
  XtX = XtX.llt().solve(Eigen::MatrixXd::Identity(X.cols(),X.cols()));
  Eigen::VectorXd beta = XtX * X.transpose() * y;
  Eigen::VectorXd fitted = X * beta;

  return List::create(_["fitted.values"] = wrap(fitted),
                      _["b"] = wrap(beta));
}
