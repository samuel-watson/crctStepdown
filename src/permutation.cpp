#define STRICT_R_HEADERS
#include <string>

#include "eigenext.h"
#include "glm.h"
#include "openmpheader.h"

using namespace Rcpp;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

Rcpp::List fast_glm_impl(Rcpp::NumericMatrix Xs,
                   Rcpp::NumericVector ys,
                   Rcpp::NumericVector weightss,
                   Rcpp::NumericVector offsets,
                   Rcpp::NumericVector starts,
                   Rcpp::NumericVector mus,
                   Rcpp::NumericVector etas,
                   Function var,
                   Function mu_eta,
                   Function linkinv,
                   Function dev_resids,
                   Function valideta,
                   Function validmu,
                   int type,
                   double tol,
                   int maxit)
{
  const Map<MatrixXd>  X(as<Map<MatrixXd> >(Xs));
  const Map<VectorXd>  y(as<Map<VectorXd> >(ys));
  const Map<VectorXd>  weights(as<Map<VectorXd> >(weightss));
  const Map<VectorXd>  offset(as<Map<VectorXd> >(offsets));
  const Map<VectorXd>  beta_init(as<Map<VectorXd> >(starts));
  const Map<VectorXd>  mu_init(as<Map<VectorXd> >(mus));
  const Map<VectorXd>  eta_init(as<Map<VectorXd> >(etas));
  Index                n = X.rows();
  if ((Index)y.size() != n) throw invalid_argument("size mismatch");

  // instantiate fitting class
  GlmBase<Eigen::VectorXd, Eigen::MatrixXd> *glm_solver = NULL;

  bool is_big_matrix = false;

  glm_solver = new glm(X, y, weights, offset,
                       var, mu_eta, linkinv, dev_resids,
                       valideta, validmu, tol, maxit, type,
                       is_big_matrix);

  // initialize parameters
  glm_solver->init_parms(beta_init, mu_init, eta_init);
  int iters = glm_solver->solve(maxit);
  VectorXd eta       = glm_solver->get_eta();

  delete glm_solver;
  return List::create(_["linear.predictors"] = eta);
}

double gaussian_cdf(double x){
  return R::pnorm(x, 0, 1, true, false);
}

vec gaussian_cdf_vec(const vec& v){
  vec res(v.size());
  for (int i = 0; i < v.size(); ++i)
    res(i) = gaussian_cdf(v(i));
  return res;
}

double gaussian_pdf(double x){
  return R::dnorm(x, 0, 1, false);
}

vec gaussian_pdf_vec(const vec& v){
  vec res(v.size());
  for (int i = 0; i < v.size(); ++i)
    res(i) = gaussian_pdf(v(i));
  return res;
}

inline vec get_G(const vec &x,
                 const Rcpp::String &family2){
  vec dx(x.size());
  if(family2 == "identity"){
    dx.setConstant(1);
  } else if(family2 == "log"){
    dx = x;
    dx = dx.array().exp().matrix();
  } else if(family2 == "logit"){
    Eigen::ArrayXd expdx = x.array().exp();
    Eigen::ArrayXd expdx2 = expdx * (1-expdx);
    expdx = (1+expdx)*(1+expdx);
    dx = (expdx2 * expdx.inverse()).matrix();
  } else if(family2 == "probit"){
    dx = (gaussian_pdf_vec(x)).array().inverse().matrix();
    dx *= -1;
  }
  return dx;
}



//' The quasi-score statistic for a generalised linear mixed model
//'
//' Generates the quasi-score statistic for a generalised linear mixed model
//'
//' @param resids A numeric vector of generalised residuals
//' @param tr A numeric vector of 1s (treatment group) and -1s (control group)
//' @param xb A numeric vector of fitted linear predictors
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family2 A string naming the link function
//' @param Z A matrix with columns indicating cluster membership
//' @param weight Logical value indicating whether to use the weighted statistic (TRUE) or the unweighted statistic (FALSE)
//' @return A scalar value with the value of the statistic
// [[Rcpp::export]]
double qscore_impl(const Eigen::VectorXd &resids,
                   Eigen::VectorXd tr,
                   const Eigen::VectorXd &xb,
                   const Eigen::MatrixXd &invS,
                   const std::string &family2,
                   const Eigen::ArrayXXd &Z,
                   bool weight=true) {
  vec g(xb.size());
  vec q(Z.cols());
  if (weight){
    g = get_G(xb, family2);
//#pragma omp parallel for
    for(int j=0; j<Z.cols(); j++){
      Eigen::ArrayXd zcol = Z.col(j);
      Eigen::ArrayXi idx = Eigen_ext::find<double>(zcol,1);
      mat invScl = Eigen_ext::mat_indexing(invS,idx,idx);
      vec gcl(idx.size());
      for(int k = 0; k<idx.size(); k++){
        gcl(k) = g(idx(k));
      }
      vec gS = invScl * gcl;
      for(int k = 0; k<idx.size(); k++){
        gcl(k) = tr(idx(k))*resids(idx(k));
      }
      q(j) = gS.transpose()*gcl;
    }
  } else {
    vec tres = (tr.array()*resids.array()).matrix();
    q = Z.matrix().transpose() * tres;
  }
  double denom = q.transpose()*q;
  double numer = q.sum();
  return std::abs(numer/pow(denom,0.5));
}

//' Generates realisations of the permutational test statistic distribution
//'
//' Generates realisations of the permutational test statistic distribution from a given matrix of permutations
//'
//' @param resids A numeric vector of generalised residuals
//' @param tr_mat A matrix. Each column is a new random treatment allocation with 1s (treatment group) and 0s (control group)
//' @param xb A numeric vector of fitted linear predictors
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family2 A string naming the link function
//' @param Z A matrix with columns indicating cluster membership
//' @param weight Logical value indicating whether to use the weighted statistic (TRUE) or the unweighted statistic (FALSE)
//' @param iter Integer. Number of permutation test iterations.
//' @param verbose Logical indicating whether to report detailed output
//' @return A numeric vector of quasi-score test statistics for each of the permutations
// [[Rcpp::export]]
Eigen::VectorXd permutation_test_impl(const Eigen::VectorXd &resids,
                                const Eigen::MatrixXd &tr_mat,
                                const Eigen::VectorXd &xb,
                                const Eigen::MatrixXd &invS,
                                const std::string &family2,
                                const Eigen::ArrayXXd &Z,
                                bool weight,
                                int iter = 1000,
                                bool verbose = true) {
  if (verbose) Rcpp::Rcout << "Starting permutations\n" << std::endl;

  vec qtest(iter);
  qtest.setZero();
#pragma omp parallel for
  for (int i = 0; i < iter; ++i) {
    vec tr = tr_mat.col(i);
    Eigen_ext::replaceVec(tr,0,-1);
    qtest(i) = qscore_impl(resids, tr, xb, invS, family2, Z, weight);
  }
  return qtest;
}

//' Confidence interval search procedure
//'
//' Search for the bound of a confidence interval using permutation test statistics
//'
//' @param start Numeric value indicating the starting value for the search procedure
//' @param b Numeric value indicating the parameter estimate
//' @param n Integer indicating the sample size
//' @param nmodel Integer. The number of models
//' @param Xnull_ Numeric matrix. The covariate design matrix with the treatment variable removed
//' @param y Numeric vector of response variables
//' @param tr_ Numeric vector. The original random allocation (0s and 1s)
//' @param new_tr_mat A matrix. Each column is a new random treatment allocation with 1s (treatment group) and 0s (control group)
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family A \link{stats}[family] object
//' @param family2 A string naming the link function
//' @param Z Matrix. Random effects design matrix describing cluster membership
//' @param type String. Either "rw" for Romano-Wolf, "b" or "br" for bonferroni, "h" or "hr" for Holm, or "none"
//' @param nsteps Integer specifying the number of steps of the search procedure
//' @param weight Logical indicating whether to use the weighted (TRUE) or unweighted (FALSE) test statistic
//' @param alpha The function generates (1-alpha)*100% confidence intervals. Default is 0.05.
//' @param verbose Logical indicating whether to provide detailed output.
//' @return The estimated confidence interval bound
// [[Rcpp::export]]
Rcpp::List confint_search(Eigen::VectorXd start,
                      Eigen::VectorXd b,
                      int n,
                      int nmodel,
                      const Rcpp::List &Xnull_,
                      const Rcpp::List &y,
                      const Rcpp::NumericVector &tr_,
                      const Eigen::MatrixXd &new_tr_mat,
                      const Rcpp::List &invS,
                      Rcpp::List family,
                      Rcpp::List family2,
                      const Eigen::ArrayXXd &Z,
                      const Rcpp::String &type,
                      int nsteps = 1000,
                      bool weight = true,
                      double alpha = 0.05,
                      bool verbose = true) {

  vec tr = as<vec>(tr_);
  Eigen_ext::replaceVec(tr,0,-1);

  Rcpp::NumericVector weights_(n);
  weights_.fill(1);

  vec bound = start;
  mat boundvals(nsteps,nmodel);
  vec qstat(nmodel);
  qstat.setZero();
  vec qtest(nmodel);
  vec step(nmodel);

  for (int i = 1; i <= nsteps; ++i) {
    boundvals.row(i-1) = bound.transpose();
    vec new_tr = new_tr_mat.col(i-1);
    Eigen_ext::replaceVec(new_tr,0,-1);
    for(int j = 0; j < nmodel; j++){
      vec tr2 = as<vec>(tr_);
      Eigen_ext::replaceVec(tr2,1,bound(j));

      //tr2.replace(1,bound(j));
      Rcpp::List faml = family[j];
      Rcpp::NumericMatrix Xs = Rcpp::as<Rcpp::NumericMatrix>(Xnull_[j]);
      Rcpp::NumericVector ys = Rcpp::as<Rcpp::NumericVector>(y[j]);
      Rcpp::NumericVector startval(Xs.ncol());
      Function linkfun = faml["linkfun"];
      Function linkinv = faml["linkinv"];
      Function varf = faml["variance"];
      Function mueta = faml["mu.eta"];
      Function devres = faml["dev.resids"];
      Function valideta = faml["valideta"];
      Function validmu = faml["validmu"];
      std::string fname = faml["family"];
      Rcpp::NumericVector eta(n);
      if(fname == "poisson"){
        eta = linkfun(ys+0.1);
      } else if(fname == "gaussian"){
        eta = ys;
      } else if(fname == "binomial"){
        eta = linkfun((ys + 0.5)/2);
      }

      Rcpp::NumericVector mu = linkinv(eta);
      Rcpp::NumericVector offs_ =  Rcpp::as<Rcpp::NumericVector>(wrap(tr2));
      Rcpp::List result = fast_glm_impl(Xs,
                                        ys,
                                        weights_,
                                        offs_,
                                        startval,
                                        mu,
                                        eta,
                                        varf,
                                        mueta,
                                        linkinv,
                                        devres,
                                        valideta,
                                        validmu,
                                        0,
                                        1e-7,
                                        100);
      const Map<VectorXd>  xb(as<Map<VectorXd> >(result["linear.predictors"]));
      Rcpp::NumericVector xbvec(wrap(xb));
      vec ypred = Rcpp::as<vec>(linkinv(xb));
      vec resids = as<vec>(y[j]) - ypred;
      qstat(j) = qscore_impl(resids, tr, xb, as<mat>(invS[j]), as<std::string>(family2[j]), Z, weight);
      qtest(j) = qscore_impl(resids, new_tr, xb, as<mat>(invS[j]), as<std::string>(family2[j]), Z, weight);
    }
    // Rcpp::Rcout << "\n qstat: " << qstat.t();
    // Rcpp::Rcout << "\n qtest: " << qtest.t();
    //Rcpp::Rcout << "\nDone fits and stats" << qstat.t() << "and " << qtest.t();
    //arma::uvec pos_t = arma::sort_index(qstat);
    Eigen::ArrayXi pos_t = Eigen_ext::sort_indexes(qstat.array());
    Eigen::ArrayXi rjct(nmodel);
    double k;
    double k_tmp;
    double stat_test;
    vec Jval(nmodel);
    vec step(nmodel);

    if(type=="rw"){
      k_tmp = gaussian_cdf(1-alpha);
      k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));
      Jval.setConstant(1);
      bool reject_rest = false;
      for(int j = 0; j < nmodel; j++){
        if(!reject_rest){
          stat_test = Eigen_ext::max_val_subvec(qtest,pos_t.head(nmodel-j));
          rjct(pos_t(nmodel - 1 - j)) = qstat(pos_t(nmodel-1-j)) >= stat_test;
        } else {
          rjct(pos_t(nmodel - 1 - j)) = 0;
        }
        step(pos_t(nmodel - 1 - j)) = k *(b(pos_t(nmodel - 1 - j)) - bound(pos_t(nmodel - 1 - j)));
        if(!rjct(pos_t(nmodel - 1 - j))){
          reject_rest = true;
        }
      }
    } else {
      if(type == "b" || type == "br"){
        k_tmp = gaussian_cdf(1-alpha/nmodel);
        k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));
        Jval.setConstant(nmodel);
        for(int j = 0; j < nmodel; j++){
          rjct(j) = qstat(j) > qtest(j);
          step(j) = k *(b(j) - bound(j));
        }
      }
      if(type == "h"|| type=="hr"){
        bool reject_rest = false;
        for(int j = 0; j < nmodel; j++){
          if(!reject_rest){
            rjct(pos_t(nmodel - 1 - j)) = qstat(pos_t(nmodel-1-j)) >= qtest(pos_t(nmodel - 1 - j));
          } else {
            rjct(pos_t(nmodel - 1 - j)) = 0;
          }
          k_tmp = gaussian_cdf(1-alpha/(nmodel-j));
          k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));
          step(pos_t(nmodel - 1 - j)) = k *(b(pos_t(nmodel - 1 - j)) - bound(pos_t(nmodel - 1 - j)));
          Jval(pos_t(nmodel - 1 - j)) = nmodel-j;
          if(!rjct(pos_t(nmodel - 1 - j))){
            reject_rest = true;
          }
        }
      }
      if(type=="none"){
        k_tmp = gaussian_cdf(1-alpha);
        k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));
        Jval.setConstant(1);
        for(int j = 0; j < nmodel; j++){
          rjct(j) = qstat(j) > qtest(j);
          step(j) = k *(b(j) - bound(j));
        }
      }
    }

    for(int j = 0; j < nmodel; j++){
      if(rjct(j)){
        bound(j) += step(j)*(alpha/Jval(j))/i;
      } else {
        bound(j) -= step(j)*(1-alpha/Jval(j))/i;
      }
    }

    if(verbose && (i % 50 == 0 || i == 1)){
      Rcpp::Rcout << "\rStep = " << i << " bound: " << bound.transpose() << std::endl;
    }
  }
  return List::create(_["bound"] = bound,
                      _["values"] = boundvals);
}
