#define STRICT_R_HEADERS
#include <string>

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "glm.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Map;

typedef MatrixXd::Index Index;

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


  // maximize likelihood
  int iters = glm_solver->solve(maxit);

  //VectorXd beta      = glm_solver->get_beta();
  //VectorXd se        = glm_solver->get_se();
  //VectorXd mu        = glm_solver->get_mu();
  VectorXd eta       = glm_solver->get_eta();
  //VectorXd wts       = glm_solver->get_w();
  //VectorXd pweights  = glm_solver->get_weights();

  //double dev         = glm_solver->get_dev();
  //int rank           = glm_solver->get_rank();
  //bool converged     = glm_solver->get_converged();

  // int df = X.rows() - rank;

  delete glm_solver;
  return List::create(_["linear.predictors"] = eta);
  // return List::create(_["coefficients"]      = beta,
  //                     _["se"]                = se,
  //                     _["fitted.values"]     = mu,
  //                     _["linear.predictors"] = eta,
  //                     _["deviance"]          = dev,
  //                     _["weights"]           = wts,
  //                     _["prior.weights"]     = pweights,
  //                     _["rank"]              = rank,
  //                     _["df.residual"]       = df,
  //                     _["iter"]              = iters,
  //                     _["converged"]         = converged);
}

double gaussian_cdf(double x){
  return R::pnorm(x, 0, 1, true, false);
}

arma::vec gaussian_cdf_vec(const arma::vec& v){
  arma::vec res = arma::zeros<arma::vec>(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i)
    res[i] = gaussian_cdf(v[i]);
  return res;
}

double gaussian_pdf(double x){
  return R::dnorm(x, 0, 1, false);
}

arma::vec gaussian_pdf_vec(const arma::vec& v){
  arma::vec res = arma::zeros<arma::vec>(v.n_elem);
  for (arma::uword i = 0; i < v.n_elem; ++i)
    res[i] = gaussian_pdf(v[i]);
  return res;
}

inline arma::vec get_G(const arma::vec &x,
                       const Rcpp::String &family2){
  arma::vec dx;
  if(family2 == "identity"){
    dx = arma::ones<arma::vec>(x.n_elem);
  } else if(family2 == "log"){
    dx = exp(x);
  } else if(family2 == "logit"){
    dx = (exp(x)/(1+exp(x))) % (1-exp(x)/(1+exp(x)));
  } else if(family2 == "probit"){
    dx = -1/gaussian_pdf_vec(x);
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
double qscore_impl(const arma::vec &resids,
                   arma::vec tr,
                   const arma::vec &xb,
                   const arma::mat &invS,
                   const std::string &family2,
                   const arma::mat &Z,
                   bool weight=true) {
  arma::rowvec g;
  arma::vec q(Z.n_cols);
  if (weight){
    g = arma::trans(get_G(xb, family2));
#pragma omp parallel for
    for(arma::uword j=0; j<Z.n_cols; j++){
      arma::uvec idx = arma::find(Z.col(j)==1);
      arma::rowvec gS = g(idx) * invS(idx,idx);
      q(j) = arma::dot((Z.col(j) % arma::trans(gS)),(tr(idx) % resids(idx)));
    }
    // arma::rowvec gS = g * invS;
    // arma::mat Z_ = Z;
    // for(arma::uword i=0; i < Z.n_cols; i++){
    //   Z_.col(i) = Z.col(i) % arma::trans(gS);
    // }
    // //Zt.each_row() % g;
    // q = Z_.t() * (tr % resids);
    //q = arma::as_scalar((g * invS) * (tr % resids));
  } else {
    q = Z.t() * (tr % resids);//arma::as_scalar(arma::dot(tr, resids));
  }
  double denom = arma::as_scalar(arma::dot(q.t(), q));
  double numer = arma::sum(q);
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
//' @param verbose Logical indicating whether to report detailed output
//' @return A numeric vector of quasi-score test statistics for each of the permutations
// [[Rcpp::export]]
arma::vec permutation_test_impl(const arma::vec &resids,
                                const arma::mat &tr_mat,
                                const arma::vec &xb,
                                const arma::mat &invS,
                                const std::string &family2,
                                const arma::mat &Z,
                                bool weight,
                                int iter = 1000,
                                bool verbose = true) {
  if (verbose) Rcpp::Rcout << "Starting permutations\n" << std::endl;

  arma::vec qtest = arma::zeros<arma::vec>(iter);
#pragma omp parallel for //uncomment for build
  for (int i = 0; i < iter; ++i) {
    arma::vec tr = tr_mat.col(i);
    tr.replace(0, -1);
    qtest[i] = qscore_impl(resids, tr, xb, invS, family2, Z, weight);
  }
  return qtest;
}

//' Confidence interval search procedure
//'
//' Search for the bound of a confidence interval using permutation test statistics
//'
//' @param start Numeric value indicating the starting value for the search procedure
//' @param b Numeric value indicating the parameter estimate
//' @param Xnull_ Numeric matrix. The covariate design matrix with the treatment variable removed
//' @param y_ Numeric vector of response variables
//' @param tr_ Numeric vector. The original random allocation (0s and 1s)
//' @param new_tr_mat A matrix. Each column is a new random treatment allocation with 1s (treatment group) and 0s (control group)
//' @param xb A numeric vector of fitted linear predictors
//' @param invS A matrix. If using the weighted statistic then it should be the inverse covariance matrix of the observations
//' @param family A \link{stats}[family] object
//' @param family2 A string naming the link function
//' @param nsteps Integer specifying the number of steps of the search procedure
//' @param weight Logical indicating whether to use the weighted (TRUE) or unweighted (FALSE) test statistic
//' @param alpha The function generates (1-alpha)*100% confidence intervals. Default is 0.05.
//' @param verbose Logical indicating whether to provide detailed output.
//' @return The estimated confidence interval bound
// [[Rcpp::export]]
Rcpp::List confint_search(arma::vec start,
                      arma::vec b,
                      int n,
                      arma::uword nmodel,
                      const Rcpp::List &Xnull_,
                      const Rcpp::List &y,
                      const Rcpp::NumericVector &tr_,
                      const arma::mat &new_tr_mat,
                      const Rcpp::List &invS,
                      Rcpp::List family,
                      Rcpp::List family2,
                      const arma::mat &Z,
                      const Rcpp::String &type,
                      int nsteps = 1000,
                      bool weight = true,
                      double alpha = 0.05,
                      bool verbose = true) {

  // arma::mat Xull = Rcpp::as<arma::mat>(Xnull_);
  // arma::vec y = Rcpp::as<arma::vec>(y_);
  arma::vec tr = as<arma::vec>(tr_);
  tr.replace(0,-1);
  // arma::vec xb(n,arma::fill::zeros);
  Rcpp::NumericVector weights_(n);
  weights_.fill(1);

  arma::vec bound = start;
  arma::mat boundvals(nsteps,nmodel);
  arma::vec qstat(nmodel,arma::fill::zeros);
  arma::vec qtest(nmodel);
  arma::vec step(nmodel);



  for (arma::uword i = 1; i <= nsteps; ++i) {
    boundvals.row(i-1) = bound.t();
    arma::vec new_tr = new_tr_mat.col(i-1);
    new_tr.replace(0, -1);
    for(arma::uword j = 0; j < nmodel; j++){
      arma::vec tr2 = as<arma::vec>(tr_);
      tr2.replace(1,bound(j));
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
      // const Map<VectorXd>  coefs(as<Map<VectorXd> >(result["coefficients"]));
      // Rcpp::NumericVector coefvec(wrap(coefs));
      // if(j==0)Rcpp::Rcout << "\n Xs: " << Xs;
      // if(j==0)Rcpp::Rcout << "\n cef: " << coefvec;
      // if(j==0)Rcpp::Rcout << "\n x: " << xbvec;
      arma::vec ypred = Rcpp::as<arma::vec>(linkinv(xb));
      arma::vec resids = as<arma::vec>(y[j]) - ypred;
      qstat(j) = qscore_impl(resids, tr, as<arma::vec>(xbvec), as<arma::mat>(invS[j]), as<std::string>(family2[j]), Z, weight);
      qtest(j) = qscore_impl(resids, new_tr, as<arma::vec>(xbvec), as<arma::mat>(invS[j]), as<std::string>(family2[j]), Z, weight);
    }
    // Rcpp::Rcout << "\n qstat: " << qstat.t();
    // Rcpp::Rcout << "\n qtest: " << qtest.t();
    //Rcpp::Rcout << "\nDone fits and stats" << qstat.t() << "and " << qtest.t();
    arma::uvec pos_t = arma::sort_index(qstat);

    arma::uvec rjct(nmodel);
    double k;
    double k_tmp;
    arma::vec Jval(nmodel);
    arma::vec step(nmodel);

    if(type=="rw"){
      k_tmp = gaussian_cdf(1-alpha);
      k = 2*(sqrt(M_2PI))/(k_tmp*exp((-k_tmp*k_tmp)/2));
      Jval.fill(1);
      bool reject_rest = false;
      for(arma::uword j = 0; j < nmodel; j++){
        if(!reject_rest){
          rjct(pos_t(nmodel - 1 - j)) = qstat(pos_t(nmodel-1-j)) >= arma::max(qtest(pos_t.subvec(0,(nmodel - 1 - j))));
        } else {
          rjct(pos_t(nmodel - 1 - j)) = false;
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
        Jval.fill(nmodel);
        for(arma::uword j = 0; j < nmodel; j++){
          rjct(j) = qstat(j) > qtest(j);
          step(j) = k *(b(j) - bound(j));
        }
      }
      if(type == "h"|| type=="hr"){
        bool reject_rest = false;
        for(arma::uword j = 0; j < nmodel; j++){
          if(!reject_rest){
            rjct(pos_t(nmodel - 1 - j)) = qstat(pos_t(nmodel-1-j)) >= qtest(pos_t(nmodel - 1 - j));
          } else {
            rjct(pos_t(nmodel - 1 - j)) = false;
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
        Jval.fill(1);
        for(arma::uword j = 0; j < nmodel; j++){
          rjct(j) = qstat(j) > qtest(j);
          step(j) = k *(b(j) - bound(j));
        }
      }
    }

    for(arma::uword j = 0; j < nmodel; j++){
      if(rjct(j)){
        bound(j) += step(j)*(alpha/Jval(j))/i;
      } else {
        bound(j) -= step(j)*(1-alpha/Jval(j))/i;
      }
    }

    if(verbose){
      Rcpp::Rcout << "\rStep = " << i << " bound: " << bound.t() << std::endl;
    }
  }
  return List::create(_["bound"] = bound,
                      _["values"] = boundvals);
}
