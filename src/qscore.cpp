#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double qscore(NumericVector &y,
              NumericVector &x,
              NumericVector &Tr,
              IntegerVector &cl,
              int ncl) {

  int n = y.size();
  NumericVector pr1(n);
  NumericVector cl_sum (ncl);
  double tot_sum = 0;
  double cl_total=0;

  pr1 = y - x;

  for(int i=0; i<n; i++){
    if(Tr[i]==1){
      tot_sum += pr1[i];
    } else {
      tot_sum += -1*pr1[i];
    }
    cl_sum[cl[i]] += pr1[i];
  }


  for(int j=0;j<ncl;j++){
    cl_total += pow(cl_sum[j],2);
  }

  return tot_sum/sqrt(cl_total);
}

// [[Rcpp::export]]
double qscorew(const arma::vec &y,
               const arma::vec &x,
               const arma::mat &Tr,
               const arma::vec &g,
               const arma::mat &sigma,
               const arma::uvec &cl,
               arma::uword ncl){

    arma::vec gvu(ncl);
    arma::vec pr1 = y - x;
    double denom = 0;
    double total = 0;
    arma::vec gv = sigma * g;

    for(arma::uword j = 0; j <ncl; ++j){
      arma::uvec ids = find(cl == j);
      gvu(j) = arma::conv_to<double>::from(gv(ids).t() * Tr(ids,ids) * pr1(ids));
      total += gvu(j);
      denom += gvu(j)*gvu(j);
    }

    //Rcout << "gvu : " << gvu << "\n";

    return total/sqrt(denom);

}
