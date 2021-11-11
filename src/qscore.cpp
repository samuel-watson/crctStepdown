#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double qscore(NumericVector y,
              NumericVector x,
              NumericVector T,
              IntegerVector cl,
              int ncl) {

  int n = y.size();
  NumericVector pr1(n);
  NumericVector cl_sum (ncl);
  double tot_sum = 0;
  double cl_total=0;

  pr1 = y - x;

  for(int i=0; i<n; i++){
    if(T[i]==1){
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



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

