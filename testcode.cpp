#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]

int some_function(int x) {
  Rcpp::NumericVector rep_vec1 = Rcpp::rep(1.0,x);
  
  // print('the size is ');
  int num = rep_vec1.size();
  
  return num;
    
}

//= Rcpp::rep(1.0,11);

