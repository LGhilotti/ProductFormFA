#include<numeric>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<int> cpp_rbern(int n, std::vector<double> prob) { 
  std::vector<int> v(n);
  std::transform( prob.begin(), prob.end(), v.begin(), [=](double p){ return R::rbinom(1, p); }); 
  return(v);
}