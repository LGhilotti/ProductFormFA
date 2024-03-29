#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> stable_sum_M_all_gamma_IBP(double alpha,double theta,int m, bool only_last, int n){
  
  if (only_last){
    
    std::vector<double> sum_M; 
    sum_M.resize(1);
    
    double par = 0;
    for (int h =1; h<m+1; h++){
      par += exp(lgamma(alpha + theta + n + h - 1) -
        lgamma(alpha + theta) -
        lgamma(theta + n + h) +
        lgamma(theta + 1) ) ;
    }
    sum_M[0] = par;
    
    return sum_M;
    
    
  } else {
    
    std::vector<double> sum_M; 
    sum_M.resize(m);
    
    sum_M[0] = exp(lgamma(alpha + theta + n) -
      lgamma(alpha + theta) -
      lgamma(theta + n + 1) +
      lgamma(theta + 1) )  ;
    
    for (int j=2; j < m+1; j++){
      double par = 0;
      for (int h =1; h<j+1; h++){
        par += exp(lgamma(alpha + theta + n + h - 1) -
          lgamma(alpha + theta) -
          lgamma(theta + n + h) +
          lgamma(theta + 1) ) ;
      }
      sum_M[j-1] = par;
    }
    
    return  sum_M;
    
  }
  
}