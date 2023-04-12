#ifndef CPP_PROB_GEN
#define CPP_PROB_GEN

#include<numeric>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<int> cpp_rbern(int n, std::vector<double> prob);
#endif
