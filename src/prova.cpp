#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
std::vector<int> prova() {
  std::vector<int> counts = {1,4,2};
  counts.resize(5);
  std::iota(counts.begin()+3, counts.end() ,6);
  
  
  return counts;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
L = prova()
*/
