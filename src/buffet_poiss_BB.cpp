#include<numeric>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
IntegerVector buffet_poiss_BB(double alpha,double theta,int n,double lambda) {
  
  
  
  int n_dish = R::rpois( -lambda*alpha/theta );
  //std::vector<int> features_first(n_dish);
  //std::iota(features_first.begin(), features_first.end(),1);
  IntegerVector features_first(n_dish);
  for(int j = 1; j < n_dish+1; j++) {
    // Set the jth element of sequence to j
    int k=j-1;
    features_first[k] =  j;
  }

  // build DataFrame (Rcpp) with each element a customer, associated to a vector of indexes of the features
  //std::vector<std::vector<int>> feat(n);
  //feat[0] = features_first;
  
  std::vector<int> a = {1,2};
  
  //feat.push_back(a);
  
  return features_first;
}


/*** R
alpha=-1
theta=2
n=5
lambda=10

buffet_poiss_BB(alpha,theta,n,lambda)
*/