#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List buffet_negbin_BB(double alpha,double theta,int n,int nstar, double p) {
  
  // n_dishes: it stores the current total number of selected dishes
  int n_dishes = R::rnbinom(nstar, 1+alpha*(1-p)/(p*theta -alpha*(1-p)) );
  
  // dishes_first: it stores the dishes selected by the first customer
  std::vector<int> dishes_first(n_dishes);
  std::iota(dishes_first.begin(), dishes_first.end(),1);
  
  // counts: vector of counts of the served dishes
  std::vector<int> counts(n_dishes,1);
  
  // build List (Rcpp) with each element a customer, associated to a vector of dishes
  List dishes(n);
  dishes[0] = dishes_first;
  
  // useful quantities repeatedly used
  double l_g_theta = lgamma(theta);
  double l_g_alpha_theta = lgamma(alpha+theta);
  double l_g_theta_i = lgamma(theta+1);
  double l_g_theta_alpha_i_m1 = l_g_alpha_theta;
  
  // for each customer from the second to the n-th
  for (int i=2; i<n+1; i++){
    l_g_theta_i += log(theta+i-1);
    l_g_theta_alpha_i_m1 += log(theta+alpha+i-2);
    
    double par_0 = exp(l_g_alpha_theta + l_g_theta_i);
    double par_1 = exp(l_g_theta_alpha_i_m1 + l_g_theta);
    double pbar = 1 + alpha*(1-p)*par_1/(par_0 - (1-p)*(theta+alpha+i-1)*par_1);

    // n_new_dishes: it stores the number of new dishes selected by the i-th customer
    int n_new_dishes = R::rnbinom(nstar+n_dishes, pbar);
    
    // prob_old_i: it stores the probs that the old dishes are served to the i-th customer
    std::vector<double> prob_old_i(n_dishes);
    std::transform(counts.cbegin(), counts.cend(), prob_old_i.begin(), [i,alpha,theta](int c){return (c - alpha)/(theta+i-1);});
    // old_observed_i: it stores 0/1 for the old dishes, indicating if i-th tries the old dishes
    std::vector<int> old_observed_i = cpp_rbern(n_dishes, prob_old_i);
    
    // n_old_observed_i: it stores the number of old dishes also tried by i-th customer
    int n_old_observed_i = std::accumulate(old_observed_i.begin(), old_observed_i.end(), 0);
    
    // dishes_i: it stores the dishes served to the i-th customer #(n_old_observed_i + n_new_dishes)
    std::vector<int> dishes_i(n_old_observed_i + n_new_dishes);
    
    std::vector<int>::iterator it = old_observed_i.begin();
    int j=0;
    while ((it = std::find_if(it, old_observed_i.end(), [](int x){return x == 1; })) != old_observed_i.end())
    {
      dishes_i[j] = std::distance(old_observed_i.begin(), it)+1;
      it++;
      j++;
    }
    std::iota(dishes_i.begin()+n_old_observed_i, dishes_i.end(),n_dishes+1);
    
    // Update total number of served dishes after i-th customer
    n_dishes += n_new_dishes;
    // Update counts of served dishes after i-th customer
    std::transform (counts.begin(), counts.end(), old_observed_i.begin(), counts.begin(), std::plus<int>());
    counts.reserve(n_dishes);
    counts.insert(counts.end(), n_new_dishes, 1);
    // Update the list of customers with i-th customer and his dishes
    dishes[i-1] = dishes_i;
    
  }
  
  return dishes;
}


// [[Rcpp::export]]
std::vector<double> p_kmn_all_negbin_BB(double alpha,double theta,int m, int n,double p){
  
  std::vector<double> pbar;
  pbar.reserve(m);
  
  // useful quantities repeatedly used
  double l_g_theta = lgamma(theta);
  double l_g_theta_alpha = lgamma(alpha+theta);
  double l_g_theta_n = lgamma(theta+n);
  double l_g_theta_alpha_n = lgamma(alpha+theta+n);
  
  double par_0 = (1-p)*exp(l_g_theta - l_g_theta_alpha);
  double par_1 = exp(l_g_theta_alpha_n - l_g_theta_n);
  
  // To track the updates in the parameters
  double l_g_theta_n_i = l_g_theta_n;
  double l_g_theta_alpha_n_i = l_g_theta_alpha_n;
  
  for (int i=1; i<m+1; i++){
    l_g_theta_alpha_n_i += log(theta+alpha+n+i-1);
    l_g_theta_n_i += log(theta+n+i-1);
    double par_2 = exp(l_g_theta_alpha_n_i - l_g_theta_n_i);
    
    pbar.push_back(1- par_0 *(par_1 - par_2)/(1-par_0*par_2));
  }
  
  return pbar;
}


// [[Rcpp::export]]
double p_kmn_negbin_BB(double alpha,double theta,int m, int n,double p){
  
  // useful quantities repeatedly used
  double l_g_theta = lgamma(theta);
  double l_g_theta_alpha = lgamma(alpha+theta);
  double l_g_theta_n = lgamma(theta+n);
  double l_g_theta_alpha_n = lgamma(alpha+theta+n);
  double l_g_theta_n_m = lgamma(theta+n+m);
  double l_g_theta_alpha_n_m = lgamma(alpha+theta+n+m); 
  
  double par_0 = (1-p)*exp(l_g_theta - l_g_theta_alpha);
  double par_1 = exp(l_g_theta_alpha_n - l_g_theta_n);
  double par_2 = exp(l_g_theta_alpha_n_m - l_g_theta_n_m);
  
  return 1- par_0 *(par_1 - par_2)/(1-par_0*par_2);
}

