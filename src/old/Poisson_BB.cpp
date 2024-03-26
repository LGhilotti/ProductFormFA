#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;

//' Random generation from BB with Poisson(lambda) mixture
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param lambda Poisson hyperparameter 
//' @param n dimension of the sample to generate
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer
//' $counts contains the counts for the observed features
//' 
//' @export
// [[Rcpp::export]]
List rPoisson_BB(double alpha,double theta,double lambda, int n){
  
  // vector containing the number of new dishes for each customer
  std::vector<int> vec_n_new_dishes;
  vec_n_new_dishes.reserve(n);
  
  // n_dishes: it stores the current total number of selected dishes
  int n_dishes = R::rpois( -lambda*alpha/theta );
  vec_n_new_dishes.push_back(n_dishes);
  
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
    double l_par = l_g_theta + l_g_theta_alpha_i_m1 - l_g_alpha_theta - l_g_theta_i;
    
    // n_new_dishes: it stores the number of new dishes selected by the i-th customer
    int n_new_dishes = R::rpois(-lambda*alpha*exp(l_par));
    vec_n_new_dishes.push_back(n_new_dishes);
    
    // prob_old_i: it stores the probs that the old dishes are served to the i-th customer
    std::vector<double> prob_old_i(n_dishes);
    std::transform(counts.cbegin(), counts.cend(), prob_old_i.begin(), [i,alpha,theta](int c){return (c - alpha)/(theta+i-1);});
    // old_observed_i: it stores 0/1 for the old dishes, indicating if i-th tries the old dishes
    std::vector<int> old_observed_i = cpp_rbern(n_dishes, prob_old_i);
    
    // n_old_observed_i: it stores the number of old dishes also tried by i-th customer
    int n_old_observed_i = std::accumulate(old_observed_i.begin(), old_observed_i.end(), 0);
    
    if (n_old_observed_i + n_new_dishes > 0){
      // dishes_i: it stores the dishes served to the i-th customer #(n_old_observed_i + n_new_dishes)
      std::vector<int> dishes_i(n_old_observed_i + n_new_dishes);
  
      if (n_old_observed_i >0){
        std::vector<int>::iterator it = old_observed_i.begin();
        int j=0;
        while ((it = std::find_if(it, old_observed_i.end(), [](int x){return x == 1; })) != old_observed_i.end())
        {
          dishes_i[j] = std::distance(old_observed_i.begin(), it)+1;
          it++;
          j++;
        }
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

  }

  return List::create(Named("features") = dishes, Named("num_new") = vec_n_new_dishes,
                      Named("counts") = counts);
}

// SISTEMARE

///////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
std::vector<double> mean_kmn_all_poiss_BB(double alpha,double theta,int m, int n,double lambda){
  
  std::vector<double> means;
  means.reserve(m);
  
  // useful quantities repeatedly used
  double l_g_theta = lgamma(theta);
  double l_g_theta_alpha = lgamma(alpha+theta);
  double l_g_theta_n = lgamma(theta+n);
  double l_g_theta_alpha_n = lgamma(alpha+theta+n);
  
  double par_0 = lambda*exp(l_g_theta - l_g_theta_alpha);
  double par_1 = exp(l_g_theta_alpha_n - l_g_theta_n);
  
  // To track the updates in the parameters
  double l_g_theta_n_i = l_g_theta_n;
  double l_g_theta_alpha_n_i = l_g_theta_alpha_n;
  
  for (int i=1; i<m+1; i++){
    l_g_theta_alpha_n_i += log(theta+alpha+n+i-1);
    l_g_theta_n_i += log(theta+n+i-1);
    
    means.push_back(par_0 *(par_1 - exp(l_g_theta_alpha_n_i - l_g_theta_n_i)));
  }
  
  return means;
}

////////////////////////////////////////////////////////////////////

//' @export
// [[Rcpp::export]]
double mean_kmn_poiss_BB(double alpha,double theta,int m, int n,double lambda){
  
  // useful quantities repeatedly used
  double l_g_theta = lgamma(theta);
  double l_g_theta_alpha = lgamma(alpha+theta);
  double l_g_theta_n = lgamma(theta+n);
  double l_g_theta_alpha_n = lgamma(alpha+theta+n);
  double l_g_theta_n_m = lgamma(theta+n+m);
  double l_g_theta_alpha_n_m = lgamma(alpha+theta+n+m); 
  
  double par_0 = lambda*exp(l_g_theta - l_g_theta_alpha);
  double par_1 = exp(l_g_theta_alpha_n - l_g_theta_n);
  double par_2 = exp(l_g_theta_alpha_n_m - l_g_theta_n_m);
  
  return par_0*(par_1 - par_2);
}

/////////////////////////////////////////////////////////////////////


//' Negative Log-EFPF for BB with Poisson mixture with reparametrization
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//'
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta+alpha in product-form feature allocation,
//' pars[2] =  value of lambda - Poisson hyperparameter
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_poiss_BB_rep(int n, std::vector<int> counts,
                                 std::vector<double> pars){
  
  double alpha = pars[0];
  double s = pars[1]; 
  double lambda = pars[2];
  
  int k = counts.size();
  
  double l_g_s = lgamma(s);
  double l_g_s_malpha = lgamma(s-alpha);
  double l_g_1_malpha = lgamma(1-alpha);
  double l_g_s_malpha_n = lgamma(s-alpha +n);
  
  double par_0 = exp(lgamma(s+n) - l_g_s + l_g_s_malpha - l_g_s_malpha_n);
  
  double log_EFPF = - lgamma(k+1) - lambda*(1- par_0 ) + 
    k*(log(-lambda*alpha) - l_g_s_malpha_n + l_g_s_malpha - l_g_1_malpha - l_g_s);
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
  }
  
  return - log_EFPF;
  
}


