#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;

//' Buffet procedure for BB with Negative-Binomial mixture from beginning
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param n dimension of the sample to simulate
//' @param nstar Negative-Binomial hyperparameter (number of successes)
//' @param p Negative-Binomial hyperparameter (success probability)
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer,
//' $counts contains the counts for the observed features
//' 
//' @export
// [[Rcpp::export]]
List buffet_negbin_BB(double alpha,double theta,int n,int nstar, double p) {
  
  // vector containing the number of new dishes for each customer
  std::vector<int> vec_n_new_dishes;
  vec_n_new_dishes.reserve(n);
  
  // n_dishes: it stores the current total number of selected dishes
  int n_dishes = R::rnbinom(nstar, 1+alpha*(1-p)/(p*theta -alpha*(1-p)) );
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
    
    double par_0 = exp(l_g_alpha_theta + l_g_theta_i);
    double par_1 = exp(l_g_theta_alpha_i_m1 + l_g_theta);
    double pbar = 1 + alpha*(1-p)*par_1/(par_0 - (1-p)*(theta+alpha+i-1)*par_1);

    // n_new_dishes: it stores the number of new dishes selected by the i-th customer
    int n_new_dishes = R::rnbinom(nstar+n_dishes, pbar);
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

//////////////////////////////////////////////////////////////////////

//' Buffet procedure for BB with Negative-Binomial mixture given initial sample
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param m dimension of the new sample to be observed
//' @param n dimension of the already observed sample
//' @param counts vector of counts for the already observed features
//' @param nstar Negative-Binomial hyperparameter (number of successes)
//' @param p Negative-Binomial hyperparameter (success probability)
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer
//' $counts contains the counts for the observed features
//' 
//' @export
// [[Rcpp::export]]
List buffet_negbin_BB_initial_sample(double alpha,double theta,int m, int n, std::vector<int> counts, int nstar, double p) {
  
  // vector containing the number of new dishes for each customer
  std::vector<int> vec_n_new_dishes;
  vec_n_new_dishes.reserve(m);
  
  // n_dishes: it stores the current total number of selected dishes so far
  int n_dishes = counts.size();
  
  // build List (Rcpp) with each element a new customer, associated to a vector of dishes
  List dishes(m);
  
  // useful quantities repeatedly used 
  double l_g_theta = lgamma(theta);
  double l_g_alpha_theta = lgamma(alpha+theta);
  double l_g_theta_n_i = lgamma(theta+n); //here i=0 ideally
  double l_g_theta_alpha_n_i_m1 = lgamma(theta+alpha+n-1); //here i=0 ideally
  
  // for each new customer from the n+1-th to the m-th
  for (int i=1; i<m+1; i++){
    l_g_theta_n_i += log(theta+n+i-1);
    l_g_theta_alpha_n_i_m1 += log(theta+alpha+n+i-2);
    
    double par_0 = exp(l_g_alpha_theta + l_g_theta_n_i);
    double par_1 = exp(l_g_theta_alpha_n_i_m1 + l_g_theta);
    double pbar = 1 + alpha*(1-p)*par_1/(par_0 - (1-p)*(theta+alpha+n+i-1)*par_1);
    
    // n_new_dishes: it stores the number of new dishes selected by the i-th customer
    int n_new_dishes = R::rnbinom(nstar+n_dishes, pbar);
    vec_n_new_dishes.push_back(n_new_dishes);
    
    // prob_old_i: it stores the probs that the old dishes are served to the i-th customer
    std::vector<double> prob_old_i(n_dishes);
    std::transform(counts.cbegin(), counts.cend(), prob_old_i.begin(), [i,n,alpha,theta](int c){return (c - alpha)/(theta+n+i-1);});
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
                      Named("counts") = counts) ;
}

//////////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////


//' Negative Log-EFPF for BB with Negative-Binomial mixture with reparametrization
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//'
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta+alpha in product-form feature allocation,
//' pars[2] =  value of nstar - NegBin hyperparameter,
//' pars[3] = value of p - NegBin hyperparameter
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_BB_rep(int n, std::vector<int> counts,
                                 std::vector<double> pars){
  
  double alpha = pars[0];
  double s = pars[1]; 
  double nstar = pars[2];
  double p = pars[3];
  
  int k = counts.size();
  
  double l_g_s = lgamma(s);
  double l_g_s_malpha = lgamma(s-alpha);
  double l_g_s_malpha_n = lgamma(s-alpha +n);
  
  double par_0 = exp(lgamma(s+n) - l_g_s + l_g_s_malpha - l_g_s_malpha_n);
  
  double log_EFPF = lgamma(k+nstar) - lgamma(k+1) - lgamma(nstar) + nstar*log(p) + 
    k*(log(1-p) + log(-alpha) + l_g_s_malpha - l_g_s_malpha_n - lgamma(1-alpha) - 
    l_g_s - log(1-(1-p)*par_0)) - nstar*log(1-(1-p)*par_0) ;
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
  }
  
  return - log_EFPF;
  
}

