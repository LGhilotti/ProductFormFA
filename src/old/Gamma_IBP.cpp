#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;


//' Random generation from IBP with Gamma(a,b) mixture
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param a Gamma hyperparameter (shape)
//' @param b Gamma hyperparameter (rate)
//' @param n dimension of the sample to simulate
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer
//' $counts contains the counts for the observed features
//' 
//' @export
// [[Rcpp::export]]
List rGamma_IBP(double alpha, double theta, double a, double b, int n){
  
  // vector containing the number of new dishes for each customer
  std::vector<int> vec_n_new_dishes;
  vec_n_new_dishes.reserve(n);
  
  // n_dishes: it stores the current total number of selected dishes
  int n_dishes = R::rnbinom(a, b/(b+1) );
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
  double l_g_theta_1 = lgamma(theta+1);
  double l_g_alpha_theta = lgamma(alpha+theta);
  double par_0 = exp(l_g_theta_1 - l_g_alpha_theta);
  
  double l_g_theta_i = l_g_theta_1; //here i=1 ideally
  double l_g_theta_alpha_i_m1 = l_g_alpha_theta; //here i=1 ideally
  double gamma_a_t_i_m1 = par_0*exp(l_g_theta_alpha_i_m1 - l_g_theta_i); //just for convenience
  double gamma_a_t_i = gamma_a_t_i_m1; //here i=1 ideally
  
  // for each customer from the second to the n-th
  for (int i=2; i<n+1; i++){
    l_g_theta_i += log(theta+i-1);
    l_g_theta_alpha_i_m1 += log(theta+alpha+i-2);
    
    gamma_a_t_i += par_0*exp(l_g_theta_alpha_i_m1 - l_g_theta_i);
    
    // n_new_dishes: it stores the number of new dishes selected by the i-th customer
    int n_new_dishes = R::rnbinom(a+n_dishes, (gamma_a_t_i_m1 +b)/(gamma_a_t_i+b));
    gamma_a_t_i_m1 = gamma_a_t_i;
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

//////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
std::vector<double> stable_sum_M_all_gamma_IBP(double alpha,double theta,int m, int n){
  
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
// [[Rcpp::export]]
std::vector<double> p_kmn_all_gamma_IBP(double alpha,double theta,int m, int n,double b){
  
  std::vector<double> pbar;
  pbar.reserve(m);
  
  // useful quantities repeatedly used
  double l_g_theta_1 = lgamma(theta+1);
  double l_g_theta_alpha = lgamma(alpha+theta);
  double par_0 = exp(l_g_theta_1 - l_g_theta_alpha);
  
  // To track the updates in the parameters
  double l_g_theta_i = l_g_theta_1; //here i=1 ideally
  double l_g_theta_alpha_i_m1 = l_g_theta_alpha; //here i=1 ideally
  double sum_n = 0;
  for (int i=1; i<n+1; i++){
    sum_n += exp(l_g_theta_alpha_i_m1 - l_g_theta_i);
    l_g_theta_alpha_i_m1 += log(theta+alpha+i-1);
    l_g_theta_i += log(theta+i);
  }
  
  // parameter gamma_{alpha,theta}^(n) 
  double gamma_a_t_n = par_0*sum_n;
  
  double l_g_theta_n_i = l_g_theta_i - log(theta+n); //here i=0 ideally
  double l_g_theta_alpha_n_i_m1 = l_g_theta_alpha_i_m1 - log(theta+alpha+n-1); //here i=0 ideally
  double r_m = 0;
  
  for (int i=1; i<m+1; i++){
    l_g_theta_n_i += log(theta+n+i-1);
    l_g_theta_alpha_n_i_m1 += log(theta+alpha+n+i-2);
    r_m += exp(l_g_theta_alpha_n_i_m1 - l_g_theta_n_i);
    
    pbar.push_back((gamma_a_t_n + b)/(gamma_a_t_n+b+par_0*r_m));
    
  }
  
  return pbar;
}

////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double p_kmn_gamma_IBP(double alpha,double theta,int m, int n,double b){
  
  // useful quantities repeatedly used
  double l_g_theta_1 = lgamma(theta+1);
  double l_g_theta_alpha = lgamma(alpha+theta);
  double par_0 = exp(l_g_theta_1 - l_g_theta_alpha);
  
  // To track the updates in the parameters
  double l_g_theta_i = l_g_theta_1;
  double l_g_theta_alpha_i_m1 = l_g_theta_alpha;
  double sum_n = 0;
  for (int i=1; i<n+1; i++){
    sum_n += exp(l_g_theta_alpha_i_m1 - l_g_theta_i);
    l_g_theta_alpha_i_m1 += log(theta+alpha+i-1);
    l_g_theta_i += log(theta+i);
  }
  
  // parameter gamma_{alpha,theta}^(n) 
  double gamma_a_t_n = par_0*sum_n;
  
  double l_g_theta_n_i = l_g_theta_i;
  double l_g_theta_alpha_n_i_m1 = l_g_theta_alpha_i_m1; 
  double sum_m = 0;
  
  for (int i=1; i<m+1; i++){
    l_g_theta_n_i += log(theta+n+i-1);
    l_g_theta_alpha_n_i_m1 += log(theta+alpha+n+i-2);
    sum_m += exp(l_g_theta_alpha_n_i_m1 - l_g_theta_n_i);
  }
  
  return (gamma_a_t_n + b)/(gamma_a_t_n+b+par_0*sum_m);
}

/////////////////////////////////////////////////////////////////////


//' Negative Log-EFPF for IBP with Gamma mixture with reparametrization
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//'
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta+alpha in product-form feature allocation,
//' pars[2] =  value of a - Gamma hyperparameter,
//' pars[3] =  value of b - Gamma hyperparameter
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_gamma_IBP_rep(int n, std::vector<int> counts,
                                 std::vector<double> pars){
  
  double alpha = pars[0];
  double s = pars[1]; 
  double a = pars[2];
  double b = pars[3];
  
  int k = counts.size();
  
  double l_g_s = lgamma(s);
  double l_g_s_malpha_1 = lgamma(s-alpha+1);
  double par_0 = exp(l_g_s_malpha_1 - l_g_s);
  
  // To track the updates in the parameters
  double l_g_s_malpha_i = l_g_s_malpha_1; //here i=1 ideally
  double l_g_s_i_m1 = l_g_s; //here i=1 ideally
  double sum_n = 0;
  for (int i=1; i<n+1; i++){
    sum_n += exp(l_g_s_i_m1 - l_g_s_malpha_i);
    l_g_s_i_m1 += log(s+i-1);
    l_g_s_malpha_i += log(s-alpha+i);
  }
  
  double gamma_a_s_n = par_0*sum_n;
  
  double log_EFPF = a*log(b) - lgamma(k+1) + lgamma(a+k) -lgamma(a) + 
    k*(l_g_s_malpha_1 - lgamma(s-alpha+n) - log(gamma_a_s_n + b) - 
    lgamma(1-alpha) - l_g_s) - a*log(gamma_a_s_n + b);
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
  }
  
  return - log_EFPF;
  
}



