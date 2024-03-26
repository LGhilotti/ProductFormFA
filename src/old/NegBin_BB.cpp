#include<numeric>
#include <Rcpp.h>
#include "cpp_probability_generators.hpp"
using namespace Rcpp;

//' Random generation from BB with NB(nstar,p) mixture
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param nstar Negative-Binomial hyperparameter (number of successes)
//' @param p Negative-Binomial hyperparameter (success probability)
//' @param n dimension of the sample to simulate
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer,
//' $counts contains the counts for the observed features
//' 
//' @export
// [[Rcpp::export]]
List rNegBin_BB(double alpha,double theta,int nstar, double p, int n){
   
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
   
   // for each customer from the second to the n-th
   for (int i=2; i<n+1; i++){
     double l_g_theta_i = lgamma(theta+i);
     double l_g_theta_alpha_i_m1 = lgamma(theta+alpha+i-1);
     double l_g_theta_alpha_i = lgamma(theta+alpha+i);
     
     double inv_l_1 = l_g_alpha_theta - log(-alpha*(1-p)) - l_g_theta + l_g_theta_i - l_g_theta_alpha_i_m1 ;
     double inv_l_2 = l_g_theta_alpha_i - l_g_theta_alpha_i_m1 - log(-alpha);
     
     double pbar = 1 - 1/(exp(inv_l_1) - exp(inv_l_2));
     
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


// SISTEMARE

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
  
  
  //' Negative Log-EFPF for BB with Negative-Binomial mixture NO reparametrization
//' 
//' @param pars vector of parameters to optimize
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param nopt_par vector of parameters not to optimize
//' @param opt boolean vector indicating if the variable is optimized
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
//' @export 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_BB(std::vector<double> pars,
                              int n, std::vector<int> counts, 
                              std::vector<double> nopt_par,
                              std::vector<bool> opt){
  
  double alpha; double theta; double nstar; double p;
  
  int i_T=0; int i_F=0;
  // Assign alpha
  if (opt[0]==TRUE){
    alpha = pars[0];
    i_T ++;
  }
  else {
    alpha = nopt_par[0]; 
    i_F ++;
  }
  // Assign theta
  if (opt[1]==TRUE){
    theta = pars[i_T];
    i_T ++;
  }
  else{
    theta = nopt_par[i_F];
    i_F ++;
  }
  // Assign nstar
  if (opt[2]==TRUE){
    nstar = pars[i_T];
    i_T ++;
  }
  else {
    nstar = nopt_par[i_F]; 
    i_F ++;
  }
  // Assign p
  if (opt[3]==TRUE){
    p = pars[i_T];
  }
  else{
    p = nopt_par[i_F];
  }
  
  
  int k = counts.size();
  
  double l_g_alpha_theta = lgamma(alpha+theta);
  double l_g_theta = lgamma(theta);
  double l_g_theta_n = lgamma(theta +n);
  
  double par_0 = exp(lgamma(alpha+theta+n) - l_g_alpha_theta + l_g_theta - l_g_theta_n);
  
  double log_EFPF = lgamma(k+nstar) - lgamma(k+1) - lgamma(nstar) + nstar*log(p) + 
    k*(log(1-p) + log(-alpha) + l_g_theta - l_g_theta_n - lgamma(1-alpha) - 
         l_g_alpha_theta - log(1-(1-p)*par_0)) - nstar*log(1-(1-p)*par_0) ;
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(alpha+theta+n-counts[l]);
  }
  
  return - log_EFPF;
  
}


/////////////////////////////////////////////////////////////////////
  
  
  //' Negative Log-EFPF for BB with Negative-Binomial mixture with reparametrization
//' 
//' @param pars vector of parameters to optimize
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param nopt_par vector of parameters not to optimize
//' @param opt boolean vector indicating if the variable is optimized
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
//' @export 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_BB_rep(std::vector<double> pars,
                                  int n, std::vector<int> counts, 
                                  std::vector<double> nopt_par,
                                  std::vector<bool> opt){
  
  double alpha = pars[0];
  double s = pars[1];
  double nstar; double p;
  
  int i_T=0; int i_F=0;
  // Assign nstar
  if (opt[2]==TRUE){
    nstar = pars[2];
    i_T ++;
  }
  else {
    nstar = nopt_par[0]; 
    i_F ++;
  }
  // Assign p
  if (opt[3]==TRUE){
    p = pars[i_T+2];
  }
  else{
    p = nopt_par[i_F];
  }
  
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


/////////////////////////////////////////////////////////////////////
  
  
  //' Negative Log-EFPF for BB with Negative-Binomial mixture with reparametrization
//' (when all the parameters are optimized)
//' 
//' @param pars vector of parameters to optimize
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
//' @export 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_BB_rep_all(std::vector<double> pars,
                                  int n, std::vector<int> counts){
  
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

///////////////////////////////////////////////////////////////////////////
///////// functions with NegBin parametrized via mean and variance ////////
///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////


//' Negative Log-EFPF for BB with Negative-Binomial mixture NO reparametrization -
//' - (mean-variance parametrization)
//' 
//' @param pars vector of parameters to optimize
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param nopt_par vector of parameters not to optimize
//' @param opt boolean vector indicating if the variable is optimized
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
//' @export 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_mv_BB(std::vector<double> pars,
                              int n, std::vector<int> counts, 
                              std::vector<double> nopt_par,
                              std::vector<bool> opt){
  
  double alpha; double theta; double mu; double sigma2;
  
  int i_T=0; int i_F=0;
  // Assign alpha
  if (opt[0]==TRUE){
    alpha = pars[0];
    i_T ++;
  }
  else {
    alpha = nopt_par[0]; 
    i_F ++;
  }
  // Assign theta
  if (opt[1]==TRUE){
    theta = pars[i_T];
    i_T ++;
  }
  else{
    theta = nopt_par[i_F];
    i_F ++;
  }
  // Assign the mu 
  if (opt[2]==TRUE){
    mu = pars[i_T];
    i_T ++;
  }
  else {
    mu = nopt_par[i_F]; 
    i_F ++;
  }
  // Assign sigma2
  if (opt[3]==TRUE){
    sigma2 = pars[i_T];
  }
  else{
    sigma2 = nopt_par[i_F];
  }
  
  
  int k = counts.size();
  
  double l_g_alpha_theta = lgamma(alpha+theta);
  double l_g_theta = lgamma(theta);
  double l_g_theta_n = lgamma(theta +n);
  
  double par_0 = exp(lgamma(alpha+theta+n) - l_g_alpha_theta + l_g_theta - l_g_theta_n);
  
  double log_EFPF = lgamma(k+(pow(mu,2))/(sigma2-mu)) - lgamma(k+1) - 
    lgamma((pow(mu,2))/(sigma2-mu)) + (pow(mu,2))/(sigma2-mu)*(log(mu) - log(sigma2)) + 
    k*(log(sigma2-mu) - log(sigma2) + log(-alpha) + l_g_theta - 
    l_g_theta_n - lgamma(1-alpha) - l_g_alpha_theta - 
    log(1-(sigma2 -mu)/sigma2*par_0)) - 
    (pow(mu,2))/(sigma2-mu)*log(1-(sigma2 -mu)/sigma2*par_0) ;
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(alpha+theta+n-counts[l]);
  }
  
  return - log_EFPF;
  
}


//' Negative Log-EFPF for BB with Negative-Binomial mixture with reparametrization
//' (and mean-variance parametrization)
//' 
//' @param pars vector of parameters to optimize
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param nopt_par vector of parameters not to optimize
//' @param opt boolean vector indicating if the variable is optimized
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
//' @export 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_mv_BB_rep(std::vector<double> pars,
                                     int n, std::vector<int> counts, 
                                     std::vector<double> nopt_par,
                                     std::vector<bool> opt){
  
  double alpha; double s; double mu; double t;
  
  int i_T=0; int i_F=0;
  // Assign alpha
  if (opt[0]==TRUE){
    alpha = pars[0];
    i_T ++;
  }
  else {
    alpha = nopt_par[0]; 
    i_F ++;
  }
  // Assign s
  if (opt[1]==TRUE){
    s = pars[i_T];
    i_T ++;
  }
  else{
    s = nopt_par[i_F];
    i_F ++;
  }
  // Assign mu
  if (opt[2]==TRUE){
    mu = pars[2];
    i_T ++;
  }
  else {
    mu = nopt_par[0]; 
    i_F ++;
  }
  // Assign t
  if (opt[3]==TRUE){
    t = pars[i_T+2];
  }
  else{
    t = nopt_par[i_F];
  }
  
  int k = counts.size();
  
  double l_g_s = lgamma(s);
  double l_g_s_malpha = lgamma(s-alpha);
  double l_g_s_malpha_n = lgamma(s-alpha +n);
  
  double par_0 = exp(lgamma(s+n) - l_g_s + l_g_s_malpha - l_g_s_malpha_n);
  
  double log_EFPF = lgamma(k+ pow(mu,2)/t) - lgamma(k+1) - 
    lgamma( pow(mu,2)/ t) + pow(mu,2)/ t *(log(mu) - log(mu + t)) + 
    k*(log(t) - log(t+mu) + log(-alpha) + l_g_s_malpha - 
    l_g_s_malpha_n - lgamma(1-alpha) - l_g_s - 
    log(1- t/(mu+t)*par_0)) - 
    pow(mu,2)/t*log(1-t/(mu+t)*par_0) ;
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
  }
  
  return - log_EFPF;
  
}


/////////////////////////////////////////////////////////////////////
  
  
  //' Negative Log-EFPF for BB with Negative-Binomial mixture with reparametrization
//' (when all the parameters are optimized -  mean-variance parametrization)
//' 
//' @param pars vector of parameters to optimize
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
//' @export 
// [[Rcpp::export]]
double neg_log_EFPF_negbin_mv_BB_rep_all(std::vector<double> pars,
                                      int n, std::vector<int> counts){
  
  double alpha = pars[0];
  double s = pars[1];
  double mu = pars[2]; 
  double t = pars[3];
  
  int k = counts.size();
  
  double l_g_s = lgamma(s);
  double l_g_s_malpha = lgamma(s-alpha);
  double l_g_s_malpha_n = lgamma(s-alpha +n);
  
  double par_0 = exp(lgamma(s+n) - l_g_s + l_g_s_malpha - l_g_s_malpha_n);
  
  double log_EFPF = lgamma(k+ pow(mu,2)/t) - lgamma(k+1) - 
    lgamma( pow(mu,2)/ t) + pow(mu,2)/ t *(log(mu) - log(mu + t)) + 
    k*(log(t) - log(t+mu) + log(-alpha) + l_g_s_malpha - 
    l_g_s_malpha_n - lgamma(1-alpha) - l_g_s - 
    log(1- t/(mu+t)*par_0)) - 
    pow(mu,2)/t*log(1-t/(mu+t)*par_0) ;
  
  for (int l=0; l<k;l++){
    log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
  }
  
  return - log_EFPF;
  
}

