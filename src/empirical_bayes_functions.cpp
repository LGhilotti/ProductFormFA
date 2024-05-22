#include<numeric>
#include <Rcpp.h>
using namespace Rcpp;

//' Negative Log-EFPF for the BB model  
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta + alpha in product-form feature allocation,
//' pars[2] =  value of Nhat_prime = Nhat - k 
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_BB(int n, std::vector<int> counts,
                       std::vector<double> pars ){
 
 double alpha = pars[0];
 double s = pars[1]; 
 double Nhat_prime = pars[2];
 
 int k = counts.size();
 double Nhat = Nhat_prime + k;
 
 
 double l_g_s = lgamma(s);
 double l_g_s_n = lgamma(s + n);
 double l_g_s_malpha = lgamma(s-alpha);
 double l_g_1_malpha = lgamma(1-alpha);
 double l_g_s_malpha_n = lgamma(s-alpha +n);
 
 double log_EFPF = lgamma(Nhat + 1) - lgamma(k + 1) - lgamma(Nhat - k + 1) +
   k*log(-alpha) + (Nhat - k)*l_g_s_n - 
   Nhat*(l_g_s + l_g_s_malpha_n - l_g_s_malpha) - k*l_g_1_malpha;
 
 for (int l=0; l<k;l++){
   log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
 }
 
 return - log_EFPF;
 
}
 
//' Negative Log-EFPF for the IBP model  
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta + alpha in product-form feature allocation,
//' pars[2] =  value of Gamma 
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_IBP(int n, std::vector<int> counts,
                      std::vector<double> pars ){
 
 double alpha = pars[0];
 double s = pars[1]; 
 double Gam = pars[2];
 
 int k = counts.size();

 double l_g_s = lgamma(s);
 double l_g_1_malpha = lgamma(1-alpha);
 double l_g_s_malpha_n = lgamma(s-alpha +n);
 double l_g_s_malpha_1 = lgamma(s-alpha + 1);
 double gn = 0;
 for (int i=1; i<n + 1; i++){
   gn += exp(lgamma(s + i - 1) - lgamma(s - alpha + i));
 }
 gn = gn* exp(l_g_s_malpha_1 - l_g_s);
 
 double log_EFPF = - lgamma(k+ 1) + 
   k*(log(Gam) + l_g_s_malpha_1 - l_g_s_malpha_n - l_g_1_malpha - l_g_s) -
   Gam*gn;
 
 for (int l=0; l<k;l++){
   log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
 }
 
 return - log_EFPF;
 
}


//' Negative Log-EFPF for the BB model with Poisson(lambda) mixture 
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta + alpha in product-form feature allocation,
//' pars[2] =  value of lambda - Poisson hyperparameter
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_PoissonBB(int n, std::vector<int> counts,
                              std::vector<double> pars ){
 
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


//' Negative Log-EFPF for BB with NB(n0,p) mixture 
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta + alpha in product-form feature allocation,
//' pars[2] = value of n0 - NB hyperparameter, first
//' pars[3] = value of p - NB hyperparameter, second
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_NegBinBB(int n, std::vector<int> counts,
                             std::vector<double> pars){
 
 double alpha = pars[0];
 double s = pars[1];
 double n0 = pars[2];
 double p = pars[3];
 
 int k = counts.size();
 
 double l_g_s = lgamma(s);
 double l_g_s_malpha = lgamma(s-alpha);
 double l_g_s_malpha_n = lgamma(s-alpha +n);
 
 double par_0 = exp(lgamma(s+n) - l_g_s + l_g_s_malpha - l_g_s_malpha_n);
 
 double log_EFPF = lgamma(k+n0) - lgamma(k+1) - lgamma(n0) + n0*log(p) + 
   k*(log(1-p) + log(-alpha) + l_g_s_malpha - l_g_s_malpha_n - lgamma(1-alpha) - 
   l_g_s - log(1-(1-p)*par_0)) - n0*log(1-(1-p)*par_0) ;
 
 for (int l=0; l<k;l++){
   log_EFPF += lgamma(counts[l]-alpha) + lgamma(s+n-counts[l]);
 }
 
 return - log_EFPF;
 
}


//' Negative Log-EFPF for IBP with Gamma(a,b) mixture
//' 
//' @param n dimension of the observed sample
//' @param counts vector of cardinalities for the observed features
//' @param pars pars[0] = value of alpha in product-form feature allocation,
//' pars[1] = value of s = theta + alpha in product-form feature allocation,
//' pars[2] =  value of a - Gamma hyperparameter,
//' pars[3] =  value of b - Gamma hyperparameter
//' 
//' @return value of the negative logarithm of the EFPF for the sample of 
//' dimensionality n described by counts
//' 
// [[Rcpp::export]]
double neg_log_EFPF_GammaIBP(int n, std::vector<int> counts,
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
