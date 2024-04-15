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
//' 
// [[Rcpp::export]]
List rPoissonBB(double alpha,double theta,double lambda, int n){
 
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


//' Random generation from BB with NB(n0,mu0) mixture
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param n0 hyperparameter "n0" of NB(n0,mu0) for N
//' @param mu0 hyperparameter "mu0" of NB(n0,mu0) for N
//' @param n dimension of the sample to simulate
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer,
//' $counts contains the counts for the observed features
//' 
//' @export
//' 
// [[Rcpp::export]]
List rNegBinBB(double alpha,double theta,int n0, double mu0, int n){
  
  int nstar = n0;
  double p = 1/(mu0/n0 + 1);

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



//' Random generation from IBP with Gamma(a,b) mixture
//' 
//' @param alpha value of alpha in product-form feature allocation
//' @param theta value of theta in product-form feature allocation
//' @param a hyperparameter "a" of Gamma(a,b) for gamma
//' @param b hyperparameter "b" of Gamma(a,b) for gamma
//' @param n dimension of the sample to simulate
//' 
//' @return list: $features contains the simulated features for each customer,
//' $num_new contains the number of new features selected for each customer
//' $counts contains the counts for the observed features
//' 
//' @export
//' 
// [[Rcpp::export]]
List rGammaIBP(double alpha, double theta, double a, double b, int n){
 
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
