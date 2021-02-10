// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(StanHeaders)]]
#include <stan/math.hpp>
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

/**
This function should compute the gradient of the KPR likelihood at vector theta using AD.
**/

 // [Rcpp::export]
 VectorXd gradient_auto(VectorXd theta_) // this needs to also take arguments corresponding to X, Y, H, and Q
 {

   std::vector<stan::math::var> theta(theta_.size()); // create a new vector of type stan::math::var and copy the contents of theta_ to it
   for (int i = 0; i < theta.size(); i++)
     { theta[i] = theta_[i];}

   VectorXd grad(theta_.size()); // value of the gradient of likelihood at theta

   // now we need to compute the likelihood using all the given data, save it as type stan::math::var, and then call the .grad() method to compute the gradient


   // a quick test:
   // this compiles without issue
   stan::math::var x = theta[0];
   stan::math::var y = x + 10;
   y.grad();





   return(grad);
 }
