#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix computeErrorMatrix(NumericMatrix Zrand,
                                 NumericMatrix Erand,
                                 NumericVector Yrand,
                                 NumericMatrix Hrand,
                                 NumericMatrix Q,
                                 NumericVector lambda,
                                 int K)
{
    Rcout << Zrand.ncol() << "\n";
    return Zrand;
}
