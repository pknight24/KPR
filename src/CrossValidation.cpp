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
    int n = Zrand.nrow();
    int p = Zrand.ncol();
    int n_lambda = lambda.size();

    NumericMatrix errors(K, n_lambda);
    for (int j = 0; j < n_lambda; j++)
    {
      for (int k = 0; k < K; k++)
      {
        // NumericVector
        // NumericVector Ytrain = Yrand[]

      }
    }

    return errors;
}
