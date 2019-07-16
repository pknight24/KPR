#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix computeErrorMatrixCpp(NumericMatrix Zrand,
                                 NumericMatrix Erand,
                                 NumericVector Yrand,
                                 NumericMatrix Hrand,
                                 NumericMatrix Q,
                                 NumericVector lambda,
                                 int K,
                                 bool cov_missing)
{
    unsigned int n = Zrand.nrow();
    unsigned int p = Zrand.ncol();
    int n_lambda = lambda.size();
    IntegerVector allIdx = Range(0, n - 1);

    arma::mat Zrand_arma = as<arma::mat>(Zrand);
    arma::mat Erand_arma = as<arma::mat>(Erand);
    arma::vec Yrand_arma = as<arma::vec>(Yrand);
    arma::mat Hrand_arma = as<arma::mat>(Hrand);
    arma::mat Q_arma = as<arma::mat>(Q);


    arma::mat errors(K, n_lambda);

    for (int j = 0; j < n_lambda; j++)
    {
      for (int k = 0; k < K; k++)
      {
        // the test index isn't perfect; may want to look at this
        IntegerVector testIdx = Range(n / K * k, n / K * (k + 1) - 1);
        IntegerVector trainIdx = setdiff(allIdx, testIdx);
        arma::uvec testIdx_arma = as<arma::uvec>(testIdx);
        arma::uvec trainIdx_arma = as<arma::uvec>(trainIdx);

        arma::vec Ytrain = Yrand_arma.elem(trainIdx_arma);
        arma::vec Ytest =  Yrand_arma.elem(testIdx_arma);
        arma::mat Ztrain = Zrand_arma.rows(trainIdx_arma);
        arma::mat Ztest = Zrand_arma.rows(testIdx_arma);
        arma::mat Etrain = Erand_arma.rows(trainIdx_arma);
        arma::mat Etest = Erand_arma.rows(testIdx_arma);
        arma::mat Htrain = Hrand_arma.submat(trainIdx_arma, trainIdx_arma);
        arma::mat Htest = Hrand_arma.submat(testIdx_arma, testIdx_arma);

        int n_train = trainIdx.size();
        arma::mat P = arma::eye(n_train,n_train);
        if(!cov_missing)
        {
          P = P - (Etrain * arma::inv(Etrain.t() * Htrain * Etrain) * Etrain.t() * Htrain);
        }

        arma::mat Y_p = P * Ytrain;
        arma::mat Z_p = P * Ztrain;

        arma::vec beta_hat = Q_arma * arma::inv(Z_p.t() * Htrain * Z_p * Q_arma + lambda[j] * arma::eye(p,p)) * Z_p.t() * Htrain * Y_p;
        arma::vec eta_hat = arma::zeros<arma::vec>(1);
        if (!cov_missing) eta_hat = arma::inv(Etrain.t() * Htrain * Etrain) * Etrain.t() * Htrain * (Ytrain - Ztrain * beta_hat);

        arma::vec coefficients = arma::join_vert(beta_hat, eta_hat);
        arma::mat Xtest = join_horiz(Ztest, Etest);

        arma::vec yhat = Xtest * coefficients;
        errors(k, j) = arma::as_scalar((Ytest - yhat).t() * Htest * (Ytest - yhat));

      }
    }

    return wrap(errors);
}
