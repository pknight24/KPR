#' Aitchisen's Variation Matrix
#'
#' Computes Aitchisen's variation matrix from a matrix X
#' @param X An n x p data matrix.
#' @return A p x p Aitchisen variation matrix.
#' @export
aitchisenVariation <- function(X)
{
  p <- ncol(X)
  A <- matrix(nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      A[i, j] = var(log(X[, i] / X[, j])) / sqrt(2)  #
    }
  }
  return(A)
}
