#' Generate a similarity kernel from a distance matrix
#'
#' Computes Gower's centered similarity matrix from a distance matrix.
#' @param D A square distance matrix. This may also be an object of class "dist".
#' @param squareValues Logical, indicating whether to square the values in `D` before building the similarity kernel.
#' @return A similarity kernel matrix of the same dimension as the input matrix.
#' @export
generateSimilarityKernel <- function(D, squareValues = TRUE)
{
  if (class(D) == "dist") D <- as.matrix(D)
  if (dim(D)[1] != dim(D)[2]) stop("D is not a square matrix")

  n <- dim(D)[1]
  J <- diag(n) - 1 / n
  if (squareValues) M <- -J %*% (D^2) %*% J / 2
  else M <- -J %*% D %*% J / 2

  eigen.M <- eigen(M)
  if (all(eigen.M$values >= 10^-5 )) return(M)
  else
  {
    cat("Correcting small and negative eigenvalues\n")
    smallest <- min(eigen.M$values[eigen.M$values > 10^-5])
    eigen.M$values[eigen.M$values < 10^-5 ] <- smallest / 2
    M <- eigen.M$vectors %*% diag(eigen.M$values) %*% t(eigen.M$vectors)
    return(M)
  }

}
