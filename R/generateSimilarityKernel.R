#' Generate a similarity kernel from a distance matrix
#'
#' Computes Gower's centered similarity matrix from a distance matrix.
#' @param D A square distance matrix. This may also be an object of class "dist".
#' @return A similarity kernel matrix of the same dimension as the input matrix.
#' @export
generateSimilarityKernel <- function(D)
{
  if (class(D) == "dist") D <- as.matrix(D)
  if (dim(D)[1] != dim(D)[2]) stop("D is not a square matrix")

  n <- dim(D)[1]
  J <- diag(n) - 1 / n
  M <- -J %*% D %*% J / 2

  eigen.M <- eigen(M)
  if (all(eigen.M$values >= 0 )) return(M)
  else
  {
    cat("\nSetting negative eigenvalues to 0\n")
    eigen.M$values[eigen.M$values < 0 ] <- 0
    M <- eigen.M$vectors %*% diag(eigen.M$values) %*% t(eigen.M$vectors)
    return(M)
  }

}
