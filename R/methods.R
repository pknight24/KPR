#' @export
print.KPR <- function(x, ...)
{
  cat("A fitted KPR model of the form\n\n")
  cat("\t y = Zb + ")
  if (!is.null(x$E)) cat("Eh + ")
  cat("e\n\n")
  q.trivial <- ifelse(all(x$Q == diag(ncol(x$Z))), "trivial","nontrivial" )
  h.trivial <- ifelse(all(x$H == diag(nrow(x$Z))), "trivial","nontrivial" )
  cat("with a", h.trivial, "H matrix and a", q.trivial, "Q matrix.\n\n")
  if (x$REML)
  {
    cat("Lambda found with REML estimation: ", x$lambda, "\n")
  }
  else
  {
    if (length(x$lambda) == 1) cat("Given value of lambda: ", x$lambda, "\n")
    else
    {
        cat("Lambda values evaluated with", nrow(x$cv.errors), "fold cross validation.\n")
        cat("    lambda.min:", x$lambda.min, "at index",x$lambda.min.index,"of", length(x$lambda),"\n")
        cat("    lambda.1se:", x$lambda.1se,"at index",x$lambda.1se.index,"of", length(x$lambda),"\n")
    }
  }

}

#' @export
summary.KPR <- function(object, ...)
{
  if (missing(inference.method)) inference.method <- "GMD"

  cat("Kernel Penalized Regression results, using", inference.method, "inference.\n\n")


  infer.out <- inference(object, method = inference.method)

  if (length(object$lambda) == 1)
  {
    b <- object$beta.hat
    l <- c(round(object$lambda), rep("", length(b) - 1))
    p <- infer.out

    data.frame(lambda = l, betahat = b, pvalue = p)

  }
  else
  {
    b <- object$beta.hat[,object$lambda.min.index]
    l <- c(round(object$lambda.min), rep("", length(b) - 1))
    p <- infer.out[,object$lambda.min.index]

    data.frame(lambda = l, betahat = b, pvalue = p)
  }

}
