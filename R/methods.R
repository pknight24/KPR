#' @export
print.KPR <- function(kpr.out)
{
  cat("A fitted KPR model of the form\n\n")
  cat("\t y = Zb + ")
  if (!is.null(kpr.out$E)) cat("Eh + ")
  cat("e\n\n")
  q.trivial <- ifelse(all(kpr.out$Q == diag(ncol(kpr.out$Z))), "trivial","nontrivial" )
  h.trivial <- ifelse(all(kpr.out$H == diag(nrow(kpr.out$Z))), "trivial","nontrivial" )
  cat("with a", h.trivial, "H matrix and a", q.trivial, "Q matrix.\n\n")
  if (kpr.out$REML)
  {
    cat("Lambda found with REML estimation: ", kpr.out$lambda, "\n")
  }
  else
  {
    if (length(kpr.out$lambda) == 1) cat("Given value of lambda: ", kpr.out$lambda, "\n")
    else
    {
        cat("Lambda values evaluated with", nrow(kpr.out$cv.errors), "fold cross validation.\n")
        cat("    lambda.min:", kpr.out$lambda.min, "\n")
        cat("    lambda.1se:", kpr.out$lambda.1se,"\n")
    }
  }

}
