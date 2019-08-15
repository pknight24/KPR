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
  cat("Kernel Penalized Regression results, using the GMD inference.\n\n")


  infer.out <- object$p.vals

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


#' A "supervised" GMD biplot
#'
#' This plots a GMD biplot of an object of class \code{KPR}. The axes correspond to first and second
#' GMD components of the data, as described by Wang et al. \code{biplot.KPR}
#' performs inference on the \code{KPR} object, and plots arrows corresponding to
#' penalized variables deemed significant. The points represent samples and are
#' colored with respect to the response vector \code{Y}.
#'
#' @param x An object of class \code{KPR}
#' @param ... Additional parameters, used by the \code{inference} function call.
#' @importFrom viridis plasma
#' @importFrom graphics arrows axTicks axis plot points text legend
#' @importFrom stats quantile
#' @references Wang et al. Technical report.
#' @export
biplot.KPR <- function(x, ...)
{
  Z <- x$Z
  H <- x$H
  Q <- x$Q
  K <- 10
  gmd.out <- GMD(X = Z, H = H, Q = Q, K = K)
  U <- gmd.out$U
  S <- gmd.out$S
  V <- gmd.out$V

  S.order <- order(S, decreasing = TRUE)

  U = U[,S.order]
  V = V[,S.order]

  k1 = S.order[1]
  k2 = S.order[2]

  S = S[S.order]

  eta = U%*%diag(S)
  total.variation <- sum(S)
  percent.variation <- (S / total.variation) * 100


  max.xlab = max(abs(eta[,1]))
  max.ylab = max(abs(eta[,2]))

  ycolor <- x$Y
  order <- findInterval(ycolor, sort(ycolor))
  sample.col = plasma(length(order))[order]

  arrow.col = 'gray50'
  legend.col = 'black'
  sample.pch <- 19

  plot(eta[,1], eta[,2],
    xlab = paste0('GMD Component ',k1, " (",round(percent.variation[1], digits = 2),"% explained)"),
    ylab = paste0("GMD Component ",k2, " (",round(percent.variation[2], digits = 2),"% explained)"),
    pch = sample.pch,
    xlim = c(-1.1*max.xlab, 1.1*max.xlab ), ylim =  c(-1.1*max.ylab, 1.1*max.ylab) ,
    col = sample.col)
  xaxp = axTicks(1)
  yaxp = axTicks(2)

  if (length(x$lambda) > 1) p.values <- x$p.values[,1]
  else p.values <- x$p.values
  signif <- which(p.values < 0.05) # we keep track of the significant coefficients

  #calculate coordinates
  V.plot = Q%*%V
  arrow.x = V.plot[,1]
  arrow.y = V.plot[,2]

  arrow.mat <- cbind(arrow.x, arrow.y)

  norms <- apply(arrow.mat, 1, function(xx) norm(xx, type="2"))

  big.norms <- which(norms > quantile(norms, 0.25)) # we only want the arrows with an L2 norm past the .25 quantile

  index  <- intersect(signif, big.norms)

  if (is.null(colnames(x$Z))) names = paste0("V", index)
  else names <- colnames(x$Z)[S.order]
  iter = 1

  max.xarrow = max(abs(arrow.x))
  max.yarrow = max(abs(arrow.y))
  xratio = max.xarrow/max.xlab
  yratio = max.yarrow/max.ylab


  xsci = as.numeric(unlist(strsplit(formatC(xratio, format = 'e'),"e")))
  xlab.arrow = round(xaxp*xsci[1]*10^(xsci[2]), digits = 2)


  ysci = as.numeric(unlist(strsplit(formatC(yratio, format = 'e'),"e")))
  ylab.arrow = round(yaxp*ysci[1]*10^(ysci[2]), digits = 2)

  points(eta[,1], eta[,2], col = sample.col, pch = sample.pch)

  legend("topleft", title="Outcome (Y)",
         legend=c(round(min(x$Y), digits = 3), round(max(x$Y), digits = 3)),
         col=c(sample.col[which.min(x$Y)], sample.col[which.max(x$Y)]),
         inset=0.01, pch = sample.pch, cex = 0.75)

  for(i in index){


      arrows(x0 = 0,y0 = 0,x1 = arrow.x[i]/xratio, y1 = arrow.y[i]/yratio, length = 0.05, col = arrow.col)
      text(arrow.x[i]/xratio, arrow.y[i]/yratio*1.1, names[iter], cex = 1, col = legend.col)

    iter = iter + 1
  }

  # add new axis
  axis(3, at = xaxp, labels = as.character(xlab.arrow))
  axis(4, at = yaxp, labels = as.character(ylab.arrow))


}
