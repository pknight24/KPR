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
  if (!exists("inference.method")) inference.method <- "GMD"

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

#' @importFrom viridis viridis
#' @export
biplot.KPR <- function(x, ...)
{
  Z <- x$Z
  H <- x$H
  Q <- x$Q
  if (!exists("K")) K <- 10
  gmd.out <- GMD(X = Z, H = H, Q = Q, K = K)
  U <- gmd.out$U
  S <- gmd.out$S
  V <- gmd.out$V


  U = U[,order(S, decreasing = TRUE)]
  V = V[,order(S, decreasing = TRUE)]

  k1 = order(S, decreasing = TRUE)[1]
  k2 = order(S, decreasing = TRUE)[2]

  S = sort(S, decreasing = TRUE)

  eta = U%*%diag(S)


  max.xlab = max(abs(eta[,1]))
  max.ylab = max(abs(eta[,2]))

  ycolor <- x$Y
  order <- findInterval(ycolor, sort(ycolor))
  sample.col = viridis(length(order))[order]

  arrow.col = 'gray50'
  legend.col = 'black'
  sample.pch <- 19

  plot(eta[,1], eta[,2],
    xlab = paste0('PC',k1), ylab = paste0('PC',k2),
    pch = sample.pch,
    xlim = c(-1.1*max.xlab, 1.1*max.xlab ), ylim =  c(-1.1*max.ylab, 1.1*max.ylab) ,
    col = sample.col)
  xaxp = axTicks(1)
  yaxp = axTicks(2)

  infer.out <- inference(x)
  if (length(x$lambda) > 1) p.values <- infer.out[,x$lambda.min.index]
  else p.values <- infer.out
  index = which(p.values < 0.05)

  #calculate coordinates
  V.plot = Q%*%V
  arrow.x = V.plot[,1]
  arrow.y = V.plot[,2]

  if (is.null(colnames(x$Z))) names = paste0("V", index)
  else names <- colnames(x$Z)
  iter = 1

  max.xarrow = max(abs(arrow.x))
  max.yarrow = max(abs(arrow.y))
  xratio = max.xarrow/max.xlab
  yratio = max.yarrow/max.ylab


  xsci = as.numeric(unlist(strsplit(formatC(xratio, format = 'e'),"e")))
  xlab.arrow = round(xaxp*xsci[1]*10^(xsci[2]), digits = 2)


  ysci = as.numeric(unlist(strsplit(formatC(yratio, format = 'e'),"e")))
  ylab.arrow = round(yaxp*ysci[1]*10^(ysci[2]), digits = 2)

  for(i in index){


      arrows(x0 = 0,y0 = 0,x1 = arrow.x[i]/xratio, y1 = arrow.y[i]/yratio, length = 0.05, col = arrow.col)
      text(arrow.x[i]/xratio, arrow.y[i]/yratio*1.1, names[iter], cex = 1, col = legend.col)

    iter = iter + 1
  }

  # add new axis
  axis(3, at = xaxp, labels = as.character(xlab.arrow))
  axis(4, at = yaxp, labels = as.character(ylab.arrow))

  points(eta[,1], eta[,2], col = sample.col, pch = sample.pch)
}
