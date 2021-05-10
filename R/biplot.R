#' A "supervised" GMD biplot
#'
#' This plots a GMD biplot where the axes correspond to first and second
#' GMD components of the data, as described by Wang et al. 
#'
#' @param X An n by p data matrix.
#' @param Y An optional outcome vector. If provided, the points in the biplot will be colored according to the corresponding value in Y.
#' @param H An n by n sample-wise similarity kernel.
#' @param Q A p by p variable-wise similarity kernel.
#' @param K An integer, corresponding to the number of GMD components to compute.
#' @importFrom viridis plasma
#' @importFrom graphics arrows axTicks axis plot points text legend
#' @importFrom stats quantile
#' @references Wang et al. Technical report.
#' @export
GMD.biplot <- function(X, Y = NULL, H = diag(nrow(X)), Q = diag(ncol(X)), K = 10)
{
  gmd.out <- GMD(X = X, H = H, Q = Q, K = K)
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

  ycolor <- Y
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
  

  # if (is.null(x$p.values)) signif <- rep(1, ncol(X)) # if there are no pvalues in the object, use all the coefficients
  # else signif <- which(x$p.values < 0.05) # we keep track of the significant coefficients

  #calculate coordinates
  V.plot = Q%*%V
  arrow.x = V.plot[,1]
  arrow.y = V.plot[,2]

  arrow.mat <- cbind(arrow.x, arrow.y)

  norms <- apply(arrow.mat, 1, function(xx) norm(xx, type="2"))

  big.norms <- which(norms > quantile(norms, 0.25)) # we only want the arrows with an L2 norm past the .25 quantile

  # index  <- intersect(signif, big.norms)
  index <- big.norms

  if (is.null(colnames(X))) names = paste0("V", index)
  else names <- colnames(X)[index]
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
  
  
  if (!is.null(Y))
    legend("topleft", title="Outcome (Y)",
         legend=c(round(min(Y), digits = 3), round(max(Y), digits = 3)),
         col=c(sample.col[which.min(Y)], sample.col[which.max(Y)]),
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
