#' Flaxseed lignan extract intervention microbiome data
#'
#' Baseline microbiome data (X) from a study on the effects of flaxseed lignan extract intervention on enterolactone production (Y).
#' The H matrix is a similarity kernel generated from unweighted Unifrac distances between the samples. The Q matrix is a similarity kernel generated from an Aitchisen variation matrix computed between the p variables.
#' The Y vector is already centered, and the X matrix has already been normalized with the centered log ratio and then column-centered.
#' @usage data(flax)
#' @format A list of the four matrices descibed above.
#' @references Lampe et al. (2019) The American Journal of Clinical Nutrition
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/31175806}{PubMed})
#' @examples
#' data(flax)
#' KPR(X = flax$X, Y = flax$Y, Q = flax$Q, H = flax$H)
"flax"
