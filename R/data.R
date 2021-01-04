#' Flaxseed lignan extract intervention microbiome data
#'
#' Baseline microbiome data (X) from a study on the effects of flaxseed lignan extract intervention on enterolactone production (Y).
#' The H matrix is a similarity kernel generated from unweighted Unifrac distances between the samples. The Q matrix is a similarity kernel generated from an Aitchisen variation matrix computed between the p variables.
#' The Y vector is already centered and scaled, and the X matrix has already been normalized with the centered log ratio and then column-centered and scaled.
#' @usage data(flax)
#' @format A list of the four matrices descibed above.
#' @references Lampe et al. (2019) The American Journal of Clinical Nutrition
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/31175806}{PubMed})
#' @examples
#' data(flax)
#' KPR(X = flax$X, Y = flax$Y, Q = flax$Q, H = flax$H)
"flax"


#' Age-related microbiome data
#'
#' Raw microbiome abundance data (n = 100, p = 149 genera) from the geographical study by Yatsunenko et al (2012).
#' @usage data(yatsunenko)
#' @format A list with the following entries:
#' \describe{
#'   \item{\code{raw.counts}}{The raw abundance counts for each sample, across 149 genera.}
#'   \item{\code{age}}{The age (in years) of each subject.}
#'   \item{\code{patristic}}{Patristic differences between each of the 149 genera.}
#'   \item{\code{unifrac}}{Unweighted UniFrac distances bewteen each of the 100 samples.}
#'   \item{\code{ec}}{1378 Enzyme commission (ec) summaries of bacterial genome data for each sample; see the paper cited below for reference.}
#'   \item{\code{geography}}{Country of origin for each sample.}
#' }
#' @references Yatsunenko et al. (2012) Nature (\href{https://www.ncbi.nlm.nih.gov/pubmed/22699611}{PubMed})
#'
#' See the KPR vignette for an in-depth example.
"yatsunenko"
