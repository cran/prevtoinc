#' Calculate mean of length-unbiased distribution from discrete length-biased distribution 
#' starting at 1
#' @param dist vector of probabilities of length-biased distribution
#' @return mean of length-unbiased distribution 

#' 
#' @examples 
#' 
#' # geometric distribution starting in 1 and cutoff at 70 with mean at about 8.
#' geom.dist <- create_dist_vec(geom_dist_fct, max.dist = 70)
#' # calculate mean of distribution
#' sum(1:length(geom.dist)*geom.dist)
#' # create length-biased distribution in same format
#' geom.dist.lb <- length_biased_dist(geom.dist)
#' 
#' # recalculate mean of original distribution based on length-biased distribution
#' length_unbiased_mean(geom.dist.lb)
#' 
#' @export
length_unbiased_mean <- function(dist) {
  1/sum(1/1:length(dist)*dist)
}

#' Calculate length-biased distribution from discrete length-unbiased distribution
#' starting at 1
#' @param dist vector of probabilities of distribution to transform
#' @return vector of probabilities of transformed distribution

#' 
#' @examples 
#' 
#' # geometric distribution starting in 1 and cutoff at 70 with mean at about 8.
#' geom.dist <- create_dist_vec(geom_dist_fct, max.dist = 70)
#' # calculate mean
#' sum(1:length(geom.dist)*geom.dist)
#' # plot original distribution
#' plot(geom.dist)
#' geom.dist.lb <- length_biased_dist(geom.dist)
#' # plot length biased distribution
#' plot(geom.dist.lb)
#' 
#' @export
length_biased_dist <- function(dist) {
  X.mean <- sum(1:length(dist)*dist)
  dist*1:length(dist)/X.mean
}

#' Calculate empirical probability mass function
#' for a discrete positive distribution starting at 1
#' 
#' @param values used for the calculation of the empirical pmf
#' @return vector of probabilities for epmf for the range 1:length(values)
#' 
#' @examples
#' 
#' # generate random sample of independent draws from Poisson distribution
#' x <- rpois(200,4)
#' # calculate empirical probability mass function and true probability mass function
#' y.emp  <- epmf(x)
#' y.theo <- dpois(1:max(x), 4)
#' plot(y.emp)
#' points(y.theo, col = "red")
#' 
#' @export
epmf <- function(values) {
  n <- length(values)
  sapply(1:max(values), function(i) sum(values==i)/n)
}

#' Transform a distribution of times of stay to a
#' distribution of staying-time up to observation point under
#' assumption of steady state.
#'
#' @param dist.stays vector of probabilities of being at the
#' hospital for 1:length(dist.stays) days at random time of
#' observation
#' @return vector of probabilities of staying
#' 1:length(dist.point) days
#' 
#' @examples 
#' 
#' # generate vector of probabilities for truncated Poisson distribution for 
#' # distribution of times of stay X
#' dist.X <- dpois(1:70, 4)
#' plot(dist.X)
#' # transform to distribution of distribution of staying-time up to observation point under
#' # assumption of steady state
#' dist.A <- X_to_A_dist(dist.X)
#' plot(dist.A)
#' # transform back to get original distribution
#' dist.X.2 <- A_to_X_dist(dist.A)
#' plot(dist.X.2)
#' 
#' @export
X_to_A_dist <- function(dist.stays) {
  numerators <- rev(cumsum(rev(dist.stays)))
  denominator <- sum(numerators)
  return(numerators/denominator)
}

#' function to transform the distribution of stays to a fixed
#' point to the distribution of the staying times
#'
#' @param dist.point vector of probabilities of staying
#' 1:length(dist.point) days
#'
#' @return  vector of probabilities
#' of being at the hospital for 1:length(dist.point) days at
#' random time of observation
#'
#' @examples 
#' 
#' # generate vector of probabilities for truncated Poisson distribution for 
#' # distribution of times of stay X
#' dist.X <- dpois(1:70, 4)
#' plot(dist.X)
#' # transform to distribution of distribution of staying-time up to observation point under
#' # assumption of steady state
#' dist.A <- X_to_A_dist(dist.X)
#' plot(dist.A)
#' # transform back to get original distribution
#' dist.X.2 <- A_to_X_dist(dist.A)
#' plot(dist.X.2)
#' 
#' @export
A_to_X_dist <- function(dist.point) {
  numerators <- dist.point - c(dist.point[2:length(dist.point)],
                               0)
  denominator <- sum(numerators)
  return(numerators/denominator)
}

#' Calculate a monotone probability mass function estimate
#' using a  rearrangement or a Grenander estimator as described in Jankoswski,
#' Wellner, 2009 <doi:10.1214/09-EJS526>
#'
#' @param values observed values of distribution
#' @param method method of estimation "rear" rearrangement or "gren" Grenander
#' @param range boundaries of the support of the distribution
#'
#' @return vector of estimated pmf (support of distribution is by default assumed to be min(values):max(values) )
#'
#' @examples
#' 
#' # generate sample from geometric distribution
#' A <- rgeom(50, 0.2)
#' # plot empirical probability mass function
#' plot(epmf(A))
#' dist.A.gren <- monotone_smoother(A, method = "gren")
#' # plot estimated probability mass function
#' points(dist.A.gren, col = "red")
#' 
#' @export
monotone_smoother <- function(values, method = "rear", 
                              range = c(1, max(values))) {
  p <- sapply(1:range[2], function(i) sum(values == i))
  p <- p/length(values)
  if (method == "gren") {
    return(.find_gren(p))
  } else if (method == "rear") {
    return(.find_rear(p)) 
    } else {
    stop("monotone smoother method not known")
  }
}

.find_rear <- function(p) {
  sort(p, decreasing = TRUE)
}


# R function to apply Grenander estimator
.find_gren <- function(p) {
  -stats::isoreg(1:length(p),-p)$yf
}

