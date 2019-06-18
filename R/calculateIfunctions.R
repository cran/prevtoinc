#' Estimate the incidence based on PPS data using monotone estimators for the distribution of A.
#' 
#' Estimate incidence from PPS by the method proposed in the companion paper.
#' `data` should be supplied as a data frame with at least a column named `A.loi`
#' giving lengths of infection up to date of PPS.
#' Values of zero for `A.loi` indicate absence of a HAI. Optionally,
#' the data frame can also contain a column `A.los` supplying lengths
#' of stay up to PPS to estimate x.los with the same method as well.
#' If `correct.one` is `TRUE`, the number infections on their first day
#' will be augmented to be at least as high as the number of infections
#' on their second day for the estimation of x.loi .
#'
#' @param data data frame which contains a column `A.loi` 
#' with lengths of nosocomial infections up to survey point ( zero if none) and
#' possibly a column `A.los` with length of stay up to survey point
#' @param method method to use for smoothing ("gren" ( Grenander ) or "rear" (rearrangement))
#' @param correct.one.loi use correction for underreporting of one day LOIs:
#'  "no" if none, "fill.ones" to set the one-day cases to be at least the number
#'   of two-day cases, "start.two" to only use P(A=2| A > 1) as a proxy for P(A=1)
#' @param correct.one.los use correction for underreporting of one day LOSs:
#'  "no" if none, "fill.ones" to set the one-day cases to be at least the number
#'   of two-day cases, "start.two" to only use P(A=2| A > 1) as a proxy for P(A=1)
#'
#' @return one-row data frame with following columns \itemize{
#'    \item{n - }{number of patients sampled} 
#'    \item{n.noso - }{number of HAIs}
#'    \item{P.hat - }{estimate of prevalence P}
#'    \item{I.hat - }{estimate of incidence rate I}
#'    \item{I.pp.hat - }{estimate of incidence proportion per admission I.pp}
#'    \item{x.loi.hat - }{estimate of x.loi}
#'    \item{x.los.hat - }{estimate of x.los}
#'    \item{method - }{name of the method}
#' }
#' 
#' @examples
#' 
#' # create example data for PPS
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
#'                                     max.dist = 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=200,
#'                                    P=0.05,
#'                                    dist.X.loi = example.dist,
#'                                    dist.X.los = example.dist.los)
#' head(data.pps.fast)
#' 
#' # estimate of incidence
#' calculate_I_smooth(data = data.pps.fast,
#'                    method = "gren")
#' 
#' 
#' @export
calculate_I_smooth <- function(data,
                               method = "gren",
                               correct.one.loi="no",
                               correct.one.los="no") {
  # vectors of length of stay up to date of survey and length
  # of nosocomial infection up to survey
  noso.vec <- na.omit(data$A.loi)
  n <- length(noso.vec)

  P.hat <- sum(noso.vec > 0)/length(noso.vec)
  P.odds <- P.hat/(1 - P.hat)
  if (length(noso.vec[noso.vec>0]) == 0) {
    I.hat <- 0; x.loi.hat <- NA
  } else {
    if (correct.one.loi == "fill.ones") { 
      noso.vec <- c(noso.vec,
                    rep(1, max(0,sum(noso.vec==2) -
                                 sum(noso.vec==1))))
    } else if (correct.one.loi == "start.two") {
      noso.vec <- noso.vec - 1
      noso.vec[noso.vec == -1] <- 0
    }
    
    n.noso <- length(noso.vec[noso.vec>0])
    
    p.A <-
      monotone_smoother(
        noso.vec[noso.vec > 0], method = method)
    x.loi.hat <- 1/p.A[1]
    if (correct.one.loi == "start.two") x.loi.hat <- x.loi.hat + 1
    I.hat <- P.odds/x.loi.hat
  }
  
  if ("A.los" %in% colnames(data) && any(!is.na(data$A.los))) {
    los.vec <- na.omit(data$A.los)
    if (correct.one.los == "fill.ones") { 
      los.vec <- c(los.vec,
                   rep(1, max(0,sum(los.vec==2) -
                                sum(los.vec==1))))
    } else if (correct.one.los == "start.two") {
      los.vec <- los.vec - 1
      los.vec[los.vec == -1] <- 0
    }
    
    p.A.los <-
      monotone_smoother(
        los.vec[los.vec > 0],
        method = method)
    x.los.hat <- 1/p.A.los[1]
    if (correct.one.los == "start.two") x.los.hat <- x.los.hat + 1
  } else {
    x.los.hat <- NA
  }
  I.pp.hat <- I.hat*(1-P.hat)*x.los.hat
  
  return(tibble::tibble(
    n = n,
    n.noso = n.noso,
    P.hat = P.hat,
    I.hat = I.hat,
    I.pp.hat = I.pp.hat,
    x.loi.hat = x.loi.hat,
    x.los.hat = x.los.hat,
    method = method))
}


#' Function to calculate incidence from PPS data using a Rhame-Sudderth like approach 
#' with estimates for x.loi and x.los supplied.
#'
#'
#' @param data one-row data frame which contains a column A.loi (only used to calculate P.hat) 
#' with lengths of nosocomial infections up to survey (a 0 indicates no HAI present)
#' @param x.loi.hat value for estimated expected length of infection x_loi
#' @param x.los.hat value for estimated expected length of stay x_los (optional)
#' @param method a string with associated name for method
#'
#' @return one-row data frame with following columns \itemize{
#'    \item{n - }{number of patients sampled} 
#'    \item{n.noso - }{number of HAIs}
#'    \item{P.hat - }{estimate of prevalence P}
#'    \item{I.hat - }{estimate of incidence rate I}
#'    \item{I.pp.hat - }{estimate of incidence proportion per admission I.pp}
#'    \item{x.loi.hat - }{estimate of x.loi}
#'    \item{x.los.hat - }{estimate of x.los}
#'    \item{method - }{name of the method}
#' }
#' 
#' @examples
#' 
#' # create example data for PPS
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
#'                                     max.dist = 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=200,
#'                                    P=0.05,
#'                                    dist.X.loi = example.dist,
#'                                    dist.X.los = example.dist.los)
#' head(data.pps.fast)
#' 
#' # estimate incidence based on Rhame-Sudderth formula
#' calculate_I_rhame(data = data.pps.fast,
#'                    x.loi.hat = 8,
#'                    x.los.hat = 13)
#' 
#' @export
calculate_I_rhame <- function(data, 
                              x.loi.hat, 
                              x.los.hat = NA, 
                              method = "rhame") {
  # vectors of length of stay up to date of survey and length
  # of nosocomial infection up to survey
  n = nrow(data)
  noso.vec <- na.omit(data$A.loi)
  n.noso <- length(noso.vec[noso.vec>0])
  
  P.hat <- sum(noso.vec > 0)/length(noso.vec)
  
  I.hat <- P.hat/(1-P.hat)/x.loi.hat
  I.pp.hat <- I.hat*(1-P.hat)*x.los.hat
  
  return(tibble::tibble(
    n = n,
    n.noso = n.noso,
    P.hat = P.hat,
    I.hat = I.hat,
    I.pp.hat = I.pp.hat,
    x.loi.hat = x.loi.hat,
    x.los.hat = x.los.hat,
    method = method)
  )
}


#' Function to calculate incidence from PPS data using a mix of two estimators 
#' 
#' A sigmoid function with parameters a and b (see below) is used to get weights
#' for a combination of the two estimator for x.loi and x.los.
#'
#'  is achieved in the following way for estimation of x.loi
#' alpha = exp(a*(n.noso-b))/(1+exp(a*(n.noso-b)))
#' x.loi.hat.mixed = alpha*x.loi.hat.1 + (1-alpha)*x.loi.hat.2
#' 
#' alpha = exp(a*(n-b))/(1+exp(a*(n-b)))
#' x.los.hat.mixed = alpha*x.los.hat.1 + (1-alpha)*x.los.hat.2
#'
#' @param I.pps.1 resulting data frame for first estimator
#' @param I.pps.2 resulting data frame for second estimator
#' @param a parameter a for the sigmoid function 
#' @param b parameter b for the sigmoid function
#' @param method name of the method
#'
#' @return one-row data frame with following columns \itemize{
#'    \item{n - }{number of patients sampled} 
#'    \item{n.noso - }{number of HAIs}
#'    \item{P.hat - }{estimate of prevalence P}
#'    \item{I.hat - }{estimate of incidence rate I}
#'    \item{I.pp.hat - }{estimate of incidence proportion per admission I.pp}
#'    \item{x.loi.hat - }{estimate of x.loi}
#'    \item{x.los.hat - }{estimate of x.los}
#'    \item{method - }{name of the method}
#' }
#' 
#' @examples 
#' 
#' # create example data for PPS
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
#'                                     max.dist = 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=200,
#'                                    P=0.05,
#'                                    dist.X.loi = example.dist,
#'                                    dist.X.los = example.dist.los)
#' head(data.pps.fast)
#' 
#' # estimate of incidence
#' I.1 <- calculate_I_smooth(data = data.pps.fast,
#'                    method = "gren")
#'                    
#' # estimate incidence based on Rhame-Sudderth formula
#' I.2 <- calculate_I_rhame(data = data.pps.fast,
#'                    x.loi.hat = 8,
#'                    x.los.hat = 13)
#' 
#' # mixed estimator                                      
#' calculate_I_mixed(I.1, I.2)
#' 
#' @export
calculate_I_mixed <- function(I.pps.1, I.pps.2, a = 0.01, b = 500, 
                              method = "pps.mixed") {
  method.name <- method
  alpha.loi <- exp(a*(I.pps.1$n.noso-b))/(1+exp(a*(I.pps.1$n.noso-b)))
  alpha.los <- exp(a*(I.pps.1$n-b))/(1+exp(a*(I.pps.1$n-b)))
  I.new.mixed <- I.pps.1 %>%
    dplyr::mutate(x.los.hat = (1-alpha.los)*.data$x.los.hat+
                    alpha.los*I.pps.2$x.los.hat,
                  x.loi.hat = (1-alpha.loi)*.data$x.loi.hat+
                    alpha.loi*I.pps.2$x.loi.hat) %>%
    dplyr::mutate(
      I.hat = .data$P.hat/(1-.data$P.hat)/.data$x.loi.hat,
      I.pp.hat = .data$P.hat/.data$x.loi.hat*.data$x.los.hat) %>%
    dplyr::mutate(method = method.name)
  I.new.mixed
}