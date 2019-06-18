globalVariables(c("n",
                  "n.noso", 
                  "P.hat",
                  "x.los.hat",
                  "x.loi.hat"))



#' Calculate confidence intervals for Grenander estimator of Ipp based on asymptotics
#' 
#' Asymptotic or bootstrap approximation of confidence intervals for estimates of Ipp with gren method
#' Can use output of calculate_I_* functions as input. The asymptotic method uses 
#' the asymptotic normality of the estimator of I.pp to calculate the confidence interval and
#' the method "bs" uses a parametric bootstrap approximation based on the "naive" estimator. 
#'
#' @param data data frame which contains at least the following columns
#'  \itemize{
#'    \item{n - }{number of patients sampled} 
#'    \item{n.noso - }{number of HAIs}
#'    \item{P.hat - }{estimate of prevalence P}
#'    \item{x.loi.hat - }{estimate of x.loi}
#'    \item{x.los.hat - }{estimate of x.los}
#'    \item{I.pp.hat - }{estimate of incidence proportion per admission I.pp}
#' }
#' @param method either "asymptotic" for asymptotic confidence interval or "bs" for bootstrap-based confidence interval
#' @param alpha confidence level
#' @param n_bs number of bootstrap replications if method is "bs"
#'
#' @return tibble with columns CI.lower.Ipp and CI.upper.Ipp
#' 
#' @examples 
#' 
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
#'                                     max.dist = 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=5000,
#'                                    P=0.05,
#'                                    dist.X.loi = example.dist,
#'                                    dist.X.los = example.dist.los)
#' gren_est <- calculate_I_smooth(data = data.pps.fast, method = "gren")
#' gren_est
#' calculate_CI_I_pp(gren_est, method = "asymptotic", alpha = 0.05)
#' 
#' @export
calculate_CI_I_pp <- function(data, method = "asymptotic", alpha = 0.05, n_bs = 10000) {
  if( method == "asymptotic") {
    data <- data %>% 
      dplyr::mutate(sd = sqrt(
        P.hat^2*x.los.hat^2/x.loi.hat^2*(
          1/n     *    (1-P.hat)/P.hat +
            1/n     *  (1-1/x.los.hat)*x.los.hat +
            1/n.noso * (1-1/x.loi.hat)*x.loi.hat)
      ))
      tibble::tibble(CI.lower.Ipp = data$I.pp.hat+stats::qnorm(alpha/2)*data$sd,
        CI.upper.Ipp = data$I.pp.hat+stats::qnorm(1-alpha/2)*data$sd)
  } else {
    if (method == "bs") {
      
      
      data_CI <- purrr::pmap_df(data,
                         function(n, n.noso, P.hat, x.los.hat, x.loi.hat,  ...) {
                           CI <- stats::quantile((stats::rbinom(n_bs, n, P.hat)/n)/
                                            (stats::rbinom(n_bs , n, 1/x.los.hat)/n)*
                                            (stats::rbinom(n_bs , n.noso, 1/x.loi.hat)/n.noso),
                                          c(alpha/2, 1- alpha/2))
                           
                           tibble::tibble(CI.lower.Ipp = CI[1],
                                  CI.upper.Ipp = CI[2])
                         } )
      data_CI
      
    } else stop("Method for calculating CIs unknown.")
  } 
}


#' Function to calculate confidence intervals I.pp for gren estimator with bootstrap method based on Grenander estimator
#' 
#' Implements a bootstrap procedure for estimation of confidence intervals for I.pp based on boostrapping 
#' from the length of stay/infection distributions estimated by the gren method.
#' 
#' @param data data frame which contains a column `A.loi` 
#' with lengths of nosocomial infections up to survey point ( zero if none) and
#'  a column `A.los` with length of stay up to survey point
#' @param n_bs number of bootstrap samples to use for calculations
#' @param alpha confidence level
#' 
#' @return single-row tibble with columns CI.lower.Ipp and CI.upper.Ipp
#' 
#' @examples 
#' 
#' example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
#' example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
#'                                     max.dist = 70)
#' data.pps.fast <- simulate_pps_fast(n.sample=5000,
#'                                    P=0.05,
#'                                    dist.X.loi = example.dist,
#'                                    dist.X.los = example.dist.los)
#' gren_est <- calculate_I_smooth(data = data.pps.fast, method = "gren")
#' gren_est
#' CI_np_bs(data.pps.fast, n_bs = 500)
#' @export
CI_np_bs <- function(data, n_bs = 1000, alpha = 0.05) {
  # create distribution to sample from
  noso.vec <- na.omit(data$A.loi)
  n <- length(noso.vec)
  n.noso <- length(noso.vec[noso.vec>0])
  P.hat <- sum(noso.vec > 0)/length(noso.vec)
  
  p.A.loi.bs <-
    monotone_smoother(
      noso.vec[noso.vec > 0], method = "gren")
  p.A.los.bs <-
    monotone_smoother(
      na.omit(data$A.los), method = "gren")
  
  p.A.loi.sample <-
    replicate(n_bs,
              monotone_smoother(sample.int(length(p.A.loi.bs),
                                           n.noso, replace = TRUE, p.A.loi.bs),
                                method = "gren")[1])
  
  p.A.los.sample <-
    replicate(n_bs,
              monotone_smoother(sample.int(length(p.A.los.bs),
                                           n, replace = TRUE, p.A.los.bs),
                                method = "gren")[1])
  
  P.sample <- stats::rbinom(n_bs, n, prob = P.hat)/n
  CI <- stats::quantile(P.sample*p.A.loi.sample/p.A.los.sample, c(alpha/2,1-alpha/2))
  tibble::tibble(CI.lower.Ipp = CI[1],
         CI.upper.Ipp = CI[2])

}
