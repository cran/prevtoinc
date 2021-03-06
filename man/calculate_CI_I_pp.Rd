% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CIfunctions.R
\name{calculate_CI_I_pp}
\alias{calculate_CI_I_pp}
\title{Calculate confidence intervals for Grenander estimator of Ipp based on asymptotics}
\usage{
calculate_CI_I_pp(data, method = "asymptotic", alpha = 0.05,
  n_bs = 10000)
}
\arguments{
\item{data}{data frame which contains at least the following columns
 \itemize{
   \item{n - }{number of patients sampled} 
   \item{n.noso - }{number of HAIs}
   \item{P.hat - }{estimate of prevalence P}
   \item{x.loi.hat - }{estimate of x.loi}
   \item{x.los.hat - }{estimate of x.los}
   \item{I.pp.hat - }{estimate of incidence proportion per admission I.pp}
}}

\item{method}{either "asymptotic" for asymptotic confidence interval or "bs" for bootstrap-based confidence interval}

\item{alpha}{confidence level}

\item{n_bs}{number of bootstrap replications if method is "bs"}
}
\value{
tibble with columns CI.lower.Ipp and CI.upper.Ipp
}
\description{
Asymptotic or bootstrap approximation of confidence intervals for estimates of Ipp with gren method
Can use output of calculate_I_* functions as input. The asymptotic method uses 
the asymptotic normality of the estimator of I.pp to calculate the confidence interval and
the method "bs" uses a parametric bootstrap approximation based on the "naive" estimator.
}
\examples{

example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
                                    max.dist = 70)
data.pps.fast <- simulate_pps_fast(n.sample=5000,
                                   P=0.05,
                                   dist.X.loi = example.dist,
                                   dist.X.los = example.dist.los)
gren_est <- calculate_I_smooth(data = data.pps.fast, method = "gren")
gren_est
calculate_CI_I_pp(gren_est, method = "asymptotic", alpha = 0.05)

}
