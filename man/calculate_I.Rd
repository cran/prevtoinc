% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulationsfunctions.R
\name{calculate_I}
\alias{calculate_I}
\title{Function to calculate different estimators for I from PPS data.}
\usage{
calculate_I(data, data.theo = NULL)
}
\arguments{
\item{data}{data frame as generated by `simulate_pps_data` or `simulate_pps_fast`}

\item{data.theo}{data frame as generated by `simulate_incidence_stats` or `simulate incidence_stats_fast``}
}
\value{
data frame with following columns \itemize{
   \item{n - number of patients sampled} 
   \item{n.noso - number of HAIs}
   \item{P.hat - estimate of prevalence P}
   \item{I.hat - estimate of incidence rate I}
   \item{I.pp.hat - estimate of incidence proportion per admission I.pp}
   \item{x.loi.hat - estimate of x.loi}
   \item{x.los.hat - estimate of x.los}
   \item{method - name of the method}
}
and rows for the estimators gren, rear, pps.median, pps.mean, pps.mixed, 
rhame.theo, L.full 
(for a description of the estimators see vignette).
}
\description{
Function to calculate different estimators for I from PPS data.
}
\examples{

# create example data for PPS
example.dist <- create_dist_vec(function(x) dpois(x-1, 7), max.dist = 70)
example.dist.los <- create_dist_vec(function(x) dpois(x-1, lambda = 12),
                                    max.dist = 70)
data.pps.fast <- simulate_pps_fast(n.sample=200,
                                   P=0.05,
                                   dist.X.loi = example.dist,
                                   dist.X.los = example.dist.los)
head(data.pps.fast)

# estimate incidence based on Rhame-Sudderth formula
calculate_I(data = data.pps.fast)

}
