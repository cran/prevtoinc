% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateIfunctions.R
\name{calculate_I_mixed}
\alias{calculate_I_mixed}
\title{Function to calculate incidence from PPS data using a mix of two estimators}
\usage{
calculate_I_mixed(I.pps.1, I.pps.2, a = 0.01, b = 500,
  method = "pps.mixed")
}
\arguments{
\item{I.pps.1}{resulting data frame for first estimator}

\item{I.pps.2}{resulting data frame for second estimator}

\item{a}{parameter a for the sigmoid function}

\item{b}{parameter b for the sigmoid function}

\item{method}{name of the method}
}
\value{
one-row data frame with following columns \itemize{
   \item{n - }{number of patients sampled} 
   \item{n.noso - }{number of HAIs}
   \item{P.hat - }{estimate of prevalence P}
   \item{I.hat - }{estimate of incidence rate I}
   \item{I.pp.hat - }{estimate of incidence proportion per admission I.pp}
   \item{x.loi.hat - }{estimate of x.loi}
   \item{x.los.hat - }{estimate of x.los}
   \item{method - }{name of the method}
}
}
\description{
A sigmoid function with parameters a and b (see below) is used to get weights
for a combination of the two estimator for x.loi and x.los.
}
\details{
is achieved in the following way for estimation of x.loi
alpha = exp(a*(n.noso-b))/(1+exp(a*(n.noso-b)))
x.loi.hat.mixed = alpha*x.loi.hat.1 + (1-alpha)*x.loi.hat.2

alpha = exp(a*(n-b))/(1+exp(a*(n-b)))
x.los.hat.mixed = alpha*x.los.hat.1 + (1-alpha)*x.los.hat.2
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

# estimate of incidence
I.1 <- calculate_I_smooth(data = data.pps.fast,
                   method = "gren")
                   
# estimate incidence based on Rhame-Sudderth formula
I.2 <- calculate_I_rhame(data = data.pps.fast,
                   x.loi.hat = 8,
                   x.los.hat = 13)

# mixed estimator                                      
calculate_I_mixed(I.1, I.2)

}
