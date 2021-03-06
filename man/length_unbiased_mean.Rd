% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimationfunctions.R
\name{length_unbiased_mean}
\alias{length_unbiased_mean}
\title{Calculate mean of length-unbiased distribution from discrete length-biased distribution 
starting at 1}
\usage{
length_unbiased_mean(dist)
}
\arguments{
\item{dist}{vector of probabilities of length-biased distribution}
}
\value{
mean of length-unbiased distribution
}
\description{
Calculate mean of length-unbiased distribution from discrete length-biased distribution 
starting at 1
}
\examples{

# geometric distribution starting in 1 and cutoff at 70 with mean at about 8.
geom.dist <- create_dist_vec(geom_dist_fct, max.dist = 70)
# calculate mean of distribution
sum(1:length(geom.dist)*geom.dist)
# create length-biased distribution in same format
geom.dist.lb <- length_biased_dist(geom.dist)

# recalculate mean of original distribution based on length-biased distribution
length_unbiased_mean(geom.dist.lb)

}
