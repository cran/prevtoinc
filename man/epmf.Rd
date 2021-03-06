% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimationfunctions.R
\name{epmf}
\alias{epmf}
\title{Calculate empirical probability mass function
for a discrete positive distribution starting at 1}
\usage{
epmf(values)
}
\arguments{
\item{values}{used for the calculation of the empirical pmf}
}
\value{
vector of probabilities for epmf for the range 1:length(values)
}
\description{
Calculate empirical probability mass function
for a discrete positive distribution starting at 1
}
\examples{

# generate random sample of independent draws from Poisson distribution
x <- rpois(200,4)
# calculate empirical probability mass function and true probability mass function
y.emp  <- epmf(x)
y.theo <- dpois(1:max(x), 4)
plot(y.emp)
points(y.theo, col = "red")

}
