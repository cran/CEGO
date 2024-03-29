% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distanceCalculation.R
\name{distanceRealEuclidean}
\alias{distanceRealEuclidean}
\title{Euclidean Distance}
\usage{
distanceRealEuclidean(x, y)
}
\arguments{
\item{x}{first real vector}

\item{y}{second real vector}
}
\value{
numeric distance value \deqn{d(x,y)}
}
\description{
The Euclidean distance for real vectors.
}
\examples{
x <- runif(5)
y <- runif(5)
distanceRealEuclidean(x,y)

}
