% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimEA.R
\name{recombinationSelfAdapt}
\alias{recombinationSelfAdapt}
\title{Self-adaptive recombination operator}
\usage{
recombinationSelfAdapt(population, parameters)
}
\arguments{
\item{population}{List of permutations}

\item{parameters}{list, contains the available single mutation functions (\code{mutationFunctions}), 
and a data.frame that collects the chosen function and mutation rate for each individual (\code{selfAdapt}).}
}
\description{
This recombination function selects an operator (provided in parameters$recombinationFunctions)
based on self-adaptive parameters chosen for each individual separately.
}
\seealso{
\code{\link{optimEA}}, \code{\link{mutationSelfAdapt}}
}
