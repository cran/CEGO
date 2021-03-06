% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutationDistances.r
\name{distancePermutationR}
\alias{distancePermutationR}
\title{R-Distance for Permutations}
\usage{
distancePermutationR(x, y)
}
\arguments{
\item{x}{first permutation (integer vector)}

\item{y}{second permutation (integer vector)}
}
\value{
numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
}
\description{
R distance or unidirectional adjacency distance. Based on count of number of times that a two element sequence in x also occurs in y, in the same order.
}
\examples{
x <- 1:5
y <- c(1,2,3,5,4)
distancePermutationR(x,y)
p <- replicate(10,sample(1:5),simplify=FALSE)
distanceMatrix(p,distancePermutationR)

}
\references{
Sevaux, Marc, and Kenneth Soerensen. "Permutation distance measures for memetic algorithms with population management." Proceedings of 6th Metaheuristics International Conference (MIC'05). 2005.

Reeves, Colin R. "Landscapes, operators and heuristic search." Annals of Operations Research 86 (1999): 473-490.
}
