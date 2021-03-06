% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutationDistances.r
\name{distancePermutationInterchange}
\alias{distancePermutationInterchange}
\title{Interchange Distance for Permutations}
\usage{
distancePermutationInterchange(x, y)
}
\arguments{
\item{x}{first permutation (integer vector)}

\item{y}{second permutation (integer vector)}
}
\value{
numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
}
\description{
The interchange distance is an edit-distance, counting how many edit operation (here: interchanges, i.e., transposition of two arbitrary elements) have to be
performed to transform permutation x into permutation y.
}
\examples{
x <- 1:5
y <- c(1,4,3,2,5)
distancePermutationInterchange(x,y)
p <- replicate(10,sample(1:5),simplify=FALSE)
distanceMatrix(p,distancePermutationInterchange)

}
\references{
Schiavinotto, Tommaso, and Thomas Stuetzle. "A review of metrics on permutations for search landscape analysis." Computers & operations research 34.10 (2007): 3143-3153.
}
