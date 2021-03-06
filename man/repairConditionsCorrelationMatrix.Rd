% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indefiniteLearning.R
\name{repairConditionsCorrelationMatrix}
\alias{repairConditionsCorrelationMatrix}
\title{Repair Conditions of a Correlation Matrix}
\usage{
repairConditionsCorrelationMatrix(mat)
}
\arguments{
\item{mat}{symmetric, PSD distance matrix. If your matrix is not CNSD, use \code{\link{correctionDefinite}} first. Or use \code{\link{correctionKernelMatrix}}.}
}
\value{
repaired correlation matrix
}
\description{
This function repairs correlation matrices, so that the following two properties are ensured:
The correlations values should be between -1 and 1, and the diagonal values should be one.
}
\examples{
x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
D <- distanceMatrix(x,distancePermutationInsert)
K <- exp(-0.01*D)
K <- correctionDefinite(K,type="PSD")$mat
K
K <- repairConditionsCorrelationMatrix(K)
}
\references{
Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer.
}
\seealso{
\code{\link{correctionDefinite}}, \code{\link{correctionDistanceMatrix}}, \code{\link{correctionKernelMatrix}}, \code{\link{correctionCNSD}}, \code{\link{repairConditionsDistanceMatrix}}
}
