% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelKrigingPredict.R
\name{computeCorrelationMatrix}
\alias{computeCorrelationMatrix}
\title{Compute Correlation Matrix}
\usage{
computeCorrelationMatrix(object, x)
}
\arguments{
\item{object}{fit of the Kriging model (settings and parameters), of class \code{modelKriging}.}

\item{x}{list of samples / data}
}
\value{
the correlation matrix
}
\description{
Compute the correlation matrix of samples x, given the model object.
}
\seealso{
\code{\link{simulate.modelKriging}}

\code{\link{predict.modelKriging}}
}
\keyword{internal}
