#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Calculate Distance Matrix
#'
#' Calculate the distance between all samples in a list, and return as matrix.
#'
#' @param X list of samples, where each list element is a suitable input for \code{distFun}
#' @param distFun Distance function of type f(x,y)=r, where r is a scalar and x and y are elements whose distance is evaluated.
#'
#' @return The distance matrix
#'
#' @examples
#' x <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2), sample(5))
#' distanceMatrix(x,distancePermutationHamming)
#'
#' @export
#' @keywords internal
###################################################################################
distanceMatrix <-function(X,distFun){
	n <- length(X)
	m <- matrix(0,nrow=n, ncol=n)
	for(i in seq_len(n - 1))
		m[seq(i+1, n),i] <- m[i,seq(i+1, n)] <- distanceVector(X[[i]],X[seq(i+1, n)],distFun)
	m
}

###################################################################################
#' Calculate Distance Vector
#'
#' Calculate the distance between a single sample and all samples in a list.
#'
#' @param a A single sample which is a suitable input for \code{distFun}
#' @param X list of samples, where each list element is a suitable input for \code{distFun}
#' @param distFun Distance function of type f(x,y)=r, where r is a scalar and x and y are elements whose distance is evaluated.
#'
#' @return A numerical vector of distances
#'
#' @examples
#' x <- 1:5
#' y <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2))
#' distanceVector(x,y,distancePermutationHamming)
#'
#' @export
#' @keywords internal
###################################################################################
distanceVector <-function(a,X,distFun){
	unlist(lapply(X,distFun,a))
}
