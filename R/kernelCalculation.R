###################################################################################
#' Calculate Kernel Matrix
#'
#' Calculate the similarities between all samples in a list, and return as matrix.
#'
#' @param X list of samples, where each list element is a suitable input for \code{kernFun}
#' @param kernFun Kernel function of type f(x,y)=r, where r is a scalar and x and y are elements whose similarity is evaluated.
#' @param ... further arguments passed to distFun
#'
#' @return The similarity / kernel matrix
#'
#' @examples
#' x <- list(5:1,c(2,4,5,1,3),c(5,4,3,1,2), sample(5))
#' kernFun <- function(x,y){
#' 		exp(-distancePermutationHamming(x,y))
#' }
#' kernelMatrix(x,distancePermutationHamming)
#'
#' @export
###################################################################################
kernelMatrix <-function(X,kernFun,...){
	n <- length(X)
	m <- matrix(0,nrow=n, ncol=n)
	for(i in seq_len(n - 1))
		m[seq(i+1, n),i] <- m[i,seq(i+1, n)] <- distanceVector(X[[i]],X[seq(i+1, n)],kernFun,...) #todo: distancevector?
	m
}#todo unfinished work