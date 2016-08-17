###################################################################################
#' Check for Conditional Negative Semi-Definiteness
#' 
#' This function checks whether a symmetric matrix is Conditionally Negative Semi-Definite (CNSD).
#' Note that this function does not check whether the matrix is actually symmetric.
#'
#' @param X a symmetric matrix
#' @param method a string, specifiying the method to be used. \code{"alg1"} is based on algorithm 1 in Ikramov and Savel'eva (2000). 
#' \code{"alg2"} is based on theorem 3.2 in Ikramov and Savel'eva (2000). \code{"eucl"} is based on Glunt (1990). 
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' Symmetric, CNSD matrices are, e.g., euclidean distance matrices, whiche are required to produce Positive Semi-Definite correlation
#' or kernel matrices. Such matrices are used in models like Kriging or Support Vector Machines.
#'
#' @return boolean, which is TRUE if X is CNSD
#'
#' @seealso \code{\link{is.NSD}}, \code{\link{is.PSD}}
#' @references Ikramov, K. and Savel'eva, N. Conditionally definite matrices, Journal of Mathematical Sciences, Kluwer Academic Publishers-Plenum Publishers, 2000, 98, 1-50
#' @references Glunt, W.; Hayden, T. L.; Hong, S. and Wells, J. An alternating projection algorithm for computing the nearest Euclidean distance matrix, SIAM Journal on Matrix Analysis and Applications, SIAM, 1990, 11, 589-600
#' @examples
#' # The following permutations will produce
#' # a non-CNSD distance matrix with Insert distance
#' # and a CNSD distance matrix with Hamming distance
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' is.CNSD(D,"alg1")
#' is.CNSD(D,"alg2")
#' is.CNSD(D,"eucl")
#' D <- distanceMatrix(x,distancePermutationHamming)
#' is.CNSD(D,"alg1")
#' is.CNSD(D,"alg2")
#' is.CNSD(D,"eucl")
#' @export
###################################################################################
is.CNSD <- function(X,method="alg1",tol=1e-8){
	n <- nrow(X)
	if(method=="alg1"){ #algorithm 1 (ikramov2000)
		P <- (diag(n)-matrix(1,n,n)/n)
		P[n,] <- c(rep(0,n-1),1)
		Xhat <- P %*% X %*% t(P)
		eigs <- eigen(Xhat[-n,-n],T,T)$values
		cnsd <- !eigs[1]>tol
	}else if(method=="eucl"){ #Glunt1990
		P <- (diag(n)-matrix(1,n,n)/n)
		Xhat <- P %*% X %*% P 
		eigs <- eigen(Xhat,T,T)$values
		cnsd <- !eigs[1]>tol
	}else if(method=="alg2"){ #algorithm 2, theorem 3.2 (ikramov2000)
		eigs <- eigen(rbind(cbind(X,rep(1,n)),c(rep(1,n),0)),T,T)$values
		cnsd <- !(eigs[1]*eigs[2])>tol
	}	
	cnsd
}

###################################################################################
#' Check for Positive Semi-Definiteness
#' 
#' This function checks whether a symmetric matrix is Positive Semi-Definite (PSD).
#' That means, it is determined whether all eigenvalues of the matrix are non-negative.
#' Note that this function does not check whether the matrix is actually symmetric.
#'
#' @param X a symmetric matrix
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' Symmetric, PSD matrices are, e.g., correlation
#' or kernel matrices. Such matrices are used in models like Kriging or Support Vector regression.
#'
#' @return boolean, which is TRUE if X is PSD
#'
#' @seealso \code{\link{is.CNSD}}, \code{\link{is.NSD}}
#' @examples
#' # The following permutations will produce
#' # a non-PSD kernel matrix with Insert distance
#' # and a PSD distance matrix with Hamming distance
#' # (for the given theta value of 0.01)
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' K <- exp(-0.01*distanceMatrix(x,distancePermutationInsert))
#' is.PSD(K)
#' K <- exp(-0.01*distanceMatrix(x,distancePermutationHamming))
#' is.PSD(K)
#' @export
###################################################################################
is.PSD <- function(X,tol=1e-8){
	eigen(X,T,T)$values[nrow(X)] >= -tol
}

###################################################################################
#' Check for Negative Semi-Definiteness
#' 
#' This function checks whether a symmetric matrix is Negative Semi-Definite (NSD).
#' That means, it is determined whether all eigenvalues of the matrix are non-positive.
#' Note that this function does not check whether the matrix is actually symmetric.
#'
#' @param X a symmetric matrix
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' Symmetric, NSD matrices are, e.g., correlation
#' or kernel matrices. Such matrices are used in models like Kriging or Support Vector regression.
#'
#' @return boolean, which is TRUE if X is NSD
#'
#' @seealso \code{\link{is.CNSD}}, \code{\link{is.PSD}}
#' @examples
#' # The following permutations will produce
#' # a non-PSD kernel matrix with Insert distance
#' # and a PSD distance matrix with Hamming distance
#' # (for the given theta value of 0.01)-
#' # The respective negative should be (non-) NSD
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' K <- exp(-0.01*distanceMatrix(x,distancePermutationInsert))
#' is.NSD(-K) 
#' K <- exp(-0.01*distanceMatrix(x,distancePermutationHamming))
#' is.NSD(-K)
#' @export
###################################################################################
is.NSD <- function(X,tol=1e-8){
	eigen(X,T,T)$values[1] <= tol
}
