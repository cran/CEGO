###################################################################################
#' Nearest CNSD matrix
#'
#' This function
#' implements the alternating projection algorithm by Glunt et al. (1990) to calculate the nearest conditionally
#' negative semi-definite (CNSD) matrix (or: the nearest Euclidean distance matrix).
#' The function is similar to the \code{\link[Matrix]{nearPD}} function from the \code{Matrix} package, 
#' which implements a very similar algorithm for finding the nearest Positive Semi-Definite (PSD) matrix.
#'
#' @param x symmetric matrix, to be turned into a CNSD matrix.
#' @param eig.tol eigenvalue torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#' @param conv.tol convergence torelance value. The algorithm stops if the norm of the difference between two iterations is below this value.
#' @param maxit maximum number of iterations. The algorithm stops if this value is exceeded, even if not converged.
#' @param conv.norm.type type of norm, by default the F-norm (Frobenius). See \code{\link[base]{norm}} for other choices.
#'
#' @return list with:
#' \describe{
#' \item{\code{mat}}{ nearestCNSD matrix}
#' \item{\code{normF}}{ F-norm between original and resulting matrices}
#' \item{\code{iterations}}{ the number of performed}
#' \item{\code{rel.tol}}{ the relative value used for the tolerance convergence criterion}
#' \item{\code{converged}}{ a boolean that records whether the algorithm}
#' }
#'
#' @seealso \code{\link[Matrix]{nearPD}}, \code{\link{correctionCNSD}}, \code{\link{correctionDistanceMatrix}}
#' 
#' @export
#' @examples
#' # example using Insert distance with permutations:
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' print(D)
#' is.CNSD(D)
#' nearD <- nearCNSD(D)
#' print(nearD)
#' is.CNSD(nearD$mat)
#' # or example matrix from Glunt et al. (1990):
#' D <- matrix(c(0,1,1,1,0,9,1,9,0),3,3)
#' print(D)
#' is.CNSD(D)
#' nearD <- nearCNSD(D)
#' print(nearD)
#' is.CNSD(nearD$mat)
#' # note, that the resulting values given by Glunt et al. (1990) are 19/9 and 76/9
#' @references Glunt, W.; Hayden, T. L.; Hong, S. and Wells, J. An alternating projection algorithm for computing the nearest Euclidean distance matrix, SIAM Journal on Matrix Analysis and Applications, SIAM, 1990, 11, 589-600
###################################################################################
nearCNSD <- function (x, eig.tol = 1e-8, conv.tol = 1e-8, maxit = 1000, conv.norm.type = "F") {
	n <- ncol(x)
	X <- -x
	iter <- 0
	converged <- FALSE
	conv <- Inf
	
  #construct transformation matrix (for Projection 1)
	v <- cbind(c(rep(1,n-1),1+sqrt(n)))
	Q <- diag(n) - 2 * (1/as.numeric(crossprod(v,v))) * tcrossprod(v,v)

  #main loop
	while (iter < maxit && !converged) {
    ## store result of last iteration
    Y <- X
    ##
    ### Compute first Projection P1: Projection to CPSD matrix
    ##
    #Compute F
    Fq <- Q %*% X %*% Q #transform
    #Eigen decomposition of submatrix F1
    F1 <- Fq[-n,-n]
    e <- eigen(F1,symmetric=TRUE)
    U <- e$vectors
    lambda <- e$values	
    p <- lambda > eig.tol * lambda[1] 
    if (!any(p)) 
      stop("Matrix seems conditionally positive semi-definite")
    U <- U[, p, drop = FALSE]		
    F1q <- tcrossprod(U * rep(lambda[p], each = nrow(U)), U) #U %*% diag(pmax(lambda,0)) %*% t(U)
    Fq[-n,-n] <- F1q #replace NON-PSD F1 in F by PSD F1q
    P1F <- Q %*% Fq %*% Q #compute D from F, yielding CPSD matrix
    
    ##
    ###Compute second Projection P2: P_2(F) = F - diag(F)
    ##
    P2P1F <- P1F
    diag(P2P1F) <- 0
            
    ## correction
    X <- Y + (P2P1F - P1F)
    
    ## calculate convergence rate based on chosen norm
    conv <- norm(Y - X, conv.norm.type)
    converged <- (conv <= conv.tol) # determine whether converged
    
    # counter of iterations
    iter <- iter + 1 
  }
	if (!converged) 
		warning(gettextf("'nearCNSD()' did not converge in %d iterations",iter))
	list(mat = -P2P1F, normF = norm(x - -P2P1F, "F"), iterations = iter, 
		rel.tol = conv, converged = converged)
}


