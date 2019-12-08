###################################################################################
#' Calculate negative log-likelihood
#' 
#' Used to determine theta/lambda/p values for the Kriging model in \code{\link{modelKriging}}
#' with Maximum Likelihood Estimation (MLE).
#'
#' @param xt vector, containing parameters like theta, p and lambda.
#' @param D matrix (or list of multiple matrices) of distances between training samples. In case of multiple distance matrices, theta (part of xt) has to be a vector, giving a weighting parameter for each matrix.
#' @param y vector of observations at sample locations.
#' @param useLambda whether to use nugget effect, i.e., lambda (FALSE at default).
#' @param corr whether to use nugget effect, i.e., lambda (fcorrGauss at default).
#' @param inverter string specifying method for inversion of correlation matrix ("chol", cholesky decomposition at default, any other string leads to using the solve function).
#' @param indefiniteMethod The specific method used for correction: spectrum \code{"clip"}, spectrum \code{"flip"}, spectrum \code{"square"}, spectrum \code{"diffusion"}, feature embedding "feature", nearest definite matrix "near". Default is no correction: \code{"none"}. See Zaefferer and Bartz-Beielstein (2016).
#' @param indefiniteType The general type of correction for indefiniteness: \code{"NSD"},\code{"CNSD"} or the default \code{"PSD"}. See Zaefferer and Bartz-Beielstein (2016).
#' @param indefiniteRepair boolean, whether conditions of the distance matrix (in case of \code{"NSD"},\code{"CNSD"} correction type) or correlation matrix (in case of \code{"PSD"} correction type) are repaired.
#' @param returnLikelihoodOnly boolean, whether the function should return only the likelihood, or a else a list (see return information below).
#' @param inverter string, defining the inverter to use. default \code{"chol"} is inversion via \code{chol2inv}. A different string will lead to use of \code{solve}.
#' @param ntheta number of kernel parameters.
#'
#' @return the numeric Likelihood value (if \code{returnLikelihoodOnly} is TRUE) or a list with elements:
#' \describe{
#' \item{\code{NegLnLike}}{ concentrated log-likelihood *-1 for minimising }
#' \item{\code{Psi}}{ correlation matrix}
#' \item{\code{Psinv}}{ inverse of correlation matrix (to save computation time in forrRegPredictor)}
#' \item{\code{mu}}{ MLE of model parameter mu }
#' \item{\code{yMu}}{ vector of observations y minus mu}
#' \item{\code{SSQ}}{ MLE of model parameter sigma^2}
#' \item{\code{a}}{ transformation vector for eigenspectrum transformation, see Zaefferer and Bartz-Beielstein (2016)}
#' \item{\code{U}}{ Matrix of eigenvectors for eigenspectrum transformation, see Zaefferer and Bartz-Beielstein (2016)}
#' \item{\code{isIndefinite}}{ whether the uncorrected correlation (kernel) matrix is indefinite}
#' }
#' 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#'
#' @seealso \code{\link{modelKriging}}
#' @keywords internal
###################################################################################
modelKrigingLikelihood <- function(xt,D,y,useLambda=FALSE,corr=fcorrGauss,
								indefiniteMethod="none",indefiniteType="PSD",indefiniteRepair=FALSE,
								returnLikelihoodOnly=TRUE,inverter = "chol",ntheta=1){	
	n <- dim(y)[1] #number of observations	
	isIndefinite <- NA
	
	if(is.list(D)){ #in case of multiple distance matrices 
		distanceWeights <- 10^xt[ntheta+(1:length(D))]
		D <- Reduce("+",mapply("*",D,distanceWeights,SIMPLIFY=FALSE)) #combine matrice by corresponding weight value, and compute sum of the matrices
		origD <- D # Fix Definiteness (NSDness, CNSDness) of the provided distance matrix	
		ret <- correctionDistanceMatrix(D,indefiniteType,indefiniteMethod,indefiniteRepair)
		D <- ret$mat
		isCNSD <- ret$isCNSD
		A <- ret$A	
	}
	
	if(ntheta<1){ #corr function has no parameters
		Psi <- corr(D)
	}else{
		theta <- xt[1:ntheta]
		Psi <- corr(D,theta)
	}
	if(any(is.infinite(Psi))){ # this is required especially if distance matrices are forced to be CNSD/NSD and hence have zero distances
	  penalty <- 1e4
		if(returnLikelihoodOnly){
			return(penalty)
		}
		return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite))
	}
	
	U <- a <- NA
	origPsi <- NA
	unrepairedPsi <- NA
	if(indefiniteType == "PSD"){
		origPsi <- Psi
		ret <- correctionKernelMatrix(Psi,indefiniteMethod,indefiniteRepair)
		a <- ret$a
		U <- ret$U
		isIndefinite <- !ret$isPSD
		Psi <- ret$mat				
		unrepairedPsi <- ret$matNoRep
		#check whether indef-correction somehow yielded malformed values
		if(any(is.na(Psi))){
			#warning("NaN or NA values due to failed indefiniteness-correction in (in modelKrigingLikelihood). Returning penalty.")
			penalty <- 1e4 
			if(returnLikelihoodOnly){
				return(penalty)
			}
			return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite))
		}
	}

	if(useLambda){
		lambda <- 10^xt[length(xt)];
		Psi <- Psi + diag(lambda,n) 
	}
		
	if(inverter=="chol"){
		## cholesky decomposition
		cholPsi <- try(chol(Psi), TRUE) 

		## give penalty if fail
		if(class(cholPsi)[1] == "try-error"){
			#warning("Correlation matrix is not positive semi-definite (in modelKrigingLikelihood). Returning penalty.")
			penalty <- 1e4 - min(eigen(Psi,symmetric=TRUE,only.values=TRUE)$values) #the minimal eigenvalue should push the search towards positive eigenvalues
			if(returnLikelihoodOnly){
				return(penalty)
			}
			return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite)) 
		}	
			
		#calculate natural log of the determinant of Psi (numerically more reliable and also faster than using det or determinant)
		LnDetPsi <- 2*sum(log(abs(diag(cholPsi))))
		
		#inverse with cholesky decomposed Psi
		Psinv <- try(chol2inv(cholPsi),TRUE)
		if(class(Psinv)[1] == "try-error"){
			#warning("Correlation matrix is not positive semi-definite (in modelKrigingLikelihood). Returning penalty.")
			penalty <- 1e4 - min(eigen(Psi,symmetric=TRUE,only.values=TRUE)$values) #the minimal eigenvalue should push the search towards positive eigenvalues
			if(returnLikelihoodOnly){
				return(penalty)
			}
			return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite)) 
		}
	}else{		
		Psinv <- try(solve(Psi),TRUE) #inverse with LU decomposition
		if(class(Psinv)[1] == "try-error"){
			#warning("Correlation matrix is not positive semi-definite (in modelKrigingLikelihood). Returning penalty.")
			penalty <- 1e4 - min(eigen(Psi,symmetric=TRUE,only.values=TRUE)$values) #the minimal eigenvalue should push the search towards larger eigenvalues
			if(returnLikelihoodOnly){
				return(penalty)
			}
			return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite)) 
		}		
		#calculate natural log of the determinant of Psi 
		LnDetPsi <- determinant(Psi,logarithm=TRUE)$modulus
	}
	
	## Check whether Psi is ill-conditioned
	kap <-  1 / ( max(colSums(abs(Psi))) *  max(colSums(abs(Psinv)))) # ==  rcond(Psi) 
	if(is.na(kap))
		kap <- 0
	if(kap < 1e-10){ 
		#warning("Correlation matrix is ill-conditioned (in modelKrigingLikelihood). Returning penalty.")
    penalty <- 1e4 
		if(returnLikelihoodOnly){
			return(penalty)
		}	
		return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite))
	}
	
	psisum <- sum(Psinv) #this sum of all matrix elements may sometimes become zero, which may be caused by inaccuracies. then, the following may help
	if(psisum==0){
		psisum <- as.numeric(rep(1,n) %*% Psinv %*% rep(1,n))
		if(psisum==0){ #if it is still zero, return penalty
			#warning("Sum of elements in inverse correlation matrix is zero (in modelKrigingLikelihood). Returning penalty.")
			penalty <- 1e4 
			if(returnLikelihoodOnly){
				return(penalty)
			}
			return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite))
		}
	}		
	mu <- sum(Psinv%*%y)/psisum
	if(is.infinite(mu)|is.na(mu)){ 
		#warning("MLE estimate of mu is infinite or NaN. Returning penalty.")
		penalty <- 1e4 
		if(returnLikelihoodOnly){
			return(penalty)
		}		
		return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite))
	}			
	yMu <- y-mu 
	#SigmaSqr <- (t(yMu)%*%Psinv%*%yMu)/n
	SigmaSqr <- crossprod(yMu,Psinv)%*%yMu/n
	if(SigmaSqr < 0){ 
		#warning("Maximum Likelihood Estimate of model parameter sigma^2 is negative (in modelKrigingLikelihood). Returning penalty. ")
		penalty <- as.numeric(1e4-SigmaSqr)
		if(returnLikelihoodOnly){
			return(penalty)
		}
		return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite)) 
	}
	NegLnLike <- n*log(SigmaSqr) + LnDetPsi
	if(is.na(NegLnLike)|is.infinite(NegLnLike)){#this may happen eg if all y are 0
		penalty <- 1e4 
		if(returnLikelihoodOnly){
			return(penalty)
		}		
		return(list(NegLnLike=penalty,origPsi=NA,Psi=NA,unrepairedPsi=NA,Psinv=NA,mu=NA,yMu=NA,SSQ=NA,a=NA,U=NA,isIndefinite=isIndefinite)) 
	}
	if(returnLikelihoodOnly){
		return(as.numeric(NegLnLike))
	}		
	ret <- list(NegLnLike=NegLnLike,origPsi=origPsi,Psi=Psi,unrepairedPsi=unrepairedPsi,Psinv=Psinv,mu=mu,yMu=yMu,SSQ=SigmaSqr,a=a,U=U,isIndefinite=isIndefinite)
	if(exists("origD")){
		ret$D <- D
		ret$origD <- origD
		ret$isCNSD <- isCNSD
		ret$A <- A
	}
	ret	
}


###################################################################################
#' Calculate negative log-likelihood
#' 
#' This is a wrapper for the Kriging likelihood function \code{\link{modelKrigingLikelihood}}.
#' It is intended for the case where parameters of the distance function are also optimized
#' during maximum likelihood estimation. Thus, the wrapper receives the data, computes the
#' parameterized distance matrix and passes it to the standard likelihood function.
#'
#' @param xt vector, containing parameters like theta, p and lambda.
#' @param xs training samples, which are the input for the distance function. Should be in list format.
#' @param ys vector of observations at training sample locations.
#' @param useLambda whether to use nugget effect, i.e., lambda (FALSE at default).
#' @param corr whether to use nugget effect, i.e., lambda (fcorrGauss at default).
#' @param inverter string specifying method for inversion of correlation matrix ("chol", cholesky decomposition at default, any other string leads to using the solve function).
#' @param indefiniteMethod The specific method used for correction: spectrum \code{"clip"}, spectrum \code{"flip"}, spectrum \code{"square"}, spectrum \code{"diffusion"}, feature embedding "feature", nearest definite matrix "near". Default is no correction: \code{"none"}. See Zaefferer and Bartz-Beielstein (2016).
#' @param indefiniteType The general type of correction for indefiniteness: \code{"NSD"},\code{"CNSD"} or the default \code{"PSD"}. See Zaefferer and Bartz-Beielstein (2016).
#' @param indefiniteRepair boolean, whether conditions of the distance matrix (in case of \code{"NSD"},\code{"CNSD"} correction type) or correlation matrix (in case of \code{"PSD"} correction type) are repaired.
#' @param returnLikelihoodOnly boolean, whether the function should return only the likelihood, or a else a list (see return information below).
#' @param distanceFunction the distance function. 
#' @param combineDistances boolean, whether to combine several distances provided as a list of distance functions.
#' @param distanceParametersLower lower boundary for the distance function(s) parameters. A vector in case of one distance, a list of vectors in case of several functions. The parameters are passed as a vector to each respective distance function.
#' @param ntheta number of kernel parameters.
#' @param scaling boolean, whether to scale the distance matrix.
#'
#' @return the numeric Likelihood value (if \code{returnLikelihoodOnly} is TRUE) or a list with elements:
#' \describe{
#' \item{\code{NegLnLike}}{ concentrated log-likelihood *-1 for minimising}
#' \item{\code{Psi}}{ correlation matrix}
#' \item{\code{Psinv}}{ inverse of correlation matrix (to save computation time in forrRegPredictor)}
#' \item{\code{mu}}{ MLE of model parameter mu }
#' \item{\code{yMu}}{ vector of observations y minus mu}
#' \item{\code{SSQ}}{ MLE of model parameter sigma^2}
#' \item{\code{a}}{ transformation vector for eigenspectrum transformation, see Zaefferer and Bartz-Beielstein (2016)}
#' \item{\code{U}}{ Matrix of eigenvectors for eigenspectrum transformation, see Zaefferer and Bartz-Beielstein (2016)}
#' \item{\code{isIndefinite}}{ whether the uncorrected correlation (kernel) matrix is indefinite}
#' } 
#'
#' @seealso \code{\link{modelKrigingLikelihood}}
#' @keywords internal
###################################################################################
modelKrigingParameterizedLikelihood <- function(xt,xs,ys,useLambda=FALSE,corr=fcorrGauss,
								indefiniteMethod="none",indefiniteType="PSD",indefiniteRepair=FALSE,returnLikelihoodOnly=TRUE,inverter = "chol",
								distanceFunction,combineDistances,distanceParametersLower,ntheta,scaling){
	#######
	nd <- length(distanceFunction) # number of distance functions
	
	if(combineDistances & nd>1){
		nweights=nd #number of weight parameters for combining distances
	}else{
		nweights=0 #number of weight parameters for combining distances
	}
	
	# parameters of the distance function(s)
	if(ntheta > 0 | nweights > 0 | useLambda)
		distanceParameters <- xt[-(1:(ntheta+nweights+useLambda))]
	else
		distanceParameters <- xt
		
	#calculate distance matrix
	ret <- modelKrigingDistanceCalculation(xs,distanceFunction=distanceFunction,parameters=distanceParameters,
							distances=NULL,scaling=scaling,combineDistances=combineDistances,indefiniteMethod=indefiniteMethod,
							indefiniteType=indefiniteType,indefiniteRepair=indefiniteRepair,distanceParametersLower)

	#Call ordinary likelihood function
	ret2	<- modelKrigingLikelihood(xt=xt[1:(ntheta+nweights+useLambda)],ret$D,ys,useLambda,corr,
				indefiniteMethod,indefiniteType,indefiniteRepair,returnLikelihoodOnly,inverter,ntheta)

	if(returnLikelihoodOnly){
		return(ret2)
	}
	
	ret2[names(ret)] <- ret
	ret2
}
#todo: testing!

