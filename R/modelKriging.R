###################################################################################
#' Kriging Model
#' 
#' Implementation of a distance-based Kriging model, e.g., for mixed or combinatorial input spaces.
#' It is based on employing suitable distance measures for the samples in input space.
#'
#' The basic Kriging implementation is based on the work of Forrester et al. (2008). 
#' For adaptation of Kriging to mixed or combinatorial spaces, as well as
#' choosing distance measures with Maximum Likelihood Estimation, see the other two references (Zaefferer et al., 2014).
#'
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric.  It can also be a list of several distance functions. In this case, Maximum Likelihood Estimation (MLE) is used 
#'		to determine the most suited distance measure.
#' @param control (list), with the options for the model building procedure:\cr
#' \code{lower} lower boundary for theta, default is \code{1e-6}\cr
#' \code{upper} upper boundary for theta, default is \code{100}\cr
#' \code{corr} function to be used for correlation modelling, default is \code{fcorrGauss}\cr
#' \code{algTheta}  algorithm used to find theta (as well as p and lambda), default is \code{\link{optimInterface}}.\cr
#' \code{algThetaControl}  list of controls passed to \code{algTheta}.\cr
#' \code{optimizeP} boolean that specifies whether the exponents (\code{p}) should be optimized. Else they will be set to two. \cr
#' \code{useLambda} whether or not to use the regularization constant lambda (nugget effect). Default is \code{FALSE}.\cr
#' \code{lambdaLower} lower boundary for lambda (log scale), default is \code{-6}\cr 
#' \code{lambdaUpper} upper boundary for lambda (log scale), default is \code{0}\cr
#' \code{distances} a distance matrix. If available, this matrix is used for model building, instead of calculating the distance matrix using the parameters \code{distanceFunction}. Default is \code{NULL}.
#' \code{scaling} If TRUE: Distances values are divided by maximum distance to avoid scale bias.\cr
#' \code{reinterpolate} If TRUE: reinterpolation is used to generate better uncertainty estimates in the presence of noise. \cr
#' \code{combineDistances} By default, several distance functions or matrices are subject to a likelihood based decision, choosing one. If this parameter is TRUE, they are instead combined by determining a weighted sum. The weighting parameters are determined by MLE.\cr
#' \code{userParameters} By default: (\code{NULL}). Else, this vector is used instead of MLE to specify the model parameters, in the following form: \code{log10(c(theta,lambda,p))}. In case of multiple combined distance functions, \code{theta} is also a vector, one element for each distance function.\cr
#' \code{indefiniteMethod} The specific method used for correction: spectrum \code{"clip"}, spectrum \code{"flip"}, spectrum \code{"square"}, spectrum \code{"diffusion"}, feature embedding "feature", nearest definite matrix "near". Default is no correction: \code{"none"}. See Zaefferer and Bartz-Beielstein (2016).\cr
#' \code{indefiniteType}  The general type of correction for indefiniteness: \code{"NSD"},\code{"CNSD"} or the default \code{"PSD"}. See Zaefferer and Bartz-Beielstein (2016). Note, that feature embedding may not work in case of multiple distance functions.\cr
#' \code{indefiniteRepair} boolean, whether conditions of the distance matrix (in case of \code{"NSD"},\code{"CNSD"} correction type) or correlation matrix (in case of \code{"PSD"} correction type) are repaired.
#'
#' @return an object of class \code{modelKriging} containing the options (see control parameter) and determined parameters for the model:\cr
#' \code{theta} activity or width parameter theta, a parameter of the correlation function determined with MLE\cr
#' \code{log10Theta} log10 \code{theta} (i.e. \code{log10(theta)})\cr
#' \code{lambda} regularization constant (nugget) lambda \cr
#' \code{log10Lambda} log10 of regularization constant (nugget) lambda (i.e. \code{log10(lambda)})\cr
#' \code{p} exponent p, parameter of the correlation function determined with MLE (if \code{optimizeP} is \code{TRUE})\cr
#' \code{yMu} vector of observations y, minus MLE of mu\cr
#' \code{SSQ} Maximum Likelihood Estimate (MLE) of model parameter sigma^2\cr
#' \code{mu} MLE of model parameter mu\cr
#' \code{Psi} correlation matrix Psi\cr
#' \code{Psinv} inverse of Psi\cr
#' \code{nevals} number of Likelihood evaluations during MLE of theta/lambda/p\cr
#' \code{distanceFunctionIndexMLE} If a list of several distance measures (\code{distanceFunction}) was given, this parameter contains the index value of the measure chosen with MLE.
#' 
#' @seealso \code{\link{predict.modelKriging}} 
#' 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#' @references Zaefferer, Martin and Bartz-Beielstein, Thomas (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#'
#' @examples
#' # Set random number generator seed
#' set.seed(1)
#' # Simple test landscape
#' fn <- landscapeGeneratorUNI(1:5,distancePermutationHamming)
#' # Generate data for training and test
#' x <- unique(replicate(40,sample(5),FALSE))
#' xtest <- x[-(1:15)]
#' x <- x[1:15]
#' # Determin true objective function values
#' y <- fn(x)
#' ytest <- fn(xtest)
#' # Build model
#' fit <- modelKriging(x,y,distancePermutationHamming,
#'     control=list(algThetaControl=list(method="L-BFGS-B"),useLambda=FALSE))
#' # Predicted obj. function values
#' ypred <- predict(fit,xtest)$y
#' # Uncertainty estimate
#' fit$predAll <- TRUE
#' spred <- predict(fit,xtest)$s
#' # Plot
#' plot(ytest,ypred,xlab="true value",ylab="predicted value",
#'     pch=20,xlim=c(0.3,1),ylim=c(min(ypred)-0.1,max(ypred)+0.1))
#' segments(ytest, ypred-spred,ytest, ypred+spred)
#' epsilon = 0.02
#' segments(ytest-epsilon,ypred-spred,ytest+epsilon,ypred-spred)
#' segments(ytest-epsilon,ypred+spred,ytest+epsilon,ypred+spred)
#' abline(0,1,lty=2)
#' # Use a different/custom optimizer (here: SANN) for maximum likelihood estimation: 
#' # (Note: Bound constraints are recommended, to avoid Inf values.
#' # This is really just a demonstration. SANN does not respect bound constraints.)
#' optimizer1 <- function(x,fun,lower=NULL,upper=NULL,control=NULL,...){
#'   res <- optim(x,fun,method="SANN",control=list(maxit=100),...)
#'   list(xbest=res$par,ybest=res$value,count=res$counts)
#' }
#' fit <- modelKriging(x,y,distancePermutationHamming,
#'                    control=list(algTheta=optimizer1,useLambda=FALSE))
#' #One-dimensional optimizer (Brent). Note, that Brent will not work when 
#' #several parameters have to be set, e.g., when using nugget effect (lambda).
#' #However, Brent may be quite efficient otherwise.
#' optimizer2 <- function(x,fun,lower,upper,control=NULL,...){
#'  res <- optim(x,fun,method="Brent",lower=lower,upper=upper,...)
#'  list(xbest=res$par,ybest=res$value,count=res$counts)
#' }
#' fit <- modelKriging(x,y,distancePermutationHamming,
#'                     control=list(algTheta=optimizer2,useLambda=FALSE))
#' @export
###################################################################################
modelKriging <- function(x, y, distanceFunction,control=list()){ 
#TODO: use of tcrossprod or crossprod for speedup??
	con<-list(lower=1e-6, upper=1e5, 
						corr=fcorrGauss, 
						algTheta= optimInterface, 
						algThetaControl= list(funEvals=200,reltol=1e-4,factr=1e12,restarts=TRUE),
						optimizeP= FALSE, 
						combineDistances=FALSE, 
						useLambda=FALSE, lambdaLower = -6, lambdaUpper = 0, 
            indefiniteMethod= "none", indefiniteType="PSD", indefiniteRepair=TRUE,
						scaling=FALSE,reinterpolate=FALSE);
	con$algThetaControl[names(control$algThetaControl)] <- control$algThetaControl
	control$algThetaControl <- con$algThetaControl
	con[names(control)] <- control;
 	control<-con;
		
	algThetaControl <- control$algThetaControl
	useLambda <- control$useLambda
	lambdaLower <- control$lambdaLower
	lambdaUpper <- control$lambdaUpper
	combineDistances <- control$combineDistances
	fcorr <- control$corr
 	fit <- control

	if(!is.matrix(y))
		y <- as.matrix(y)
	if(any(duplicated(x)) & !control$useLambda){ #duplicates HAVE to be removed, but duplicates for noisy problems are okay. 
		duplicates <- which(duplicated(x))
		x <- x[-duplicates]
		y <- as.matrix(y[-duplicates])
	}	
	
	#
	##
	#
	if(control$indefiniteMethod=="feature"){ #distances as features #TODO what to do in case of multiple distance functions ?
		if(is.null(control$distances))
			D <-distanceMatrix(x,distanceFunction) 
		else
			D <- control$distances
		fit$origx <- x
		fit$origDistanceFunction <- distanceFunction
		x <- split(D,seq(nrow(D))) #each distance vector in the distance matrix is now a feature vector
		distanceFunction <- function(x,y){  #now use positive semidefinite euclidean distance on feature vector
			sqrt(sum((x-y)^2))
		}
		control$distances <- distanceMatrix(x,distanceFunction) 
	}	
	#
	##
	#
	fit$x <- x
	fit$y <- y

	n <- length(fit$x) #number of observations
	nd <- length(distanceFunction) # number of distance functions
	
	#calculate distance matrix
	if(nd==1){ #one distance function
		if(is.null(control$distances))
			D <-distanceMatrix(x,distanceFunction) 
		else
			D <- control$distances
    maxD <- max(D) #maximum distance
    if(control$scaling){
      D <- D/maxD
    }    
	}else{ #multiple distance functions
		if(is.null(control$distances)){
			D <- list()
      maxD <- list()
			for(i in 1:nd){
        D[[i]] <-distanceMatrix(x,distanceFunction[[i]]) 
        maxD[[i]] <- max(D[[i]]) #maximum distance
        if(control$scaling){
          D[[i]] <- D[[i]]/maxD[[i]]
        } 
			}
		}else{
			D <- control$distances
      maxD <- list()
			for(i in 1:nd){
        maxD[[i]] <- max(D[[i]]) #maximum distance
        if(control$scaling){
          D[[i]] <- D[[i]]/maxD[[i]]
        } 
			}    
		}
	}
  
	# Fix Definiteness (NSDness, CNSDness) of the provided distance matrix/matrices	
	fit$origD <- D 
	if(nd==1){#in case of one distance function
		ret <- correctionDistanceMatrix(D,control$indefiniteType,control$indefiniteMethod,control$indefiniteRepair)
		D <- ret$mat
		fit$D <- ret$mat
		fit$isCNSD <- ret$isCNSD
		fit$A <- ret$A	
	}else if(!combineDistances){ #in case of multiple distances, which are not combined (but chosen from):
		fit$D <- list()
		fit$isCNSD <- list()
		fit$A <- list()	
		for(i in 1:nd){
			ret <- correctionDistanceMatrix(D[[i]],control$indefiniteType,control$indefiniteMethod,control$indefiniteRepair)
			D[[i]] <- ret$mat
			fit$D[[i]] <- ret$mat
			fit$isCNSD[[i]] <- ret$isCNSD
			fit$A[[i]] <- ret$A	
		}
	}
	
	if(is.null(control$userParameters)){ 
		# start point for theta, and bounds:	
		res <- modelKrigingInit(fit$startTheta,log10(fit$lower),log10(fit$upper),fit$optimizeP,useLambda,lambdaLower,lambdaUpper,combineDistances,nd)
		x0 <- res$x0
		lowerTheta <- res$lower
		upperTheta <- res$upper
	  
		# adapt tuning (MLE) budget to dimensionality of parameter space
		algThetaControl$funEvals <- algThetaControl$funEvals*length(x0)	
		if(combineDistances | nd==1){
			res <- control$algTheta(x=x0,fun=modelKrigingLikelihood,lower=lowerTheta,upper=upperTheta,
							control=algThetaControl,D=D,y=fit$y,optimizeP=fit$optimizeP,useLambda=useLambda,corr=fcorr,
							indefiniteMethod=control$indefiniteMethod,indefiniteType=control$indefiniteType,indefiniteRepair=control$indefiniteRepair,returnLikelihoodOnly=TRUE)	
			fit$distanceFunction <- distanceFunction
		}else{
			res <- list()
			minlik=Inf
			minlikindex=1
			for(i in 1:length(distanceFunction)){
				res[[i]] <- control$algTheta(x=x0,fun=modelKrigingLikelihood,lower=lowerTheta,upper=upperTheta,
							control=algThetaControl,D=D[[i]],y=fit$y,optimizeP=fit$optimizeP,useLambda=useLambda,corr=fcorr,
							indefiniteMethod=control$indefiniteMethod,indefiniteType=control$indefiniteType,indefiniteRepair=control$indefiniteRepair,returnLikelihoodOnly=TRUE)	
				if(res[[i]]$ybest < minlik){
					minlik <- res[[i]]$ybest
					minlikindex <- i
				}
			}
			res <- res[[minlikindex]]
			fit$distanceFunction <- distanceFunction[[minlikindex]]
			D<-D[[minlikindex]]
			maxD <- maxD[[minlikindex]]
			fit$D <- fit$D[[minlikindex]]
			fit$origD <- fit$origD[[minlikindex]]
			fit$isCNSD <- fit$isCNSD[[minlikindex]]
			fit$A <- fit$A[[minlikindex]]
			fit$distanceFunctionIndexMLE <- minlikindex
		}	
		if(is.null(res$xbest)){
			res$xbest <- x0
		}
		Params <- res$xbest
		nevals <- as.numeric(res$count[[1]])
	}else{
		Params <- control$userParameters
		nevals <- 0
		fit$distanceFunction <- distanceFunction
	}

	# extract model parameters:
	# theta
	if(combineDistances){
		fit$theta <- 10^Params[1:nd]
		fit$log10Theta <- Params[1:nd]
	}else{
		fit$theta <- 10^Params[1]
		fit$log10Theta <- Params[1]
	}
	# p
	if(fit$optimizeP){	
		fit$p <- Params[length(Params)-1]
	}	
	# lambda
	if(useLambda){
		fit$log10Lambda <- Params[length(Params)];
		fit$lambda <- 10^fit$log10Lambda
	}else{
		fit$log10Lambda <- NULL;
		fit$lambda <- 0;
	}
	
	res <- modelKrigingLikelihood(c(fit$log10Theta,fit$p, fit$log10Lambda),D,fit$y,fit$optimizeP,useLambda,fcorr,
											control$indefiniteMethod,control$indefiniteType,control$indefiniteRepair,
											returnLikelihoodOnly=FALSE)	#need to also return the correlation matrix and other elements of the model

	if(combineDistances & nd>1){
		D <- res$D
		fit$A <- res$A
		fit$D <- res$D
		fit$origD <- res$origD
		fit$isCNSD <- res$isCNSD
	}
	
	#thus, multiple distances are reduced to one: either combined, or one chosen.
	
	fit$isIndefinite <- res$isIndefinite
	fit$U <- res$U	
	fit$a <- res$a
	fit$yMu <- res$yMu
	fit$SSQ <- as.numeric(res$SSQ)
	fit$mu <- res$mu
	fit$Psi <- res$Psi 
	if(combineDistances & nd>1){
		thetatmp <- 1
	}else{
		thetatmp <- fit$theta
	}
	if(fit$optimizeP){
		fit$origPsi <- fcorr(thetatmp,D^fit$p) 
	}else{
	  fit$origPsi <- fcorr(thetatmp,D) 
	}
	fit$Psinv <- res$Psinv
	##precompute transformations
	if(control$indefiniteType=="PSD" & !fit$indefiniteRepair & fit$isIndefinite & any(control$indefiniteMethod==c("clip","flip","square","param","diffusion"))){ #RETRANSFORMATION OF THE SOLUTION ONLY  
		A <- res$U %*% diag(res$a) %*% t(res$U)
    fit$A <- A 
		#
    fit$Psinv <- t(A)%*%fit$Psinv #retransform the result	for prediction
		fit$PsinvA <- fit$Psinv %*% A #retransform the result	(for variance estimation only)
	}
	if(useLambda){ 
		PsiB <- res$Psi-diag(fit$lambda,n)+diag(.Machine$double.eps,n) 
		fit$SSQReint <- as.numeric((t(res$yMu)%*%res$Psinv%*%PsiB%*%res$Psinv%*%res$yMu)/n) #res is used intentionally, needs to be untransformed Psinv
		fit$PsinvReint <- try(chol2inv(chol(PsiB)), TRUE) 
		if(class(fit$PsinvReint) == "try-error"){
			fit$PsinvReint <- ginv(PsiB) 
		}	
		#now apply same transformations as for non-reinterpolating matrices
		if(control$indefiniteType=="PSD" & fit$isIndefinite  & !fit$indefiniteRepair & any(control$indefiniteMethod==c("clip","flip","square","param","diffusion"))){ #RETRANSFORMATION OF THE SOLUTION ONLY  
      fit$PsinvReint <- t(A)%*%fit$PsinvReint %*% A #retransform
		} 
	}
	#	
	##
	fit$nevals <- nevals
	fit$like <- res$NegLnLike
  fit$maximumDistance <- maxD
  fit$predAll <- FALSE
	class(fit)<- "modelKriging"
	if(is.na(fit$Psinv[1])){ #model building failed. no invertible correlation matrix was found. return NA fit
		stop("Building the Kriging model failed, no invertible correlation matrix was found. This may be due to the specific data-set or distance function used.")
	}else{
		return(fit)
	}	
}

###################################################################################
#' Kriging Model
#'
#' DEPRECATED version of the Kriging model, please use \code{\link{modelKriging}}
#' 
#' @param x list of samples in input space
#' @param y column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value
#' @param control options for the model building procedure
#'
#' @keywords internal
#' @export
###################################################################################
combinatorialKriging <- function(x, y, distanceFunction, control = list()){
	.Deprecated("modelKriging")
	modelKriging(x,y,distanceFunction,control)
}


###################################################################################
#' Gaussian Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrGauss <- function(theta,D){
	exp(-theta * D)
}

###################################################################################
#' Cubic Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrCubic <- function(theta,D){
	Psi <- pmin(D * theta,1)
	1 - Psi^2 * (3 - 2*Psi)
}

###################################################################################
#' Linear Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrLinear <- function(theta,D){
	pmax(1- D * theta,0)
}

###################################################################################
#' Spherical Kernel for Kriging
#'
#' @param theta kernel parameter
#' @param D distance matrix
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrSphere <- function(theta,D){
	Psi <- pmin(D * theta,1)
	 1 - Psi * (1.5 - 0.5*Psi^2)
}	

###################################################################################
#' Kriging: Initial guess and bounds
#'
#' Initialize parameter tuning for the Kriging model, setting the initial guess
#' as well as bound constraints.
#'
#' @param startTheta user provided start guess (optional).
#' @param lowerTheta lower boundary for theta values (log scale).
#' @param upperTheta upper boundary for theta values (log scale).
#' @param optimizeP boolean, whether parameter p is optimized.
#' @param useLambda boolean, whether nugget effect (lambda) is used.
#' @param lambdaLower lower boundary for lambda (log scale).
#' @param lambdaUpper upper boundary for lambda (log scale).
#' @param combineDistances boolean, whether multiple distances are combined.
#' @param nd number of distance function.
#'
#' @return a list with elements \code{x0} (start guess), \code{lower} (lower bound), \code{upper} (upper bound).
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @keywords internal
###################################################################################
modelKrigingInit	<- function(startTheta=NULL,lowerTheta,upperTheta,optimizeP,useLambda, lambdaLower, lambdaUpper, combineDistances,nd){
	if(is.null(startTheta)){
		x1 <- 0
		if(combineDistances){
			x1 <- rep(x1,nd)
			lowerTheta <- rep(lowerTheta,nd)
			upperTheta <- rep(upperTheta,nd)
		}
	}else{
		x1 <- startTheta
	}
	
	if(optimizeP){ # optimize p
		lowerTheta <- c(lowerTheta, 0.01)
		upperTheta <- c(upperTheta, 2)		
		x3 <- 1 #start values for p
		x0 <- c(x1,x3)
	}else{ # p  is fixed to 1 
		x0 <- c(x1)
	}
	if(useLambda){
		# start value for lambda:
		x2 <- lambdaLower + (lambdaUpper - lambdaLower)*runif(1)
		x0 <- c(x0,x2)
		#append regression constant lambda (nugget)
		lowerTheta <- c(lowerTheta,lambdaLower)
		upperTheta <- c(upperTheta, lambdaUpper)
	}	

	#force x0 into bounds
	x0 <- pmin(x0,upperTheta)
	x0 <- pmax(x0,lowerTheta)
	
	list(x0=x0,lower=lowerTheta,upper=upperTheta)
}

	