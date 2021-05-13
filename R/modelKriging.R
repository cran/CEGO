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
#'		The distance function may have additional parameters. For that case, see distanceParametersLower/Upper in the controls.
#'    If distanceFunction is missing, it can also be provided in the control list.
#' @param control (list), with the options for the model building procedure:
#' \describe{
#' \item{\code{lower}}{ lower boundary for theta, default is \code{1e-6}}
#' \item{\code{upper}}{ upper boundary for theta, default is \code{100}}
#' \item{\code{corr}}{ function to be used for correlation modelling, default is \code{fcorrGauss}}
#' \item{\code{algTheta}}{ algorithm used to find theta (as well as p and lambda), default is \code{\link{optimInterface}}.}
#' \item{\code{algThetaControl}}{ list of controls passed to \code{algTheta}.}
#' \item{\code{useLambda}}{ whether or not to use the regularization constant lambda (nugget effect). Default is \code{FALSE}.}
#' \item{\code{lambdaLower}}{ lower boundary for lambda (log scale), default is \code{-6}}
#' \item{\code{lambdaUpper}}{ upper boundary for lambda (log scale), default is \code{0}}
#' \item{\code{distanceParametersLower}}{ lower boundary for parameters of the distance function, default is \code{NA} which means there are no distance function parameters. If several distance functions are supplied, this should be a list of lower boundary vectors for each function.}
#' \item{\code{distanceParametersUpper}}{ upper boundary for parameters of the distance function, default is \code{NA} which means there are no distance function parameters. If several distance functions are supplied, this should be a list of upper boundary vectors for each function.}
#' \item{\code{distances}}{ a distance matrix. If available, this matrix is used for model building, instead of calculating the distance matrix using the parameters \code{distanceFunction}. Default is \code{NULL}.}
#' \item{\code{scaling}}{ If TRUE: Distances values are divided by maximum distance to avoid scale bias.}
#' \item{\code{reinterpolate}}{ If TRUE: reinterpolation is used to generate better uncertainty estimates in the presence of noise. }
#' \item{\code{combineDistances}}{ By default, several distance functions or matrices are subject to a likelihood based decision, choosing one. If this parameter is TRUE, they are instead combined by determining a weighted sum. The weighting parameters are determined by MLE.}
#' \item{\code{userParameters}}{ By default: (\code{NULL}). Else, this vector is used instead of MLE to specify the model parameters, in the following order: kernel parameters, distance weights, lambda, distance parameters.}
#' \item{\code{indefiniteMethod}}{ The specific method used for correction: spectrum \code{"clip"}, spectrum \code{"flip"}, spectrum \code{"square"}, spectrum \code{"diffusion"}, feature embedding "feature", nearest definite matrix "near". Default is no correction: \code{"none"}. See Zaefferer and Bartz-Beielstein (2016).}
#' \item{\code{indefiniteType}}{  The general type of correction for indefiniteness: \code{"NSD"},\code{"CNSD"} or the default \code{"PSD"}. See Zaefferer and Bartz-Beielstein (2016). Note, that feature embedding may not work in case of multiple distance functions.}
#' \item{\code{indefiniteRepair}}{ boolean, whether conditions of the distance matrix (in case of \code{"NSD"},\code{"CNSD"} correction type) or correlation matrix (in case of \code{"PSD"} correction type) are repaired.}
######## \item{\code{conditionalSimulation}}{ boolean, whether a later performed simulation of the fitted model should be conditional on the training data.}
#' }
#'
#' @return an object of class \code{modelKriging} containing the options (see control parameter) and determined parameters for the model:
#' \describe{
#' \item{\code{theta}}{ parameters of the kernel / correlation function determined with MLE.}
#' \item{\code{lambda}}{ regularization constant (nugget) lambda}
#' \item{\code{yMu}}{ vector of observations y, minus MLE of mu}
#' \item{\code{SSQ}}{ Maximum Likelihood Estimate (MLE) of model parameter sigma^2}
#' \item{\code{mu}}{ MLE of model parameter mu}
#' \item{\code{Psi}}{ correlation matrix Psi}
#' \item{\code{Psinv}}{ inverse of Psi}
#' \item{\code{nevals}}{ number of Likelihood evaluations during MLE of theta/lambda/p}
#' \item{\code{distanceFunctionIndexMLE}}{ If a list of several distance measures (\code{distanceFunction}) was given, this parameter contains the index value of the measure chosen with MLE.}
#' }
#' 
#' @seealso \code{\link{predict.modelKriging}} 
#' 
#' @references Forrester, Alexander I.J.; Sobester, Andras; Keane, Andy J. (2008). Engineering Design via Surrogate Modelling - A Practical Guide. John Wiley & Sons.
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282
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
	con<-list(lower=-6, upper=5, 
						corr=fcorrGauss, 
						algTheta= optimInterface, 
						algThetaControl= list(funEvals=200,reltol=1e-4,factr=1e12,restarts=TRUE),#TODO: change reltol and factr defaults?
						combineDistances=FALSE, 
						distanceParametersLower= NA,
						distanceParametersUpper= NA,
						useLambda=FALSE, lambdaLower = -6, lambdaUpper = 0, 
						#conditionalSimulation=FALSE, simulationReturnAll = FALSE, lambdaUpper = 0, 
            indefiniteMethod= "none", indefiniteType="PSD", indefiniteRepair=TRUE,
						scaling=FALSE,reinterpolate=FALSE) #todo always scale, remove scaling variable?
	con$algThetaControl[names(control$algThetaControl)] <- control$algThetaControl
	control$algThetaControl <- con$algThetaControl
	con[names(control)] <- control
 	control<-con

	#
	if(missing(distanceFunction))
		distanceFunction <- control$distanceFunction
		
	if(is.null(distanceFunction))
		stop("No distanceFunction passed to modelKriging.")
	
	if(length(distanceFunction)==1)
		control$combineDistances <- FALSE
	
	#check whether distance function has parameters
	useDistanceParameters=FALSE	
	if(!any(is.na(control$distanceParametersLower))&!any(is.na(control$distanceParametersUpper)))
		useDistanceParameters=TRUE 
		
	algThetaControl <- control$algThetaControl
	useLambda <- control$useLambda
	lambdaLower <- control$lambdaLower
	lambdaUpper <- control$lambdaUpper
	distanceParametersLower <- control$distanceParametersLower
	distanceParametersUpper <- control$distanceParametersUpper
	combineDistances <- control$combineDistances
	indefiniteMethod <- control$indefiniteMethod
	indefiniteType <- control$indefiniteType
	indefiniteRepair <- control$indefiniteRepair
	scaling <- control$scaling
	fcorr <- control$corr
 	fit <- control

	fit$useDistanceParameters <- useDistanceParameters
	
	if(!is.matrix(y)) #TODO why a matrix...
		y <- as.matrix(y)
	if(any(duplicated(x)) & !control$useLambda){ #duplicates HAVE to be removed, but duplicates for noisy problems are okay. 
		duplicates <- which(duplicated(x))
		x <- x[-duplicates]
		y <- as.matrix(y[-duplicates])
	}	
	
	fit$x <- x
	fit$y <- y

	n <- length(fit$x) #number of observations
	nd <- length(distanceFunction) # number of distance functions
	ntheta <- length(fit$lower) #number of theta parameters

	if(any(is.na(fit$lower))) #no lower bound means no theta parameter.
		ntheta=0
	
	#calculate distance matrix
	if(!useDistanceParameters){ #no distance parameters, can compute distance now. else: optimize and compute during MLE.
		ret <- modelKrigingDistanceCalculation(x,distanceFunction=distanceFunction,parameters=NA,
							distances=control$distances,scaling=scaling,combineDistances=combineDistances,indefiniteMethod=indefiniteMethod,
							indefiniteType=indefiniteType,indefiniteRepair=indefiniteRepair,lower=distanceParametersLower)
		fit[names(ret)] <- ret
		D <- fit$D
		fit$D <- NULL
  }
	
	
	if(is.null(control$userParameters)){ 
		# start point for theta and other model parameters + bounds:	
		res <- modelKrigingInit(fit$startTheta,fit$lower,fit$upper,
						useLambda,lambdaLower,lambdaUpper,
						combineDistances,nd,useDistanceParameters,
						distanceParametersLower,distanceParametersUpper)
		x0 <- res$x0
		lower <- res$lower
		upper <- res$upper
		
		# adapt tuning (MLE) budget to dimensionality of parameter space
		algThetaControl$funEvals <- algThetaControl$funEvals*length(x0)	
		if(combineDistances | nd==1){
			if(!useDistanceParameters){ # if distance function has no parameters (or default parameters are used:)
				res <- control$algTheta(x=x0,fun=modelKrigingLikelihood,lower=lower,upper=upper,
								control=algThetaControl,D=D,y=fit$y,useLambda=useLambda,corr=fcorr,
								indefiniteMethod=indefiniteMethod,indefiniteType=indefiniteType,indefiniteRepair=indefiniteRepair,returnLikelihoodOnly=TRUE,inverter="chol",ntheta=ntheta)	
				fit$distanceFunction <- distanceFunction
			}else{ # parameters of the distance function optimized during MLE
				res <- control$algTheta(x=x0,fun=modelKrigingParameterizedLikelihood,lower=lower,upper=upper,
								control=algThetaControl,xs=fit$x,ys=fit$y,useLambda=useLambda,corr=fcorr,
								indefiniteMethod=indefiniteMethod,indefiniteType=indefiniteType,indefiniteRepair=indefiniteRepair,returnLikelihoodOnly=TRUE,inverter="chol",
								distanceFunction=distanceFunction,combineDistances=combineDistances,distanceParametersLower=distanceParametersLower,ntheta=ntheta,scaling=scaling)	
				fit$distanceFunction <- distanceFunction #todo?
			}	
			nevals <- as.numeric(res$count[[1]])
		}else{
			res <- list()
			minlik=Inf
			minlikindex=1
			nevals <- 0
			for(i in 1:length(distanceFunction)){
				if(!useDistanceParameters){ # if distance function has no parameters (or default parameters are used:)
					res[[i]] <- control$algTheta(x=x0,fun=modelKrigingLikelihood,lower=lower,upper=upper,
							control=algThetaControl,D=D[[i]],y=fit$y,useLambda=useLambda,corr=fcorr,
							indefiniteMethod=indefiniteMethod,indefiniteType=indefiniteType,indefiniteRepair=indefiniteRepair,returnLikelihoodOnly=TRUE,inverter="chol",ntheta=ntheta)	
				}else{ # parameters of the distance function optimized during MLE
					res[[i]] <- control$algTheta(x=x0,fun=modelKrigingParameterizedLikelihood,lower=lower,upper=upper,
									control=algThetaControl,xs=fit$x,ys=fit$y,useLambda=useLambda,corr=fcorr,
									indefiniteMethod=indefiniteMethod,indefiniteType=indefiniteType,indefiniteRepair=indefiniteRepair,returnLikelihoodOnly=TRUE,inverter="chol",
									distanceFunction=distanceFunction[[i]],combineDistances=combineDistances,distanceParametersLower=distanceParametersLower,ntheta=ntheta,scaling=scaling)
				}
				if(res[[i]]$ybest < minlik){
					minlik <- res[[i]]$ybest
					minlikindex <- i
				}
				nevals <- nevals + as.numeric(res$count[[1]])
			}
			res <- res[[minlikindex]]
			fit$distanceFunction <- distanceFunction[[minlikindex]]
			fit$maximumDistance <- fit$maximumDistance[[minlikindex]]
			D <- D[[minlikindex]]
			fit$origD <- fit$origD[[minlikindex]]
			fit$isCNSD <- fit$isCNSD[[minlikindex]]
			fit$A <- fit$A[[minlikindex]]
			fit$distanceFunctionIndex <- minlikindex
			nd <- 1
		}	
		if(is.null(res$xbest)){
			res$xbest <- x0
		}
		Params <- res$xbest
	}else{
		Params <- control$userParameters
		nevals <- 0
		fit$distanceFunction <- distanceFunction
	}

	# extract model parameters:
	# kernel parameters (theta)	
	if(ntheta>0){
		fit$theta <- Params[1:ntheta]
	}
	# weights for each distance matrix (combination)
	if(combineDistances & nd>1){
		fit$distanceWeights <- 10^Params[ntheta+(1:nd)]
		nweights=nd
	}else{
		nweights=0#number of weight parameters
	}
	# lambda
	if(useLambda){
		fit$lambda <- 10^Params[ntheta+nweights+1]
	}else{
		fit$lambda <- 0
	}
	#distance function parameters
	if(useDistanceParameters){
		fit$distanceParameters <- Params[(ntheta+nweights+useLambda+1):length(Params)]
		res <- modelKrigingParameterizedLikelihood(Params,fit$x,fit$y,useLambda,fcorr,
											indefiniteMethod,indefiniteType,indefiniteRepair,
											returnLikelihoodOnly=FALSE,inverter="chol",
											distanceFunction=fit$distanceFunction,combineDistances=combineDistances,
											distanceParametersLower=distanceParametersLower,ntheta=ntheta,scaling=scaling
											)	#need to also return the correlation matrix and other elements of the model		
	}else{	
		res <- modelKrigingLikelihood(Params,D,fit$y,useLambda,fcorr,
											indefiniteMethod,indefiniteType,indefiniteRepair,
											returnLikelihoodOnly=FALSE,inverter="chol",ntheta=ntheta)	#need to also return the correlation matrix and other elements of the model
	}
	
	if(is.na(res$Psinv[1])){ #model building failed. no invertible correlation matrix was found. return NA fit
		stop("Building the Kriging model failed, no invertible correlation matrix was found. This may be due to the specific data-set or distance function used.")
	}
	
	if(useDistanceParameters | nd>1){
		fit$A <- res$A
		#D <- res$D
		if(!is.null(res$maximumDistance))
			fit$maximumDistance <- res$maximumDistance
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
	fit$origPsi <- res$origPsi
	fit$Psinv <- res$Psinv
	#fit$unrepairedPsi <- res$unrepairedPsi
	
	##precompute transformations
	if(indefiniteType=="PSD" & !fit$indefiniteRepair & fit$isIndefinite & any(indefiniteMethod==c("clip","flip","square","diffusion"))){ #RETRANSFORMATION OF THE SOLUTION ONLY  
		A <- res$U %*% diag(res$a) %*% t(res$U)
    fit$A <- A 
		#
    fit$Psinv <- t(A) %*% fit$Psinv #retransform the result	for prediction
		fit$PsinvA <- fit$Psinv %*% A #retransform the result	(for variance estimation only)
	}
	if(indefiniteType=="PSD" & (fit$indefiniteRepair %in% c(2,3,4)) & fit$isIndefinite & any(indefiniteMethod==c("clip","flip","square","diffusion"))){
		A <- res$U %*% diag(res$a) %*% t(res$U)
    fit$A <- A 
		fit$diagUnrepairedPsi <- diag(res$unrepairedPsi)			
		if(fit$indefiniteRepair==2){# for repair with nystroem only:	
			unrepairedPsinv <- try(chol2inv(chol(res$unrepairedPsi)), TRUE) 
			if(class(unrepairedPsinv)[1] == "try-error"){
				unrepairedPsinv <- ginv(res$unrepairedPsi) 
			}			
			fit$unrepairedAPsinvA <- t(A) %*% unrepairedPsinv %*% A 
		}
		fit$ADividedSqrtDiagPsi <- t(A) %*% diag(1/sqrt(diag(res$unrepairedPsi))) # divider for repair during prediction, including A
	}
	if(useLambda){ 
		PsiB <- res$Psi-diag(fit$lambda,n)+diag(.Machine$double.eps,n) 
		fit$SSQReint <- as.numeric((t(res$yMu)%*%res$Psinv%*%PsiB%*%res$Psinv%*%res$yMu)/n) #res is used intentionally, needs to be untransformed Psinv
		fit$PsinvReint <- try(chol2inv(chol(PsiB)), TRUE) 
		if(class(fit$PsinvReint)[1] == "try-error"){
			fit$PsinvReint <- ginv(PsiB) 
		}	
		#now apply same transformations as for non-reinterpolating matrices
		if(indefiniteType=="PSD" & fit$isIndefinite  & !fit$indefiniteRepair & any(indefiniteMethod==c("clip","flip","square","diffusion"))){ #RETRANSFORMATION OF THE SOLUTION ONLY  
      fit$PsinvReint <- t(A)%*%fit$PsinvReint %*% A #retransform
		} 
	}
	#	
	##
	fit$nevals <- nevals
	fit$like <- res$NegLnLike
  fit$predAll <- FALSE #todo : should be option
  fit$D <- D
	class(fit)<- "modelKriging"
	return(fit)
}

###################################################################################
#' Gaussian Kernel for Kriging
#'
#' @param D distance matrix
#' @param theta kernel parameter
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrGauss <- function(D,theta=0){
	theta <- 10^theta
	exp(-theta * D)
}

###################################################################################
#' Cubic Kernel for Kriging
#'
#' @param D distance matrix
#' @param theta kernel parameter
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrCubic <- function(D,theta=0){
	theta <- 10^theta
	Psi <- pmin(D * theta,1)
	1 - Psi^2 * (3 - 2*Psi)
}

###################################################################################
#' Linear Kernel for Kriging
#'
#' @param D distance matrix
#' @param theta kernel parameter
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrLinear <- function(D,theta=0){
	theta <- 10^theta
	pmax(1- D * theta,0)
}

###################################################################################
#' Spherical Kernel for Kriging
#'
#' @param D distance matrix
#' @param theta kernel parameter
#'
#' @return matrix (Psi)
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @export
#' @keywords internal
###################################################################################
fcorrSphere <- function(D,theta=0){
	theta <- 10^theta
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
#' @param lowerTheta lower boundary for theta values (log scale), the kernel parameters.
#' @param upperTheta upper boundary for theta values (log scale), the kernel parameters.
#' @param useLambda boolean, whether nugget effect (lambda) is used.
#' @param lambdaLower lower boundary for lambda (log scale).
#' @param lambdaUpper upper boundary for lambda (log scale).
#' @param combineDistances boolean, whether multiple distances are combined.
#' @param nd number of distance function.
#' @param distanceParameters whether the distance function parameters should be optimized
#' @param distanceParametersLower lower boundary for parameters of the distance function, default is \code{NA} which means there are no distance function parameters. If several distance functions are supplied, this should be a list of lower boundary vectors for each function.
#' @param distanceParametersUpper upper boundary for parameters of the distance function, default is \code{NA} which means there are no distance function parameters. If several distance functions are supplied, this should be a list of upper boundary vectors for each function.
#'
#' @return a list with elements \code{x0} (start guess), \code{lower} (lower bound), \code{upper} (upper bound).
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @keywords internal
###################################################################################
modelKrigingInit	<- function(startTheta=NULL,lowerTheta=NULL,upperTheta=NULL,useLambda, lambdaLower, lambdaUpper, combineDistances,nd,distanceParameters=F,distanceParametersLower=NA,distanceParametersUpper=NA){
	#ordering of the parameters:
	#first, the kernel function parameters. number:   - length(lowerTheta)
	#second, the weights for combining several distances (optional), number:   - 0|nd
	#fourth, lambda, regression constant, number:   - 0|1
	#fifth, distance parameters, number:   - 0|length(distanceParametersLower)
	if(any(is.na(lowerTheta))){ #NA bounds -> no parameter in the correlation function (at least none to be estimated)
		lowerTheta <- NULL
		upperTheta <- NULL
	}
	if(combineDistances){
		lowerTheta <- c(lowerTheta,rep(-8,nd))
		upperTheta <- c(upperTheta,rep(6,nd))
	}
	if(useLambda){
		#append regression constant lambda (nugget)
		lowerTheta <- c(lowerTheta,lambdaLower)
		upperTheta <- c(upperTheta,lambdaUpper)
	}	
  
  #parameters of the distance function
  if(distanceParameters){
    if(is.list(distanceParametersLower)){
      distanceParametersLower <- unlist(distanceParametersLower)
      distanceParametersUpper <- unlist(distanceParametersUpper)
    }
  	lowerTheta <- c(lowerTheta,distanceParametersLower)
		upperTheta <- c(upperTheta,distanceParametersUpper)  
  }
    
	#start value for theta	
	if(is.null(startTheta)){
		x0 <- lowerTheta + (upperTheta - lowerTheta)*0.5
	}else{
		#force x0 into bounds
		x0 <- pmin(x0,upperTheta)
		x0 <- pmax(x0,lowerTheta)
	}
	list(x0=x0,lower=lowerTheta,upper=upperTheta)
}


###################################################################################
#' Kriging: Distance Matrix Calculation
#'
#' Calculate and scale the distance matrix used in a Kriging model.
#' Include definiteness correction.
#' Not to be called directly.
#'
#' @param x list of samples in input space
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric.  It can also be a list of several distance functions. In this case, Maximum Likelihood Estimation (MLE) is used 
#'		to determine the most suited distance measure.
#'		The distance function may have additional parameters.
#' @param parameters parameters passed to the distance function as a vector.
#' @param distances precomputed distances, set to NA if not available.
#' @param scaling boolean, whether to scale the distance matrix.
#' @param combineDistances boolean, whether to combine the distances of different functions.
#' @param indefiniteMethod method for handling non-conditionally-definite matrices.
#' @param indefiniteType type of handling for non-conditionally-definite matrices.
#' @param indefiniteRepair whether to further repair other conditions (beside definiteness).
#' @param lower lower boundary for distance function parameters.
#'
#' @return a list with elements \code{D} (distance matrix), \code{maxD} (maximal distance for scaling purpose).
#'
#' @seealso \code{\link{modelKriging}}
#' 
#' @keywords internal
###################################################################################
modelKrigingDistanceCalculation <- function(x,distanceFunction,parameters=NA,
	distances,scaling,combineDistances,indefiniteMethod,indefiniteType,indefiniteRepair,lower){
	nd <- length(distanceFunction) # number of distance functions
	
	#calculate distance matrix
	if(nd==1){ #one distance function
		if(is.null(distances)){
			if(any(is.na(parameters))) #no parameters given
				D <-distanceMatrix(x,distanceFunction) 
			else #parameters are given
				D <-distanceMatrix(x,distanceFunction,parameters) 
		}else{
			D <- distances
		}
    maxD <- max(D) #maximum distance
    if(scaling){
			D <- D/maxD
    }    
	}else{ #multiple distance functions
		if(is.null(distances)){
			D <- list()
			maxD <- list()
			indices <- rep(1:nd,sapply(lower,length)) #indices assigning each parameter to a distances function
			for(i in 1:nd){
				if(any(is.na(parameters))) #no parameters given
					D[[i]] <-distanceMatrix(x,distanceFunction[[i]]) 
				else
					D[[i]] <-distanceMatrix(x,distanceFunction[[i]],parameters[indices==i])  	
        maxD[[i]] <- max(D[[i]]) #maximum distance
        if(scaling){
          D[[i]] <- D[[i]]/maxD[[i]]
        } 
			}
		}else{
			D <- distances
      maxD <- list()
			for(i in 1:nd){
        maxD[[i]] <- max(D[[i]]) #maximum distance
        if(scaling){
          D[[i]] <- D[[i]]/maxD[[i]]
        } 
			}    
		}
	}
	
	# Fix Definiteness (NSDness, CNSDness) of the provided distance matrix/matrices	
	origD <- D 
	A <- NA
	isCNSD <- NA
	matNoRep <- NA
	if(nd==1){#in case of one distance function
		ret <- correctionDistanceMatrix(D,indefiniteType,indefiniteMethod,indefiniteRepair)
		D <- ret$mat
		isCNSD <- ret$isCNSD
		A <- ret$A	
		matNoRep <- ret$matNoRep		
	}else if(!combineDistances){ #in case of multiple distances, which are not combined (but chosen from):
		isCNSD <- list()
		A <- list()	
		matNoRep <- list()
		for(i in 1:nd){
			ret <- correctionDistanceMatrix(D[[i]],indefiniteType,indefiniteMethod,indefiniteRepair)
			matNoRep[[i]] <- ret$matNoRep	
			D[[i]] <- ret$mat
			isCNSD[[i]] <- ret$isCNSD
			A[[i]] <- ret$A	
		}
	}		
	list(maximumDistance=maxD,D=D,origD=origD,A=A,isCNSD=isCNSD,matNoRep=matNoRep)
}