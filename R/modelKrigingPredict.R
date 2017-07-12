###################################################################################
#' Kriging Prediction
#' 
#' Predict with a model fit resulting from \code{\link{modelKriging}}.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{modelKriging}.
#' @param x list of samples to be predicted
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$predAll}\cr
#' TRUE: list with function value (mean) \code{object$y} and uncertainty estimate \code{object$s} (standard deviation)\cr
#' FALSE:\code{object$y}only
#'
#' @seealso \code{\link{modelKriging}}
#' @seealso \code{\link{simulate.modelKriging}}
#' @export
###################################################################################
predict.modelKriging <- function(object,x,...){ 
	ret <- modelKrigingInternalPredictor(object,x)
	psi <- ret$psi
	## return value:
	res <- list(y=ret$y)
	##########################################################################
	if (object$predAll){
		Psinv <- object$Psinv 
		lambda <- object$lambda
		SigmaSqr <- object$SSQ	
		if(object$indefiniteType=="PSD" & any(object$indefiniteMethod==c("clip","flip","square","diffusion"))){ 
			if(object$isIndefinite){
				if(!object$indefiniteRepair){
					Psinv <- object$PsinvA
				}
			}
		}
		if(object$reinterpolate & lambda > 0){
			SigmaSqr <- object$SSQReint	
			Psinv <- object$PsinvReint 
			lambda <- 0
		}
		# Psinv / PsinvReint has t(A)%*%Psi%*%A included already, if necessary. else, transformation has already been done performed for psi
		SSqr <- SigmaSqr*(1+lambda-diag(psi%*%Psinv%*%t(psi))) 
		s <- sqrt(abs(SSqr))
		res$s <- as.numeric(s) #return value
	}
	res
}



###################################################################################
#' Kriging Simulation
#' 
#' (Conditional) Simulate at given locations, with a model fit resulting from \code{\link{modelKriging}}.
#' In contrast to prediction or estimation, the goal is to reproduce the covariance 
#' structure, rather than the data itself. Note, that the conditional simulation 
#' also reproduces the training data, but
#' has a two times larger error than the Kriging predictor.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{modelKriging}.
#' @param nsim number of simulations
#' @param seed random number generator seed. Defaults to NA, in which case no seed is set
#' @param xsim list of samples in input space, to be simulated
#' @param conditionalSimulation  logical, if set to TRUE (default), the simulation is conditioned with the training data of the Kriging model.
#' Else, the simulation is non-conditional.
#' @param returnAll if set to TRUE, a list with the simulated values (y) and the corresponding covariance matrix (covar)
#' of the simulated samples is returned. 
#' @param ... further arguments, not used
#'
#' @return Returned value depends on the setting of \code{object$simulationReturnAll}
#'
#' @references N. A. Cressie. Statistics for Spatial Data. JOHN WILEY & SONS INC, 1993.
#' @references C. Lantuejoul. Geostatistical Simulation - Models and Algorithms. Springer-Verlag Berlin Heidelberg, 2002.
#'
#' @seealso \code{\link{modelKriging}}, \code{\link{predict.modelKriging}}
#' @export
###################################################################################
simulate.modelKriging <- function(object,nsim=1,seed=NA,xsim,conditionalSimulation=TRUE,returnAll=FALSE,...){
  if (!is.na(seed)){
	  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
			runif(1)
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }	
	#
  len <- length(xsim) #number of simulated samples (points  
  noise <- matrix(rnorm(len*nsim),len, nsim)
  #
	res <- computeCorrelationMatrix(object,xsim)
	covar <- res$psi
	#
  if(conditionalSimulation){
		ret <- modelKrigingInternalPredictor(object,xsim)
		y <- ret$y
		psi <- ret$psi
		
    covarDifference <- covar - psi %*% object$Psinv %*% t(psi)
    eigv <- eigen(object$SSQ *covarDifference,symmetric=T) #eigen decomposition
    covarDecomposed <- eigv$vectors %*% diag(sqrt(abs(eigv$values))) %*% eigv$vectors
    ysim <- covarDecomposed %*% noise
    
    #and the following adds the simulation part to the predictor
    y <- matrix(y,len,nsim) + ysim	
  }else{
    eigv <- eigen(object$SSQ *covar,symmetric=T) #eigen decomposition
    covarDecomposed <- eigv$vectors %*% diag(sqrt(abs(eigv$values))) %*% eigv$vectors
    y <- object$mu + covarDecomposed %*% noise
  }
	res$y <- y
  if(returnAll)
    return(res)	
  else
    return(y)
}


###################################################################################
#' Kriging Prediction (internal)
#' 
#' Predict with a model fit resulting from \code{\link{modelKriging}}.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{modelKriging}.
#' @param x list of samples to be predicted
#' @param ... further arguments, not used
#'
#' @return returns a list with:
#' \describe{
#' \item{\code{y}}{predicted values}
#' \item{\code{psi}}{correlations between x and training data}
#' }
#'
#' @seealso \code{\link{simulate.modelKriging}}
#' @seealso \code{\link{predict.modelKriging}}
#' @keywords internal
###################################################################################
modelKrigingInternalPredictor <- function(object,x){
	if(!is.list(x))x<-list(x)
	xo <- object$x
	
	Psinv <- object$Psinv 
	n <- length(xo)
	#one <- rep(1,n)
	mu <- object$mu
	yMu <- object$yMu	
	psi <- matrix(1,length(x),n)
	fundist <- object$distanceFunction

	if(is.list(fundist)){ # multiple distance functions to be combined
		psi <- replicate(length(fundist),psi,simplify=FALSE)
		if(object$useDistanceParameters){
			indices <- rep(1:length(fundist),sapply(object$distanceParametersLower,length)) #indices assigning each parameter to a distances function
		}
		for(j in 1:length(fundist)){
			for (i in 1:n){
				if(!object$useDistanceParameters){
					psi[[j]][,i] <- distanceVector(xo[[i]],x,fundist[[j]])				
				}else{
					psi[[j]][,i] <- distanceVector(xo[[i]],x,fundist[[j]],object$distanceParameters[indices==j])			
				}
			}				
			if(object$scaling){
				psi[[j]] <- psi[[j]]/object$maximumDistance[[j]]
			}
		}
		psi <-  Reduce("+",mapply("*",psi,object$distanceWeights,SIMPLIFY=FALSE)) #combine result by weighted sum
	}else{ #only one distance function
		for (i in 1:n){			
			if(!object$useDistanceParameters){
				psi[,i] <- distanceVector(xo[[i]],x,fundist)	
			}else{
				psi[,i] <- distanceVector(xo[[i]],x,fundist,object$distanceParameters)			
			}
		}		
		if(object$scaling){
			psi <- psi/object$maximumDistance
		}
  }
  #
	##
  #
  if(any(object$indefiniteMethod==c("clip","flip","near","square","diffusion"))){ #Distance matrix was not CNSD, and a suitable correction method was chosen.
		if(object$indefiniteType=="NSD" & !object$indefiniteRepair){ #no repair, NSD-correction: transformation can be used directly 		
			if(!object$isCNSD)
				psi <- psi %*% t(object$A)
    }
		if(object$indefiniteType=="CNSD" | (object$indefiniteType=="NSD" & object$indefiniteRepair)){ 
			if(!object$isCNSD)
				psi <- correctionAugmentedDistanceVector(psi,object,x) #in case of repair, CNSD-correction: retransform with augmented distance matrix
		}
  }	
	if((object$indefiniteType=="CNSD" | object$indefiniteType=="NSD") & object$indefiniteMethod=="feature"){ #distances as features
		if(!object$isCNSD){
			tempx <- split(psi,seq(nrow(psi)))
			for (i in 1:n)
				psi[,i] <- distanceVector(object$origD[i,],tempx,distanceRealEuclidean) #todo choice of distance
		}
	}
  #
	##
	#
	if(is.null(object$theta)) #corr function has no parameters
		psi <- object$corr(psi)
	else
		psi <- object$corr(psi,object$theta)
 
	if(object$indefiniteType=="PSD" & any(object$indefiniteMethod==c("clip","flip","near","square","diffusion"))){ 
    if(object$isIndefinite){
			if(!object$indefiniteRepair & any(object$indefiniteMethod==c("clip","flip","square","diffusion"))){
				#psi <- psi %*% t(object$A)  #This is already included in Psinv. do nothing. 
			}else{
				psi <- correctionAugmentedKernelVector(psi,object,x)
			}
    }
  }
  y <- as.numeric(psi%*%Psinv%*%yMu)+mu  #todo: Psinv%*%yMu can be precomputed for speedup
  list(y=y,psi=psi)
}


###################################################################################
#' Compute Correlation Matrix
#' 
#' Compute the correlation matrix of samples x, given the model object.
#'
#' @param object fit of the Kriging model (settings and parameters), of class \code{modelKriging}.
#' @param x list of samples / data
#'
#' @return the correlation matrix
#'
#' @seealso \code{\link{simulate.modelKriging}}
#' @seealso \code{\link{predict.modelKriging}}
#' @keywords internal
###################################################################################
computeCorrelationMatrix <- function(object,x){
	if(!is.list(x))x<-list(x)
	if(is.null(object$distanceParameters))
		object$distanceParameters <- NA
	ret <- modelKrigingDistanceCalculation(x,object$distanceFunction,parameters=object$distanceParameters,
					NULL,object$scaling,object$combineDistances,object$indefiniteMethod,object$indefiniteType,object$indefiniteRepair,object$distanceParametersLower)	
	psi <- ret$D
  #
	##
	#
	if(is.null(object$theta)) #corr function has no parameters
		psi <- object$corr(psi)
	else
		psi <- object$corr(psi,object$theta)
  #
	##
	# 
	ret$U <- NA
	ret$a <- NA
	ret$isIndefinite <- NA
	ret$origPsi <- NA
	if(object$indefiniteType=="PSD" & any(object$indefiniteMethod==c("clip","flip","near","square","diffusion"))){ 
    #psi <- correctionKernelMatrix(psi,object$indefiniteMethod,object$indefiniteRepair)$mat
		ret$origPsi <- psi
		ret2 <- correctionKernelMatrix(psi,object$indefiniteMethod,object$indefiniteRepair)
		ret$a <- ret2$a
		ret$U <- ret2$U
		ret$A <- ret2$A
		ret$isIndefinite <- !ret2$isPSD
		psi <- ret2$mat			
		
  }
  #
	##
	#	
	if(object$useLambda){
		psi <- psi + diag(object$lambda,length(x)) 
	}
  #
	##
	#	
  ret$psi <- psi
	ret
}