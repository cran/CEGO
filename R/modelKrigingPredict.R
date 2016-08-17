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
#' @export
###################################################################################
predict.modelKriging <- function(object,x,...){ 
	if(!is.list(x))x<-list(x)
	xo <- object$x
	
	#
	##
	#
	if(object$indefiniteMethod=="feature"){ #distances as features
		xtransform <- list()
		for(j in 1:length(x)){
			xtransform[[j]] <- distanceVector(x[[j]],object$origx,object$origDistanceFunction) #compute "original" distance to training data
		}
		x <- xtransform
	}
	#
	##
  #
	
	theta <- object$theta 
	Psinv <- object$Psinv 
	n <- length(xo)
	#one <- rep(1,n)
	mu <- object$mu
	yMu <- object$yMu	
	psi <- matrix(1,length(x),n)
	fundist <-object$distanceFunction
	
	if(is.list(fundist)){ # multiple distance functions to be combined
		psi <- replicate(length(fundist),psi,simplify=FALSE)
		for(j in 1:length(fundist)){
			for (i in 1:n)
				psi[[j]][,i] <- distanceVector(xo[[i]],x,fundist[[j]])					
			if(object$scaling){
				psi[[j]] <- psi[[j]]/object$maximumDistance[[j]]
			}
		}
		psi <-  Reduce("+",mapply("*",psi,theta,SIMPLIFY=FALSE)) #combine result by weighted sum
		theta <- 1 #theta used in correlation function (corr) set to 1, because included in weighted sum.
	}else{ #only one distance function
		for (i in 1:n)
			psi[,i] <- distanceVector(xo[[i]],x,fundist)
				
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
				psi <- correctionAugmentedDistanceVector (psi,object,x) #in case of repair, CNSD-correction: retransform with augmented distance matrix
		}
  }	
  #
	##
	#
 	if(object$optimizeP){
		psi <- psi^object$p
	} 
	#
	psi <- object$corr(theta,psi)
 
	if(object$indefiniteType=="PSD" & any(object$indefiniteMethod==c("clip","flip","near","square","diffusion"))){ 
    if(object$isIndefinite){
			if(!object$indefiniteRepair & any(object$indefiniteMethod==c("clip","flip","square","diffusion"))){
				#psi <- psi %*% t(object$A)  #This is already included in Psinv. do nothing. 
			}else{
				psi <- correctionAugmentedKernelVector(psi,object,x)
			}
    }
  }
  f <- as.numeric(psi%*%Psinv%*%yMu)+mu  
  
	## return value:
	res <- list(y=f)
	##########################################################################
	if (object$predAll){
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
		res$s <- as.numeric(s)
	}
	res
}