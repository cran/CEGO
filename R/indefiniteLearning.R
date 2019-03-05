###################################################################################
#' Correcting Conditional Negative Semi-Definiteness
#'
#' Correcting, e.g., a distance matrix with chosen methods so that it becomes a CNSD matrix.
#'
#' @param mat symmetric matrix, which should be at least of size 3x3
#' @param method string that specifies method for correction: spectrum clip \code{"clip"}, spectrum flip \code{"flip"}, nearest definite matrix \code{"near"}, spectrum square\code{"square"}, spectrum diffusion \code{"diffusion"}.
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' @return the corrected CNSD matrix
#'
#' @seealso \code{\link{modelKriging}}
#' @export
#' @examples
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' is.CNSD(D) #matrix should not be CNSD
#' D <- correctionCNSD(D)
#' is.CNSD(D) #matrix should now be CNSD
#' D
#' # note: to fix the negative distances, use repairConditionsDistanceMatrix. 
#' # Or else, use correctionDistanceMatrix.
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
###################################################################################
correctionCNSD <- function(mat,method="flip",tol=1e-8){
	n <- nrow(mat)
	v <- cbind(c(rep(1,n-1),1+sqrt(n)))
	Q <- diag(n) - 2 * (1/as.numeric(crossprod(v,v))) * tcrossprod(v,v)
	Fq <- Q %*% -mat %*% Q
	F1 <- Fq[-n,-n]        
	Fq[-n,-n] <- correctionDefinite(F1,type="PSD",method=method,tol)$mat #see nearCNSD, or nearest euclidean matrix. just one step.
	-Q %*% Fq %*% Q
}

###################################################################################
#' Correcting Definiteness of a Matrix
#'
#' Correcting a (possibly indefinite) symmetric matrix with chosen approach so that it will have desired definiteness type: positive or negative semi-definite (PSD, NSD).
#'
#' @param mat symmetric matrix
#' @param type string that specifies type of correction: \code{"PSD"},\code{"NSD"} to enforce PSD or NSD matrices respectively.
#' @param method string that specifies method for correction: spectrum clip \code{"clip"}, spectrum flip \code{"flip"}, nearest definite matrix \code{"near"}, spectrum square\code{"square"}, spectrum diffusion \code{"diffusion"}.
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' @return list with
#' \describe{ 
#'  \item{\code{mat}}{ corrected matrix}
#'  \item{\code{isIndefinite}}{ boolean, whether original matrix was indefinite}
#'  \item{\code{lambda}}{ the eigenvalues of the original matrix}
#'  \item{\code{lambdanew}}{ the eigenvalues of the corrected matrix }
#'  \item{\code{U}}{ the matrix of eigenvectors}
#'  \item{\code{a}}{ the transformation vector}
#' }
#'
#' @seealso \code{\link{modelKriging}}
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @export
#' @examples
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' is.NSD(D) #matrix should not be CNSD
#' D <- correctionDefinite(D,type="NSD")$mat
#' is.NSD(D) #matrix should now be CNSD
#' # different example: PSD kernel
#' D <- distanceMatrix(x,distancePermutationInsert)
#' K <- exp(-0.01*D)
#' is.PSD(K)
#' K <- correctionDefinite(K,type="PSD")$mat
#' is.PSD(K)
###################################################################################
correctionDefinite <- function(mat,type='PSD',method="flip",tol=1e-8){
  U <- NA
	isDefinite <- NA
  if(method != "none"){	
    if(type=="NSD")
      defSign <- -1
    else if(type=="PSD")
      defSign <- 1
		eig <- eigen(mat,symmetric=T)
    defEig <- defSign * eig$values
		U <- eig$vectors
		a <- rep(1,nrow(U))
		isDefinite <- min(defEig) >= -tol
		if(!isDefinite){ # only adapt if actually needed, else use default MLE algorithm.
			if(method=="clip"){ #or denoise in wu2005
				sel <- defEig>=tol 
				m <- sum(sel)
				a <- as.numeric(sel) 
				mat <- U[,sel,drop=F] %*% diag(eig$values[sel],m,m) %*% t(U[,sel,drop=F])   
			}else if(method=="flip"){ # see chen and wu
				sel <- defEig>=tol | defEig<=-tol 
				m <- sum(sel)
				a <- sign(defEig)
				mat <- U[,sel,drop=F] %*% diag(a[sel]*eig$values[sel],m,m)%*% t(U[,sel,drop=F])    
			}else if(method=="square"){ # mat(mat^T)
				#mat <- mat %*% t(mat) 
				sel <- defEig>=tol | defEig<=-tol
				m <- sum(sel)
				a <- defEig
				mat <- U[,sel,drop=F] %*% diag(a[sel]*eig$values[sel],m,m)%*% t(U[,sel,drop=F])      
			}else if(method=="diffusion"){ # expm(mat)
				sel <- defEig>=tol | defEig<=-tol
				m <- sum(sel)
				a <- defSign * exp(defEig) / eig$values
				mat <- U[,sel,drop=F] %*% diag(defSign * exp(defEig[sel]),m,m)%*% t(U[,sel,drop=F])    
			}else if(method=="near"){
				pd <- nearPD(defSign * mat, eig.tol = tol, conv.tol = tol,corr=TRUE, do2eigen=FALSE,keepDiag=FALSE,conv.norm.type="F") #corr=T forces diagonal 1, do2eigen should not be used! ruins results., the norm type may affect speed, chosen type "F" is in line with higham2002 
				mat <- defSign * as.matrix(pd$mat)
			}
		}	
		return(list(a=a,U=U,lambda=eig$values,lambdanew=a*eig$values,isDefinite=isDefinite,mat=mat))
	}else{
		return(NA)
	}
}

###################################################################################
#' Correction of a Distance Matrix
#'
#' Convert (possibly non-euclidean or non-metric) distance matrix with chosen approach so that it becomes a CNSD matrix.
#' Optionally, the resulting matrix is enforced to have positive elements and zero diagonal, with the \code{repair} parameter.
#' Essentially, this is a combination of functions \code{\link{correctionDefinite}} or \code{\link{correctionCNSD}} with \code{\link{repairConditionsDistanceMatrix}}.
#'
#' @param mat symmetric distance matrix
#' @param type string that specifies type of correction: \code{"CNSD"},\code{"NSD"} to enforce CNSD or NSD matrices respectively.
#' @param method string that specifies method for correction: spectrum clip \code{"clip"}, spectrum flip \code{"flip"}, nearest definite matrix \code{"near"}, spectrum square\code{"square"}, spectrum diffusion \code{"diffusion"}, feature embedding \code{"feature"}.
#' @param repair boolean, whether or not to use condition repair, so that elements are positive, and diagonal is zero.
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' @return list with corrected distance matrix \code{mat}, \code{isCNSD} (boolean, whether original matrix was CNSD) and transformation matrix \code{A}.
#'
#' @seealso \code{\link{correctionDefinite}},\code{\link{correctionCNSD}},\code{\link{repairConditionsDistanceMatrix}}
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @export
#' @examples
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' is.CNSD(D) #matrix should not be CNSD
#' D <- correctionDistanceMatrix(D)$mat
#' is.CNSD(D) #matrix should now be CNSD
#' D
###################################################################################
correctionDistanceMatrix <- function(mat,type="NSD",method="flip",repair=TRUE,tol=1e-8){
	isCNSD <- NA
	A <- NA
	matNoRep <- NA
  if((type=="NSD" | type=="CNSD") & any(method==c("clip","flip","near","square","diffusion","feature"))){
    isCNSD <- is.CNSD(mat,tol=tol) # check if definite 
		A <- diag(nrow(mat))
    if(!isCNSD){# mat is not CNSD, needs correction
			if(type=="NSD"){
        ret <- correctionDefinite(mat,type="NSD",method=method,tol=tol)
				#resulting, transformed matrix
        mat <- ret$mat
        #transformation matrix for new data (predict):
        A <- ret$U %*% diag(ret$a) %*% t(ret$U)
				if(repair){ # fix diagonal and range of values
					if(repair>1)
						matNoRep <- mat # mat before repair. only needed for Nystroem approximated repair.
					mat <- repairConditionsDistanceMatrix(mat)
        }
      }else if(type=="CNSD"){
        if(method == "near"){
          mat <- nearCNSD(mat,eig.tol=tol)$mat
				}else if(method=="feature"){
					x <- split(mat,seq(nrow(mat))) #each distance vector in the distance matrix is now a feature vector
					mat <- distanceMatrix(x,distanceRealEuclidean)  #TODO options for other surrogate distances?
				}else{
					mat <- correctionCNSD(mat,method=method,tol=tol)
					if(repair){ # fix diagonal and range of values
						if(repair>1)
							matNoRep <- mat # mat before repair. only needed for Nystroem approximated repair.
						mat <- repairConditionsDistanceMatrix(mat)
          }
        }
      }
    }
  }
	return(list(mat=mat,isCNSD=isCNSD,A=A,matNoRep=matNoRep))
}  

###################################################################################
#' Correction of a Kernel (Correlation) Matrix
#'
#' Convert a non-PSD kernel matrix with chosen approach so that it becomes a PSD matrix.
#' Optionally, the resulting matrix is enforced to have values between -1 and 1 and a diagonal of 1s, with the \code{repair} parameter.
#' That means, it is (optionally) converted to a valid correlation matrix.
#' Essentially, this is a combination of \code{\link{correctionDefinite}} with \code{\link{repairConditionsCorrelationMatrix}}.
#'
#' @param mat symmetric kernel matrix
#' @param method string that specifies method for correction: spectrum clip \code{"clip"}, spectrum flip \code{"flip"}, nearest definite matrix \code{"near"}, spectrum square\code{"square"}, spectrum diffusion \code{"diffusion"}.
#' @param repair boolean, whether or not to use condition repair, so that elements between -1 and 1, and the diagonal values are 1.
#' @param tol torelance value. Eigenvalues between \code{-tol} and \code{tol} are assumed to be zero.
#'
#' @return list with corrected kernel matrix \code{mat}, \code{isPSD} (boolean, whether original matrix was PSD), transformation matrix \code{A},
#' the matrix of eigenvectors (\code{U}) and the transformation vector (\code{a})
#'
#' @seealso \code{\link{correctionDefinite}}, \code{\link{repairConditionsCorrelationMatrix}}
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @export
#' @examples
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' K <- exp(-0.01*D)
#' is.PSD(K) #matrix should not be PSD
#' K <- correctionKernelMatrix(K)$mat
#' is.PSD(K) #matrix should now be CNSD
#' K
###################################################################################
correctionKernelMatrix <- function(mat,method="flip",repair=TRUE,tol=1e-8){
	isPSD <- NA
	A <- diag(nrow(mat))
	a <- NA
	U <- NA
	matNoRep <- NA
  if(any(method==c("clip","flip","near","square","diffusion"))){
    isPSD <- is.PSD(mat,tol=tol) # check if definite 
    if(!isPSD){# mat is not PSD, needs correction
			ret <- correctionDefinite(mat,type="PSD",method=method,tol=tol)
			#resulting, transformed matrix
			mat <- ret$mat
			a <- ret$a
			U <- ret$U
			isPSD <- ret$isDefinite
			#transformation matrix for new data (predict):
			A <- ret$U %*% diag(ret$a) %*% t(ret$U)			
			if(repair){ # fix diagonal and range of values
				if(repair>1)
					matNoRep <- mat # mat before repair. only needed for Nystroem approximated repair.
				mat <- repairConditionsCorrelationMatrix(mat)
      }
    }
  }
	return(list(mat=mat,matNoRep=matNoRep,isPSD=isPSD,A=A,a=a,U=U))
}  


###################################################################################
#' Repair Conditions of a Distance Matrix
#'
#' This function repairs distance matrices, so that the following two properties are ensured:
#' The distance values should be non-zero and the diagonal should be zero.
#' Other properties (conditionally negative semi-definitene (CNSD), symmetric) are
#' assumed to be given.
#'
#' @param mat symmetric, CNSD distance matrix. If your matrix is not CNSD, use \code{\link{correctionCNSD}} first. Or use \code{\link{correctionDistanceMatrix}}.
#'
#' @return repaired distance matrix
#'
#' @seealso \code{\link{correctionDefinite}}, \code{\link{correctionDistanceMatrix}}, \code{\link{correctionKernelMatrix}}, \code{\link{correctionCNSD}}, \code{\link{repairConditionsCorrelationMatrix}}
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @export
#' @examples
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' D <- correctionCNSD(D)
#' D
#' D <- repairConditionsDistanceMatrix(D)
#' D
###################################################################################
repairConditionsDistanceMatrix <- function(mat){ 
	n <- nrow(mat)
	eps <- sqrt(.Machine$double.eps)
	if(sum(abs(diag(mat)))>eps | (min(mat) < -eps)){ #if diagonal values are non zero, or if negative distance
		mat <- -mat #make cpsd kernel. (mat HAS TO be cnsd)
		Kaa <- matrix(diag(mat),n,n)
		mat <- Kaa + t(Kaa) - 2*mat #convert to valid distance (proven, because (c)psd.) 
	}
	mat	
}

###################################################################################
#' Repair Conditions of a Correlation Matrix
#'
#' This function repairs correlation matrices, so that the following two properties are ensured:
#' The correlations values should be between -1 and 1, and the diagonal values should be one.
#'
#' @param mat symmetric, PSD distance matrix. If your matrix is not CNSD, use \code{\link{correctionDefinite}} first. Or use \code{\link{correctionKernelMatrix}}.
#'
#' @return repaired correlation matrix
#'
#' @seealso \code{\link{correctionDefinite}}, \code{\link{correctionDistanceMatrix}}, \code{\link{correctionKernelMatrix}}, \code{\link{correctionCNSD}}, \code{\link{repairConditionsDistanceMatrix}}
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @export
#' @examples
#' x <- list(c(2,1,4,3),c(2,4,3,1),c(4,2,1,3),c(4,3,2,1),c(1,4,3,2))
#' D <- distanceMatrix(x,distancePermutationInsert)
#' K <- exp(-0.01*D)
#' K <- correctionDefinite(K,type="PSD")$mat
#' K
#' K <- repairConditionsCorrelationMatrix(K)
###################################################################################
repairConditionsCorrelationMatrix <- function(mat){ #dg: diag 1 or diag 0
	s <- diag(1/sqrt(diag(mat)))
	mat <- s %*% mat %*% s
}

###################################################################################
#' Augmented Distance Correction
#'
#' Correct new (test) distances, via correcting the augmented distance matrix. Internal use only.
#'
#' @param d new distance vector
#' @param object a modelKriging fit
#' @param x new samples (belonging to distances d)
#'
#' @return vector of augmented, corrected distances
#'
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @keywords internal
###################################################################################
correctionAugmentedDistanceVector <- function(d,object,x){ 
  if(is.vector(d))
    d <- matrix(d,1)
	D <- object$origD
	
	if(is.list(object$distanceFunction)){ #in case of multiple distance matrices 
		dself <- list()
		if(object$useDistanceParameters){
			indices <- rep(1:length(object$distanceFunction),sapply(object$distanceParametersLower,length)) #indices assigning each parameter to a distances function
		}
		for(i in 1:length(object$distanceFunction)){
			if(!object$useDistanceParameters){
				dself[[i]] <- distanceMatrix(x,object$distanceFunction[[i]]) 
			}else{
				dself[[i]] <- distanceMatrix(x,object$distanceFunction[[i]],object$distanceParameters[indices==i])			
			}	
			if(object$scaling){
				dself[[i]] <- dself[[i]]/object$maximumDistance[[i]]
			}
		}
		dself <- Reduce("+",mapply("*",dself,object$distanceWeights,SIMPLIFY=FALSE)) #weight each matrix by corresponding theta value, and compute sum of the matrices
	}else{
		if(!object$useDistanceParameters)
			dself <- distanceMatrix(x,object$distanceFunction) 
		else	
			dself <- distanceMatrix(x,object$distanceFunction,object$distanceParameters)			
		if(object$scaling){
			dself <- dself/object$maximumDistance
		}
	}	
	
	daug <- cbind(d,dself)
	Daug <- rbind(D,d)
	Daug <- cbind(Daug,t(daug))
	## Fix Definiteness (NSDness, CNSDness) of the provided distance matrix
	Daugtransformed <- correctionDistanceMatrix(Daug,object$indefiniteType,object$indefiniteMethod,object$indefiniteRepair)$mat
	## extract only the new values
	dnewtransformed <- Daugtransformed[(nrow(D)+1):nrow(Daugtransformed),1:ncol(D),drop=FALSE]
	## return
	dnewtransformed
}

###################################################################################
#' Augmented Kernel Correction
#'
#' Correct new (test) kernel values, via correcting the augmented kernel matrix. Internal use only.
#'
#' @param k new kernel value vector
#' @param object a modelKriging fit
#' @param x new samples (belonging to kernel values k)
#'
#' @return vector of augmented, corrected kernel values
#'
#' @references Martin Zaefferer and Thomas Bartz-Beielstein. (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @keywords internal
###################################################################################
correctionAugmentedKernelVector <- function(k,object,x){ #todo: tolerances! here and above
  if(is.vector(k))
    k <- matrix(k,1)
	K <- object$origPsi 
	if(is.list(object$distanceFunction)){ #in case of multiple distance matrices 
		dself <- list()
		if(object$useDistanceParameters){
			indices <- rep(1:length(object$distanceFunction),sapply(object$distanceParametersLower,length)) #indices assigning each parameter to a distances function
		}
		for(i in 1:length(object$distanceFunction)){
			if(!object$useDistanceParameters){
				dself[[i]] <- distanceMatrix(x,object$distanceFunction[[i]]) 
			}else{
				dself[[i]] <- distanceMatrix(x,object$distanceFunction[[i]],object$distanceParameters[indices==i])			
			}	
			if(object$scaling){
				dself[[i]] <- dself[[i]]/object$maximumDistance[[i]]
			}
		}
		dself <- Reduce("+",mapply("*",dself,object$distanceWeights,SIMPLIFY=FALSE)) #weight each matrix by corresponding theta value, and compute sum of the matrices
	}else{
		if(!object$useDistanceParameters)
			dself <- distanceMatrix(x,object$distanceFunction) 
		else	
			dself <- distanceMatrix(x,object$distanceFunction,object$distanceParameters)
		if(object$scaling){
			dself <- dself/object$maximumDistance
		}
	}
	if(is.null(object$theta)) #corr function has no parameters
		kself <- object$corr(dself)
	else
		kself <- object$corr(dself,object$theta)	
	kaug <- cbind(k,kself)
	Kaug <- rbind(K,k)
	Kaug <- cbind(Kaug,t(kaug))
	## Fix Definiteness (PNSDness) of the provided kernel matrix
	#Kaugtransformed <- correctionDefinite(Kaug,"PSD",object$indefiniteMethod,object$a)$mat
	#Kaugtransformed <- repairConditionsCorrelationMatrix(Kaugtransformed) 
  Kaugtransformed <- correctionKernelMatrix(Kaug,object$indefiniteMethod,object$indefiniteRepair)$mat
				## The following would only be needed if the whole matrix is of interest
				#Kaugtransformed <- Kaugtransformed + diag(object$lambda,nrow(Kaugtransformed))
	## extract only the new values
	knewtransformed <- Kaugtransformed[(nrow(K)+1):nrow(Kaugtransformed),1:ncol(K),drop=FALSE]
	
	#NOTE: lambda is not added to diagonal of Kaugtransformed, because this affects only the diagonal
	# which is not part of the returned vector
	
	## return
	knewtransformed
}	
