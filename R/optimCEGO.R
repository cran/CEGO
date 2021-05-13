#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Combinatorial Efficient Global Optimization
#' 
#' Model-based optimization for combinatorial or mixed problems. Based on measures of distance or dissimilarity.
#' 
#' @param x Optional initial design as a list. If NULL (default), \code{creationFunction} (in \code{control} list) is used to create initial design. 
#' If \code{x} has less individuals than specified by \code{control$evalInit}, \code{creationFunction} will fill up the design.
#' @param fun target function to be minimized
#' @param control (list), with the options of optimization and model building approaches employed:
#' \describe{
#' \item{\code{evalInit}}{ Number of initial evaluations (i.e., size of the initial design), integer, default is \code{2}}
#' \item{\code{vectorized}}{ Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE}
#' \item{\code{verbosity}}{ Level of text output during run. Defaults to 0, no output.}
#' \item{\code{plotting}}{ Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.}
#' \item{\code{targetY}}{ optimal value to be found, stopping criterion, default is \code{-Inf}}
#' \item{\code{budget}}{ maximum number of target function evaluations, default is \code{100}}
#' \item{\code{creationRetries}}{ When a model does not predict an actually improving solution, a random exploration step is performed. \code{creationRetries} solutions are created randomly. 
#' 		For each, distance to all known solutions is calculated. The minimum distance is recorded for each random solution. 
#' 		The random solution with maximal minimum distance is chosen doe be evaluated in the next iteration.}
#' \item{\code{model}}{ Model to be used as a surrogate of the target function. Default is "K" (Kriging). Also
#'		available are: "LM" (linear, distance-based model), "RBFN" Radial Basis Function Network.}
#' \item{\code{modelSettings}}{ List of settings for \code{model} building, passed on as the \code{control} argument to the model training functions \code{\link{modelKriging}}, \code{\link{modelLinear}}, \code{\link{modelRBFN}}.}
#' \item{\code{infill}}{ This parameter specifies a function to be used for the infill criterion (e.g., the default is expected improvement \code{infillExpectedImprovement}).
#' To use no specific infill criterion this has to be set to \code{NA}, in which case the prediction of the surrogate model is used. Infill criteria are only used with models that may provide some error estimate with predictions.}
#' \item{\code{optimizer}}{ Optimizer that finds the minimum of the surrogate model. Default is \code{\link{optimEA}}, an Evolutionary Algorithm.}
#' \item{\code{optimizerSettings}}{ List of settings (\code{control}) for the \code{optimizer} function.}
#' \item{\code{initialDesign}}{ Design function that generates the initial design. Default is \code{designMaxMinDist}, which creates a design that maximizes the minimum distance between points.}
#' \item{\code{initialDesignSettings}}{ List of settings (\code{control}) for the \code{initialDesign} function.}
#' \item{\code{creationFunction}}{ Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6}
#' \item{\code{distanceFunction}}{ distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are not a problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric. With the setting \code{control$model="K"} this can also be a list of different fitness functions.
#'    Default is Hamming distance for permutations: distancePermutationHamming.}
#' }
#'
#' @return a list:
#' \describe{
#' \item{\code{xbest}}{ best solution found}
#' \item{\code{ybest}}{ fitness of the best solution}
#' \item{\code{x}}{ history of all evaluated solutions}
#' \item{\code{y}}{ corresponding target function values f(x)}
#' \item{\code{fit}}{ model-fit created in the last iteration}
#' \item{\code{fpred}}{ prediction function created in the last iteration}
#' \item{\code{count}}{ number of performed target function evaluations}
#' \item{\code{message}}{ message string, giving information on termination reason}
#' \item{\code{convergence}}{ error/status code: \code{-1} for termination due 
#' to failed model building, \code{0} for termination due to depleted budget, 
#' \code{1} if attained objective value is equal to or below target (\code{control$targetY})}
#' }
#' 
#' @examples
#' seed <- 0
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(1:5,dF)
#' #start optimization
#' set.seed(seed)
#' res1 <- optimCEGO(,lF,list(
#'				creationFunction=cF,
#'				distanceFunction=dF,
#'				optimizerSettings=list(budget=100,popsize=10,
#'				mutationFunction=mF,recombinationFunction=rF),
#'		evalInit=5,budget=15,targetY=0,verbosity=1,model=modelKriging,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' set.seed(seed)
#' res2 <- optimCEGO(,lF,list(
#'				creationFunction=cF,
#'				distanceFunction=dF,
#'				optimizerSettings=list(budget=100,popsize=10,
#'				mutationFunction=mF,recombinationFunction=rF),
#'				evalInit=5,budget=15,targetY=0,verbosity=1,model=modelRBFN,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res1$xbest 
#' res2$xbest 
#'
#' @seealso \code{\link{modelKriging}}, \code{\link{modelLinear}}, \code{\link{modelRBFN}}, \code{\link{buildModel}}, \code{\link{optimEA}}
#' 
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#'
#' @export
###################################################################################
optimCEGO <- function(x=NULL,fun,control=list()){
	## default settings
	con<-list(evalInit = 2
			, vectorized = FALSE
			, verbosity = 0
			, plotting = FALSE
			, targetY = -Inf
      , budget = 100
			, creationRetries = 100
			, distanceFunction = distancePermutationHamming
			, creationFunction = solutionFunctionGeneratorPermutation(6)
			, infill= infillExpectedImprovement			
      , model = modelKriging
			, modelSettings= list()
      , optimizer = optimEA
			, optimizerSettings = list()
			, initialDesign = designMaxMinDist
			, archiveModelInfo = NULL #TODO document
			, initialDesignSettings = list())
	con[names(control)] <- control
	control<-con
	rm(con)
	count <- control$evalInit
	archiveModelInfo <- control$archiveModelInfo
  vectorized <- control$vectorized
  verbosity <- control$verbosity
  plotting <- control$plotting
	creationFunction <- control$creationFunction
	distanceFunction <- control$distanceFunction
	
	if(is.null(control$initialDesignSettings$distanceFunction))
		control$initialDesignSettings$distanceFunction <- distanceFunction
	
	## if target function not vectorized: vectorize with lapply
	fun #lazy load
	if(!vectorized) 
		fn <- function(x)unlist(lapply(x,fun))
	else 
		fn <- fun
	
	## Create main object of this function, which will also be the return value	
	res <- list(xbest=NA, ybest=NA, x=NA,y=NA,distances=NA,modelArchive=NA,count=count,convergence=0,message="")
	
  ## Termination information:
  msg <- "Termination message:"
  
	## Create initial design of experiment
	res$x <- control$initialDesign(x,creationFunction,count,control$initialDesignSettings)

	## Calculate distances between samples. If distance function has parameters: do not calculate.
	distanceHasParam <- FALSE
	if(is.function(distanceFunction)){
		if(length(distanceFunction)==1)
			distanceHasParam <- length(formalArgs(distanceFunction))>2
		else
			distanceHasParam <- any(sapply(sapply(distanceFunction,formalArgs,simplify=FALSE),length) > 2)
		if(!distanceHasParam)
			res$distances <- distanceMatrixWrapper(res$x,distanceFunction)
	}		

	
	## Evaluate initial population	
	res$y <- fn(res$x)

	## determine best
	indbest <- which.min(res$y)
	res$ybest <- res$y[[indbest]]
	res$xbest <- res$x[[indbest]]	
	
	## build initial model
	model <- buildModel(res,distanceFunction,control)
		
	## archive desired model information
	if(!is.null(archiveModelInfo)){
		res$modelArchive <- list()
		archiveIndex <- 1
		if(identical(model,NA)){
			res$modelArchive[[archiveIndex]] <- rep(NA,length(archiveModelInfo))
			names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
		}else{ #todo!			
			res$modelArchive[[archiveIndex]] <- model$fit[archiveModelInfo]
			names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
		}		
	}

	## check whether EI infill is used
	useEI <- is.function(control$infill)
	
	## main loop
	while((res$count < control$budget) & (res$ybest > control$targetY)){
		## Optimize the surrogate:
    if(!identical(model,NA)){ 
      optimres <- optimizeModel(res,creationFunction,model,control)				
      ## Handle duplicated new candidate solution	
      duplicate <- list(optimres$xbest) %in% res$x
			## Check whether next candidate solution is better than the best observed so far (based on model prediction)
      improved <- optimres$ybest < optimres$fpredbestKnownY
    }else{# model building failed, force termination
      msg <- paste(msg,"Model building failed, optimization stopped prematurely.")		
      warning("Model building failed in optimCEGO, optimization stopped prematurely.")
      res$convergence <- -1
      break;
    }
    
  	## Update evaluation counter
		res$count <- res$count+1
    
		if(!duplicate && ((improved || useEI))){ #exploitation step		#for a new individual to be accepted, it has to have a better predicted value than the prediction for the best known individual. One could also use the actual value of the best known individual, but this would deteriorate results in case of an offset.
			res$x[[res$count]] <- optimres$xbest		#NOTE; exploitation step and exploration is both the same when EI is used., thus the "||"
		}else{ #exploration step: no promising individual found, or promising individual is a duplicate -> create a new one randomly
			if(!distanceHasParam){
				designSize <- length(res$x)+1
				if(is.list(distanceFunction)) #in case of multiple distances
					dfun <- distanceFunction[[1]]
				else 	
					dfun <- distanceFunction
				xc <- designMaxMinDist(res$x,creationFunction,designSize,control=list(budget=control$creationRetries,distanceFunction=dfun))			
				res$x[[res$count]] <- xc[[designSize]] #this maximizes distance, but may still create a duplicate if max(min(dist))==0, e.g. if all randomly created individuals are duplicates of known solutions
			}else{
				res$x[[res$count]] <- optimres$xbest
			}
		}
		res$x <- removeDuplicates(res$x, creationFunction)
		
		## evaluate with real target function
		res$y <-  c(res$y,fn(res$x[res$count]))
			
		## Logging
		indbest <- which.min(res$y)
		res$ybest <- res$y[[indbest]]
		res$xbest <- res$x[[indbest]]
		
		## Plotting and Text output
    if(verbosity > 0)
      print(paste("Evaluations:",res$count,"    Quality:",res$ybest))
		if(plotting){
			plot(res$y,type="l",xlab="number of evaluations", ylab="y")
			abline(res$ybest,0,lty=2)
		}
	
		## Update the distance matrix				#TODO what if distance parameters?
		if(!distanceHasParam & is.function(distanceFunction))
			res$distances <- distanceMatrixUpdate(res$distances,res$x,distanceFunction)
	
		## Update surrogate model and prediction function:
		model <- buildModel(res,distanceFunction,control)
		
		## archive desired model information
		if(!is.null(archiveModelInfo)){
			archiveIndex <- archiveIndex+1
			if(identical(model,NA)){
				res$modelArchive[[archiveIndex]] <- rep(NA,length(archiveModelInfo))
				names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
			}else{
				res$modelArchive[[archiveIndex]] <- model$fit[archiveModelInfo]
				names(res$modelArchive[[archiveIndex]]) <- archiveModelInfo
			}
		}
	}
  
  #stopping criteria information for user:
	if(min(res$ybest,na.rm=TRUE) <= control$targetY) {
		msg <- paste(msg,"Successfully achieved target fitness.")		
    res$convergence <- 1
	}else if(res$count >= control$budget){ #budget exceeded
    msg <- paste(msg,"Target function evaluation budget depleted.")		
  }
  res$message <- msg
  res$distances <- NULL
	res #return
}


###################################################################################
#' Model building
#' 
#' Model building support function for optimCEGO. Should not be called directly.
#' 
#' @param x list of samples in input space
#' @param y  matrix, column vector of observations for each sample
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric.  In case Kriging is chosen, it can also be a list of several distance functions. In this case, MLE is used 
#'		to determine the most suited distance measure (see the last reference).
#' @param control list with options:
#' \describe{
#' \item{\code{model}}{ Model to be used as a surrogate of the target function. Default is "K" (Kriging). Also
#'		available are: "LM" (linear, distance-based model), "RBFN" Radial Basis Function Network.}
#' \item{\code{modelSettings}}{ List of settings for model building, passed on as the control argument to the model training functions \code{\link{modelKriging}}, \code{\link{modelLinear}}, \code{\link{modelRBFN}}.}
#' \item{\code{infill}}{ This parameter specifies a function to be used for the infill criterion (e.g., the default is expected improvement \code{infillExpectedImprovement}).
#' To use no specific infill criterion this has to be set to \code{NA}. Infill criteria are only used with models that may provide some error estimate with predictions.}
#' }
#'
#' @return a list:
#' \describe{
#' \item{\code{fit}}{ model-fit }
#' \item{\code{fpred}}{ prediction function}
#' }
#' 
#' @seealso \code{\link{optimCEGO}} 
#'
#' @keywords internal
#' @export
###################################################################################
buildModel <- function(res,distanceFunction,control){ 
	y <- res$y
	x <- res$x
	control$modelSettings$distances <- res$distances
	
	y #against lazy evaluation
	if(identical(control$model,"RBFN"))
		control$model <- modelRBFN
	else if(identical(control$model,"LM"))
		control$model <- modelLinear
	else if(identical(control$model,"K"))
		control$model <- modelKriging
	
	fit<-try(control$model(x,y,distanceFunction,control$modelSettings),TRUE)
	#fit <- control$model(x,y,distanceFunction,control$modelSettings)
	if(class(fit)[1] == "try-error"){
    #warning("Model building in optimCEGO failed.") #same warning given in optimCEGO function
    return(NA)
  }
	
	fit$predAll <- is.function(control$infill)
	fit
	
	if(is.function(control$infill)){
		fpred <- function(x){	
				res=predict(fit,x)
				control$infill(res$y,res$s,min(y))
			} 
	}else{
		fpred <- function(x){predict(fit,x)$y}
	}

	list(fit=fit,fpred=fpred)
}


###################################################################################
#' Optimize Surrogate Model
#'
#' Interface to the optimization of the surrogate model
#' 
#' @param res result state of the optimization process
#' @param creationFunction Function to create individuals/solutions in search space.
#' @param model result of the buildModel function
#' @param control list of settings, from optimCEGO
#'
#' @return result list of the optimizer
#' 
#' @seealso \code{\link{optimCEGO}} 
#'
#' @keywords internal
###################################################################################
optimizeModel <- function(res,creationFunction,model,control){ 
	if(identical(control$optimizer,"EA")){
		control$optimizer=optimEA
	}
	if(is.null(control$optimizerSettings$creationFunction))
		control$optimizerSettings$creationFunction <- creationFunction
	if(is.null(control$optimizerSettings$vectorized))
		control$optimizerSettings$vectorized <- TRUE		
	optimres <- control$optimizer(NULL,model$fpred,control$optimizerSettings)
	optimres$fpredbestKnownY <- model$fpred(list(res$xbest))
	optimres
}