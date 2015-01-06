#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Combinatorial Efficient Global Optimization
#' 
#' Model-based optimization for combinatorial or mixed problems. Based on measures of distance or dissimilarity.
#' 
#' @param par Optional initial design as a list. If NULL (default), creationFunction is used to create initial design. 
#' If par has less individuals than specified by \code{control$evalInit}, creationFunction will fill up the design.
#' @param fun target function to be minimized
#' @param creationFunction Function to create individuals/solutions in search space
#' @param distanceFunction a suitable distance function of type f(x1,x2), returning a scalar distance value, preferably between 0 and 1.
#'      Maximum distances larger 1 are no problem, but may yield scaling bias when different measures are compared.
#' 		Should be non-negative and symmetric. With the setting \code{control$model="K"} this can also be a list of different fitness functions.
#' @param control (list), with the options of optimization and model building approaches employed\cr
#' \code{evalInit} Number of initial evaluations (i.e., size of the initial design), integer, default is \code{2}\cr
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
#' \code{verbosity} Level of text output during run. Defaults to 0, no output.\cr
#' \code{plotting} Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.\cr
#' \code{targetY} optimal value to be found, stopping criterion, default is \code{-Inf}\cr
#' \code{evalBudget} maximum number of target function evaluations, default is \code{100}\cr
#' \code{creationRetries} When a model does not predict an actually improving solution, a random exploration step is performed. \code{creationRetries} solutions are created randomly. 
#' 		For each, distance to all known solutions is calculated. The minimum distance is recorded for each random solution. 
#' 		The random solution with maximal minimum distance is chosen doe be evaluated in the next iteration.\cr
#' \code{model} Model to be used as a surrogate of the target function. Default is "K" (Kriging). Also
#'		available are: "LM" (linear, distance-based model), "RBFN" Radial Basis Function Network.\cr
#' \code{modelSettings} List of settings for model building, passed on as the control argument to the model training functions \code{\link{combinatorialKriging}}, \code{\link{combinatorialLM}}, \code{\link{combinatorialRBFN}}.
#' \code{infill} This parameter specifies a function to be used for the infill criterion (e.g., the default is expected improvement \code{infillExpectedImprovement}).
#' To use no specific infill criterion this has to be set to \code{NA}. Infill criteria are only used with models that may provide some error estimate with predictions.\cr
#' \code{optimizer} Optimizer that finds the minimum of the surrogate model. Default is \code{"EA"} an Evolutionary Algorithm. No alternatives implemented yet.
#' \code{optimSettings} List of settings for the method to optimize the model. 
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{x} history of all evaluated solutions\cr
#' \code{y} corresponding target function values f(x)\cr
#' \code{fit} model-fit created in the last iteration\cr
#' \code{fpred} prediction function created in the last iteration\cr
#' \code{count} number of performed target function evaluations 
#' 
#' @examples
#' seed=0
#' glgseed=1
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(ref=1:5,distFun=dF)
#' #start optimization
#' set.seed(seed)
#' res1 <- optimCEGO(,lF,cF,dF,list(optimSettings=list(
#'				mutationFunction=mF,recombinationFunction=rF),
#'		evalInit=5,budget=15,targetY=0,verbosity=1,model="K"))
#' set.seed(seed)
#' res2 <- optimCEGO(,lF,cF,dF,list(optimSettings=list(
#'				mutationFunction=mF,recombinationFunction=rF),
#'		evalInit=5,budget=15,targetY=0,verbosity=1,model="RBFN"))
#' res1$xbest 
#' res2$xbest 
#'
#' @seealso \code{\link{combinatorialKriging}}, \code{\link{combinatorialLM}}, \code{\link{combinatorialRBFN}}, \code{\link{buildModel}}, \code{\link{optimEA}} 
#' 
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282 http://doi.acm.org/10.1145/2576768.2598282 
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#'
#' @export
###################################################################################
optimCEGO <- function(par=NULL,fun,creationFunction,distanceFunction,control=list()){
	## default settings
	con<-list(evalInit=2
			, vectorized=FALSE
			, verbosity= 0
			, plotting=FALSE
			, targetY = -Inf
      , budget=100
			, creationRetries=100
      , model = "K"
			, modelSettings= list()
			, infill= infillExpectedImprovement
      , optimizer = "EA"
			, optimSettings= list())
	con[(namc <- names(control))] <- control
	control<-con

	count=control$evalInit
  vectorized = control$vectorized
  verbosity = control$verbosity
  plotting = control$plotting
  budget = control$budget
	
	## Create initial design of experiment
	x <- initializeDesign(par,creationFunction,count)
	
	## Calculate distances between samples	
	if(length(distanceFunction)==1){ # in case of a single distance function (all models)
		control$modelSettings$distances <- distanceMatrix(x,distanceFunction)
	}else{	# in case of multiple distance functions (kriging only atm.)
		control$modelSettings$distances <- list()
		for(i in 1:length(distanceFunction)){
			control$modelSettings$distances[[i]] <- distanceMatrix(x,distanceFunction[[i]]) 
		}
	}
	##
		
	## Evaluate initial population	
	if(vectorized) 
		y <- fun(x)
	else
		y <- unlist(lapply(x,fun)) 
    
	##build initial model
	res <- buildModel(x,y,distanceFunction,control)
    fit <- res$fit
	fpred <- res$fpred
	best <- which.min(y)
	bestKnownY <- y[best] #best in initial population
	bestKnownX <- x[[best]] #best in initial population
	#print(x)
	#check whether EI infill is used
	useEI <- is.function(control$infill)
	while((count < budget) & (bestKnownY > control$targetY)){
		count=count+1
		
		#Optimize the surrogate:
		if(is.list(fit)){ #only if the model is present
			if(control$optimizer=="EA"){
				bestPredicted <- optimEA(NULL,fpred,creationFunction,control$optimSettings)
			}#else if(control$optimizer=="RS"){}		
			fpredbestKnownX <- fpred(bestKnownX) 
		}else{ #this case comes into play when the model building failed for some reason. this will result into a random new individual.
			bestPredicted <- list(xbest=creationFunction(),ybest=bestKnownY+1)
			fpredbestKnownX <- bestKnownY
		}
		
		duplicate = FALSE
		if(duplicated(append(x,list(bestPredicted$xbest)))[count])
			duplicate = TRUE
		
		#todo in the following line, if useEI is true, but the predicted best has EI of zero.... ?
		if(!duplicate && ((bestPredicted$ybest < fpredbestKnownX || useEI))){ #exploitation step		#for a new individual to be accepted, it has to have a better predicted value than the prediction for the best known individual. One could also use the actual value of the best known individual, but this would deteriorate results in case of an offset.
			x[[count]] <- bestPredicted$xbest		#NOTE; exploitation step and exploration is both the same when EI is used., thus the "||"
		}else{ #exploration step: no promising individual found, or promising individual is a duplicate -> create a new one randomly
			#xplore <- creationFunction()
			#while(any(duplicated(append(x,list(xplore))))) #todo if a "simple" problem gets oversampled, the algorithm will get stuck here. solution: calculate number of possible combinations? or something similar?
				#xplore <- creationFunction()
			#x[[count]] <- xplore	
			xc <- list() #candidates
			distlist <- NULL
			for (i in 1:control$creationRetries){ #create candidates, to get one new individual that is maximally different from the existing.
				xc[[i]] <- creationFunction()
				distx <- NULL
				if(length(distanceFunction)>1)
					dfn = distanceFunction[[1]]
				else 
					dfn = distanceFunction
				for(j in 1:length(x)){#calculate distance to each sampled location				
					distx <- c(distx,dfn(xc[[i]],x[[j]]))	#todo apply?
				}
				distlist <- c(distlist,min(distx)) #maximise the minimum distance
			}
			x[[count]] <- xc[[which.max(distlist)]] # todo: this maximises distance, but may still create a duplicate if max(min(dist))==0, e.g. if all randomly created individuals are duplicates.
		}
		x <- removeDuplicates(x, creationFunction)
    if (vectorized)
      y <-  c(y,fun(x[count]))
    else
      y <-  c(y,fun(x[[count]]))
		if(y[[count]]<bestKnownY){ #new best solution
			bestKnownY <- y[[count]]
			bestKnownX <- x[[count]]
			#bestKnownPred <- fpred(bestKnownX)
		}
    if(verbosity > 0)
      print(paste("Evaluations:",count,"    Quality:",bestKnownY))#"    Solution:",paste(x[[count]],collapse="")))
		if(plotting){
			plot(y,type="l")
			abline(bestKnownY,0,lty=2)
		}
		## Update the distance matrix
		#newdist = distanceVector(x[[count]],x[-count],distanceFunction)
		#control$modelSettings$distances = cbind(rbind(control$modelSettings$distances,c(newdist)),c(newdist,0))

		## Update the distance matrix
		if(length(distanceFunction)==1){ # in case of a single distance function (all models)
			newdist = distanceVector(x[[count]],x[-count],distanceFunction)
			control$modelSettings$distances = cbind(rbind(control$modelSettings$distances,c(newdist)),c(newdist,0))
		}else{	# in case of multiple distance functions (kriging only atm.)
			for(i in 1:length(distanceFunction)){
				newdist = distanceVector(x[[count]],x[-count],distanceFunction[[i]])
				control$modelSettings$distances[[i]] <- cbind(rbind(control$modelSettings$distances[[i]],c(newdist)),c(newdist,0))
			}
		}
		##		
	
		##update model and prediction function:
		res <- buildModel(x,y,distanceFunction,control)
		fit <- res$fit
		fpred <- res$fpred		
	}
	list(xbest=bestKnownX, ybest=bestKnownY, x=x,y=y,fit=fit,fpred=fpred,count=count) #return values. do NOT return bestKnownPred, since this is prediction of the best individual, not its actual utility value.
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
#' @param control list of options \cr
#' \code{model} Model to be used as a surrogate of the target function. Default is "K" (Kriging). Also
#'		available are: "LM" (linear, distance-based model), "RBFN" Radial Basis Function Network.\cr
#' \code{modelSettings} List of settings for model building, passed on as the control argument to the model training functions \code{\link{combinatorialKriging}}, \code{\link{combinatorialLM}}, \code{\link{combinatorialRBFN}}.
#' \code{infill} This parameter specifies a function to be used for the infill criterion (e.g., the default is expected improvement \code{infillExpectedImprovement}).
#' To use no specific infill criterion this has to be set to \code{NA}. Infill criteria are only used with models that may provide some error estimate with predictions.\cr
#'
#' @return a list:\cr 	
#' \code{fit} model-fit \cr
#' \code{fpred} prediction function
#' 
#' @seealso \code{\link{optimCEGO}} 
#'
#' @keywords internal
###################################################################################
buildModel <- function(x,y,distanceFunction,control){ 
	y #against lazy evaluation
	if(control$model=="RBFN"){
		fit <- combinatorialRBFN(x,y,distanceFunction,control$modelSettings) 
		fit
		if(is.function(control$infill)){
			fpred <- function(x){	
						res=predict(fit,x,TRUE)
						resEI= control$infill(res$f,res$s,min(y))
					} 
		}else{
			fpred <- function(x){predict(fit,x)$f}
		}
	}else if(control$model=="LM"){
		fit <- combinatorialLM(x,y,distanceFunction,control$modelSettings) 
		fit
		fpred <- function(x){predict(fit,x)} 
	}else if(control$model=="K"){
		fit<-combinatorialKriging(x,as.matrix(y),distanceFunction,control$modelSettings) 
		fit
		if(is.function(control$infill)){
			fpred <- function(x){	
					res=predict(fit,x,TRUE)
					resEI= control$infill(res$f,res$s,min(y))
				} 
		}else{
			fpred <- function(x){predict(fit,x)$f}
		}
	}else if(control$model=="SVR"){
		stop("Using Support Vector Regression (SVR) as a surrogate model (see: control$model option) is not yet implemented in the current version of the CEGO package. Please contact the package maintainer for information on future developments of this feature.")
	}
	list(fit=fit,fpred=fpred)
}

###################################################################################
#' Evolutionary Algorithm for Combinatorial Optimization
#' 
#' A basic implementation of a simple Evolutionary Algorithm for Combinatorial Optimization. Default evolutionary operators
#' aim at permutation optimization problems.
#'
#' @param par Optional start individual(s) as a list. If NULL (default), creationFunction is used to create initial population. 
#' If par has less individuals than the population size, creationFunction will fill up the rest.
#' @param fun target function to be minimized
#' @param control (list), with the options\cr
#' \code{budget} The limit on number of target function evaluations (stopping criterion) (default: 1000)\cr
#' \code{popsize} Population size (default: 100)\cr
#' \code{generations} Number of generations (stopping criterion) (default: Inf)\cr
#' \code{targetY} Target function value (stopping criterion) (default: -Inf)\cr
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
#' \code{verbosity} Level of text output during run. Defaults to 0, no output.\cr
#' \code{plotting} Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.\cr
#' \code{recombinationFunction} Function that performs recombination, default: \code{\link{recombinationPermutationCycleCrossover}}, which is cycle crossover for permutations.\cr
#' \code{recombinationParameters} Parameter list for recombination (e.g., recombinationParameters$recrate => recombination rate, defaults to 0.5). List is passed to recombinationFunction. \cr
#' \code{mutationFunction} Function that performs mutation, default: \code{\link{mutationPermutationSwap}}, which is swap mutation for permutations.\cr
#' \code{mutationParameters} Parameter list for mutation (e.g., mutationParameters$mutationRate => mutation rate, defaults to 0.5).List is passed to mutationFunction. \cr
#' \code{selection} Selection process: "tournament" (default) or "truncation"\cr
#' \code{tournamentSize} Tournament size (default: 2)\cr
#' \code{tournamentProbability} Tournament probability (default: 0.9)\cr
#' \code{localSearchFunction} If specified, this function is used for a local search step. Default is NULL. \cr
#' \code{localSearchRate} Specifies on what fraction of the population local search is applied. Default is zero. Maximum is 1 (100 percent).
#' \code{localSearchSettings} List of settings passed to the local search function control parameter.
#' \code{stoppingCriterionFunction} Custom additional stopping criterion. Function evaluated on the population, receiving all individuals (list) and their fitness (vector). If the result is FALSE, the algorithm stops.
#' \code{verbosity} >0 for text output.
#' @param creationFunction Function to create individuals/solutions in search space
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{x} history of all evaluated solutions\cr
#' \code{y} corresponding target function values f(x)\cr
#' \code{count} number of performed target function evaluations 
#' \code{message} Termination message: Which stopping criterion was reached.
#' \code{population} Last population
#' \code{fitness} Fitness of last population
#' 
#' @examples
#' seed=0
#' glgseed=1
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(ref=1:5,distFun=dF)
#' #start optimization
#' set.seed(seed)
#' res <- optimEA(,lF,cF,list(mutationFunction=mF,recombinationFunction=rF,
#'		popsize=15,budget=100,targetY=0,verbosity=1))
#' res$xbest 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimRS}}, \code{\link{optim2Opt}} 
#' 
#' @export
###################################################################################
optimEA <- function(par=NULL,fun,creationFunction,control=list()){ 
	#default controls:
	con<-list(budget = 100 #default controls:
			 , popsize = 10
			 , generations = Inf
			 , targetY = -Inf
			 , vectorized=FALSE
			 , recombinationParameters = list(recrate=0.5)
			 , mutationParameters = list(mutationRate=0.5)
			 , mutationFunction = mutationPermutationSwap
			 , recombinationFunction = recombinationPermutationCycleCrossover
			 , selection = "tournament" #or "truncation
			 , tournamentSize = 2 
			 , tournamentProbability = 0.9
			 , localSearchFunction = NULL
			 , localSearchRate = 0
			 , localSearchSettings = list()
			 , stoppingCriterionFunction = NULL
			 , verbosity = 0 
       , plotting = FALSE
			 );
	con$recombinationParameters[namc <- names(control$recombinationParameters)] <- control$recombinationParameters
	con$mutationParameters[namc <- names(control$mutationParameters)] <- control$mutationParameters
	control$recombinationParameters <- con$recombinationParameters
	control$mutationParameters <- con$mutationParameters
	con[(namc <- names(control))] <- control;
	control<-con;
	
	budget <- control$budget
	vectorized <- control$vectorized
	popsize <- control$popsize
	generations <- control$generations
	targetY <- control$targetY
	tournamentSize <- control$tournamentSize
	recombinationParameters <- control$recombinationParameters
	recrate <- recombinationParameters$recrate
	recombinationFunction <- control$recombinationFunction
	mutationParameters <-  control$mutationParameters
	mutationRate <- mutationParameters$mutationRate
	mutationFunction <- control$mutationFunction
	localSearchFunction <- control$localSearchFunction
	localSearchRate <- control$localSearchRate
	localSearchSettings <- control$localSearchSettings
	localSearchSettings$vectorized <- vectorized
	selection <- control$selection
	stoppingCriterionFunction <- control$stoppingCriterionFunction
	verbosity <- control$verbosity
	plotting <- control$plotting
	tournamentProbability <-  control$tournamentProbability #probability in the tournament selection
	tournamentSize=max(tournamentSize,1)
	
	## Create initial population
	population <- initializeDesign(par,creationFunction,popsize)
	
	if(vectorized) 
		fitness <- fun(population)
	else
		fitness <- unlist(lapply(population,fun))
		
	count = popsize	
	gen=0
	
	fithist <- fitness
	xhist <- population
	if(plotting){
		besthist <- NULL
	}
	run <- TRUE
	while((count < budget) & (gen < generations) & (min(fithist,na.rm=TRUE) > targetY) & (run)){
		gen=gen+1
		#recombine
		if(selection == "tournament"){ #tournament selection
			popindex <- tournamentSelection(fitness,tournamentSize,tournamentProbability,max(popsize*recrate,2))
		}else{ #truncation selection
			popindex <- order(fitness)[1:max(popsize*recrate,2)]
		}
		offspring <- recombinationFunction(population[popindex],recombinationParameters)
		#mutate
		offspring <- mutationFunction(offspring,mutationParameters)
		#optional local search
		if(!is.null(localSearchFunction) & localSearchRate>0){
			if(localSearchRate<1){
				subsetsize = ceiling(length(offspring)*localSearchRate)
				offspringsubset <- sample(length(offspring),subsetsize)
			}else{
				offspringsubset = 1:length(offspring)
			}	
			for(i in offspringsubset){
				res <- localSearchFunction(par=offspring[i],fun=fun,creationFunction=creationFunction,control=localSearchSettings) 
				offspring[[i]] <- res$xbest #todo: local search already evaluations xbest, is reevaluated in population.
				count <- count + res$count #add local search counted evaluations to evaluation counter of main loop.
			}
		}
		#remove offspring which violate the budget
		offspring <- offspring[1:min(budget-count,length(offspring))]
		#append offspring to population, but remove duplicates first. duplicates are replaced by random, unique solutions.		
		offspring <- removeDuplicatesOffspring(xhist,offspring, creationFunction)
		population <- c(population,  offspring)
		#if any new were created:	
		if(length(population)>length(fitness)){
			xhist=append(xhist,offspring)
			#evaluate
			if(vectorized)
				newfit <- fun(offspring)
			else
				newfit <- unlist(lapply(offspring,fun))
			fithist= c(fithist, newfit)
			fitness <- c(fitness, newfit) #evaluate the new individuals after recombination and mutation
			#update count
			count=count+ length(fitness)-popsize			
			#tournament selection 
			if(selection == "tournament"){ #tournament selection
				popindex <- tournamentSelection(fitness,tournamentSize,tournamentProbability,popsize) #todo should it really be possible to select the same individual several times for the next generation?
			}else{ # truncation selection
				popindex <- order(fitness)[1:popsize]
			}
			population <- population[popindex]
			fitness <- fitness[popindex]			
		}	
		if(!is.null(stoppingCriterionFunction)) # calculation of additional stopping criteria
			run <- stoppingCriterionFunction(population,fitness)
		if(plotting){
			besthist <- c(besthist,min(fithist))
			par(mfrow=c(2,1))
			plot(log(besthist-min(besthist)+1e-16),type="l")
			plot(besthist,type="l")
		}
    if(verbosity > 0){
    	print(paste("Generations: ",gen," Evaluations: ",count, "Best Fitness: ", min(fithist)))
    }
	}
	#stopping criteria information for user:
	msg <- "Termination message:"
	if(!run) #success
		msg=paste(msg,"Custom stopping criterion satisfied.")
	if(min(fithist,na.rm=TRUE) < targetY) 
		msg=paste(msg,"Target fitness achieved.")		
	else if(count >= budget) #budget exceeded
		msg=paste(msg,"Target function evaluation budget depleted.")		
	else if(gen >= generations) #generation limit exceeded
		msg=paste(msg,"Number of generations limit reached.")		

	index <- which.min(fithist)
	list(xbest=xhist[[index]],ybest=fithist[index],x=xhist,y=fithist, count=count, message=msg, population=population, fitness=fitness)#return best individual
}

###################################################################################
#' Combinatorial Random Search
#' 
#' Random Search for mixed or combinatorial optimization.
#'
#' @param par Optional set of solution(s) as a list. If NULL (default), creationFunction is used to create all solutions randomly. 
#' If par has less individuals than the population size, creationFunction will fill up the rest. par basically defines additional user specified solutions to be evaluated.
#' @param fun target function to be minimized
#' @param control (list), with the options\cr
#' \code{budget} The limit on number of target function evaluations (stopping criterion) (default: 100)
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
##' @param creationFunction Function to create individuals/solutions in search space
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{x} history of all evaluated solutions\cr
#' \code{y} corresponding target function values f(x)\cr
#' \code{count} number of performed target function evaluations 
#' 
#' @examples
#' seed=0
#' glgseed=1
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(ref=1:5,distFun=dF)
#' #start optimization
#' set.seed(seed)
#' res <- optimRS(,lF,cF,list(budget=100))
#' res$xbest 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimEA}} 
#' 
#' @export
###################################################################################
optimRS <- function(par=NULL,fun,creationFunction,control=list()){
	con<-list(budget=100,vectorized=FALSE
			 );
	con[(namc <- names(control))] <- control;
	control<-con;

  budget <- control$budget
	vectorized <- control$vectorized
 
	## Create random solutions without duplicates, filled up with par
	x <- initializeDesign(par,creationFunction,budget)
  
  #evaluate
  if(vectorized) 
		y <- fun(x)
	else
		y <- unlist(lapply(x,fun))
  
  #best value found:  
	j <- which.min(y)
  
  #return
	list(xbest=x[[j]],ybest=y[j],x=x,y=y, count=budget)
}



###################################################################################
#' Two-Opt
#' 
#' Implementation of a Two-Opt local search.
#'
#' @param par start solution of the local search
#' @param fun function that determines cost or length of a route/permutation
#' @param control (list), with the options\cr
#' \code{budget} The limit on number of target function evaluations (stopping criterion) (default: 100)
#' \code{vectorized} Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE\cr
#' @param creationFunction Function to create individuals/solutions in search space
#'
#' @return a list:\cr 	
#' \code{xbest} best solution found\cr
#' \code{ybest} fitness of the best solution\cr
#' \code{count} number of performed target function evaluations 
#' 
#' @examples
#' seed=0
#' glgseed=1
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mF <- mutationPermutationSwap
#' #recombination
#' rF <-  recombinationPermutationCycleCrossover 
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(ref=1:5,distFun=dF)
#' #start optimization
#' set.seed(seed)
#' res <- optim2Opt(,lF,cF,list(budget=100))
#' res
#'
#' @references Wikipedia contributors. "2-opt." Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 13 Jun. 2014. Web. 21 Oct. 2014. (http://en.wikipedia.org/wiki/2-opt)
#'
#' @export
###################################################################################
optim2Opt <- function(par=NULL,fun,creationFunction,control=list()){ 
	con<-list(budget=100,vectorized=FALSE
			 );
	con[(namc <- names(control))] <- control;
	control<-con;

  budget <- control$budget
  vectorized <- control$vectorized

  if(is.null(par)){
    route=creationFunction()
  }else{ #given start population
    route=par[[1]]
  }  

	improvement=TRUE
	bestRoute=route
  if (vectorized)
    bestDist = fun(list(route))
  else
    bestDist = fun(route)
	N=length(route)
	count=1
	while(improvement){
		improvement=FALSE
		i=1
		while(i<=(N-1)){	
			for(k in (i+1):N){
				newRoute = step2Opt(bestRoute,i,k)
        if (vectorized)
          newDist = fun(list(newRoute))
        else
          newDist = fun(newRoute)
				count=count+1
				if (newDist < bestDist) {
					bestRoute = newRoute
					bestDist = newDist
					#i=N
					improvement=TRUE
					break;
				}
				if(count == budget){
					improvement=FALSE
					i=N
					break;
				}
			}
			i=i+1	
		}
	}
	#print(count)
	list(xbest=	bestRoute, ybest= bestDist, count=count)
}

###################################################################################
#' 2-Opt Step
#' 
#' Helper function: A single 2-opt step for \code{\link{optim2Opt}}
#'
#' @param route route to be partially reversed
#' @param i start of reversal
#' @param k end of reversal
#'
#' @return a new route
#' 
#' @keywords internal
###################################################################################
step2Opt <- function(route, i, k) { 
	newRoute <- route
	newRoute[i:k] <- rev(route[i:k]) #reversal
	newRoute
}

