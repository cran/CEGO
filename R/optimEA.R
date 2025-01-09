
###################################################################################
#' Evolutionary Algorithm for Combinatorial Optimization
#' 
#' A basic implementation of a simple Evolutionary Algorithm for Combinatorial Optimization. Default evolutionary operators
#' aim at permutation optimization problems.
#'
#' @param x Optional start individual(s) as a list. If NULL (default), \code{creationFunction} (in \code{control} list) is used to create initial design. 
#' If \code{x} has less individuals than the population size, creationFunction will fill up the rest.
#' @param fun target function to be minimized
#' @param control (list), with the options:
#' \describe{
#'   \item{\code{budget}}{The limit on number of target function evaluations (stopping criterion) (default: 1000).}
#'   \item{\code{popsize}}{Population size (default: 100).}
#'   \item{\code{generations}}{Number of generations (stopping criterion) (default: Inf).}
#'   \item{\code{targetY}}{Target function value (stopping criterion) (default: -Inf).}
#'   \item{\code{vectorized}}{Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE.}
#'   \item{\code{verbosity}}{Level of text output during run. Defaults to 0, no output.}
#'   \item{\code{plotting}}{Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.}
#'   \item{\code{archive}}{Whether to keep all candidate solutions and their fitness in an archive (TRUE) or not (FALSE). Default is TRUE. New solutions that are identical to an archived one, will not be evaluated. Instead, their fitness is taken from the archive.}
#'   \item{\code{recombinationFunction}}{Function that performs recombination, default: \code{\link{recombinationPermutationCycleCrossover}}, which is cycle crossover for permutations.}
#'   \item{\code{recombinationRate}}{Number of offspring, defined by the fraction of the population (popsize) that will be recombined.}
#'   \item{\code{mutationFunction}}{Function that performs mutation, default: \code{\link{mutationPermutationSwap}}, which is swap mutation for permutations.}
#'   \item{\code{parameters}}{Default parameter list for the algorithm, e.g., mutation rate, etc.}
#'   \item{\code{selection}}{Survival selection process: "tournament" (default) or "truncation".}
#'   \item{\code{tournamentSize}}{Tournament size (default: 2).}
#'   \item{\code{tournamentProbability}}{Tournament probability (default: 0.9).}
#'   \item{\code{localSearchFunction}}{If specified, this function is used for a local search step. Default is NULL. }
#'   \item{\code{localSearchRate}}{Specifies on what fraction of the population local search is applied. Default is zero. Maximum is 1 (100 percent).}
#'   \item{\code{localSearchSettings}}{List of settings passed to the local search function control parameter.}
#'   \item{\code{stoppingCriterionFunction}}{Custom additional stopping criterion. Function evaluated on the population, receiving all individuals (list) and their fitness (vector). If the result is FALSE, the algorithm stops.}
#'   \item{\code{verbosity}}{>0 for text output.}
#'   \item{\code{creationFunction}}{Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6.}
#'   \item{\code{selfAdaption}}{An optional list object, that describes parameters of the optimization (see \code{parameters}) which are subject to self-adaption. An example is given in \link{mutationSelfAdapt}. Types of the parameters can be integer, discrete, or numeric.}
#'   \item{\code{selfAdaptTau}}{Positive numeric value, that controls the learning rate of numerical/integer self-adaptive parameters.}
#'   \item{\code{selfAdaptP}}{Value in [0,1]. A probability of mutation for all categorical, self-adaptive parameters.}
#' }
#'
#' @return a list:
#' \describe{
#'   \item{\code{xbest}}{best solution found.}
#'   \item{\code{ybest}}{fitness of the best solution.}
#'   \item{\code{x}}{history of all evaluated solutions.}
#'   \item{\code{y}}{corresponding target function values f(x).}
#'   \item{\code{count}}{number of performed target function evaluations.}
#'   \item{\code{message}}{Termination message: Which stopping criterion was reached.}
#'   \item{\code{population}}{Last population.}
#'   \item{\code{fitness}}{Fitness of last population.}
#' }
#' 
#' @examples
#' #First example: permutation optimization
#' seed=0
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
#' res <- optimEA(,lF,list(creationFunction=cF,mutationFunction=mF,recombinationFunction=rF,
#'		popsize=6,budget=60,targetY=0,verbosity=1,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res$xbest 
#' #Second example: binary string optimization
#' #number of bits
#' N <- 50
#' #target function (simple example)
#' f <- function(x){
#'	 sum(x)
#' }
#' #function to create random Individuals
#' cf <- function(){
#' 		sample(c(FALSE,TRUE),N,replace=TRUE)
#' }
#' #control list
#' cntrl <- list(
#' 	budget = 100,
#' 	popsize = 5,
#' 	creationFunction = cf,
#' 	vectorized = FALSE, #set to TRUE if f evaluates a list of individuals
#' 	recombinationFunction = recombinationBinary2Point,
#' 	recombinationRate = 0.1,
#' 	mutationFunction = mutationBinaryBitFlip,
#' 	parameters=list(mutationRate = 1/N),
#' 	archive=FALSE #recommended for larger budgets. do not change.
#' )
#' #start algorithm
#' set.seed(1)
#' res <- optimEA(fun=f,control=cntrl)
#' res$xbest
#' res$ybest
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimRS}}, \code{\link{optim2Opt}}, \code{\link{optimMaxMinDist}} 
#' 
#' @export
#'
#' @import fastmatch
###################################################################################
optimEA <- function(x=NULL,fun,control=list()){ 
	#default controls:
	con<-list(budget = 1000 #default controls:
			 , popsize = 100
			 , generations = Inf
			 , targetY = -Inf
			 , vectorized=FALSE
			 , creationFunction = solutionFunctionGeneratorPermutation(6)
			 , recombinationRate=0.5
			 , parameters = list() #manually set paramters for recombination and mutation TODO
			 , mutationFunction = mutationPermutationSwap
			 , recombinationFunction = recombinationPermutationCycleCrossover
			 , selection = "tournament" #or "truncation
			 , tournamentSize = 2 
			 , tournamentProbability = 0.9
			 , localSearchFunction = NULL
			 , localSearchRate = 0
			 , localSearchSettings = list()
       , archive = TRUE
			 , stoppingCriterionFunction = NULL
			 , verbosity = 0 
       , plotting = FALSE		
			 , selfAdaptTau = 1/sqrt(2)
			 , selfAdaptP = 0.5
			 );
	con[names(control)] <- control;
	control<-con;

  archive <- control$archive
	creationFunction <- control$creationFunction	
	budget <- control$budget
	vectorized <- control$vectorized
	popsize <- control$popsize
	generations <- control$generations
	targetY <- control$targetY
	tournamentSize <- control$tournamentSize
	parameters <- control$parameters
	recombinationRate <- control$recombinationRate
	recombinationFunction <- control$recombinationFunction
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
	selfAdaption <-  control$selfAdaption 
	selfAdaptTau <-  control$selfAdaptTau 
	selfAdaptP <-  control$selfAdaptP 
	tournamentSize <- max(tournamentSize,1)
	
	#for backwards compatibility only:
	if(!is.null( control$mutationParameters))
		parameters[names(control$mutationParameters)] <- control$mutationParameters
	if(!is.null( control$recombinationParameters))
		parameters[names(control$recombinationParameters)] <- control$recombinationParameters
	
	## Create initial population
	population <- designRandom(x,creationFunction,popsize)
	
	## Create initial parameters (if self-adaptive)
	populationSelfAdapt <- NA
	if(!is.null(selfAdaption)){
		populationSelfAdapt <- t(replicate(popsize,lapply(lapply(selfAdaption,'[[',"default"),eval))) 
		#t(replicate(popsize,getDefaults(old)))
		types <- sapply(selfAdaption,'[[',"type") #getParamTypes(selfAdaption)
		lower <- unlist(sapply(selfAdaption,'[[',"lower")) #getLower(selfAdaption)
		upper <- unlist(sapply(selfAdaption,'[[',"upper")) #getUpper(selfAdaption)
		values <- lapply(selfAdaption,'[[',"values") 
		valnotnull <- !sapply(values,is.null)
		values <- values[valnotnull]
		  #values <- sapply(getValues(selfAdaption,dict=list()),unlist,simplify=FALSE)
		inum <- types=="numeric"
		icat <- types=="discrete"
		iint <- types=="integer"
		nnum <- sum(inum)
		ncat <- sum(icat)
		nint <- sum(iint)			
	}	
	
	if(vectorized) 
		fitness <- fun(population)
	else
		fitness <- unlist(lapply(population,fun))
		
	count <- popsize	
	gen <- 1
	
  fitbest <- min(fitness,na.rm=TRUE)
  xbest <- population[[which.min(fitness)]]
  if(archive){
    fithist <- fitness
    xhist <- population
  }
	besthist <- fitbest # initialization for plotting
	run <- TRUE
	while((count < budget) & (gen < generations) & (fitbest > targetY) & (run)){
		gen <- gen+1
		## parent selection
		if(selection == "tournament"){ #tournament selection
			parents <- tournamentSelection(fitness,tournamentSize,tournamentProbability,max(floor(popsize*recombinationRate)*2,2))
		}else{ #truncation selection
			parents <- rep(order(fitness),2)[1:max(floor(popsize*recombinationRate)*2,2)]
		}
		## shuffle parents
		parents <- sample(parents)
		##  self-adapt parameters (recombine, apply learning)
		if(!is.null(selfAdaption)){
			offspringSelfAdapt <- selfAdapt(populationSelfAdapt[parents,],inum,icat,iint,nnum,ncat,nint,lower,upper,values,selfAdaptTau,selfAdaptP) #adapt parameters (learn)
			parameters$selfAdapt <- offspringSelfAdapt #get parameters of the offspring individuals
		}				
		## recombine parents
		offspring <- recombinationFunction(population[parents],parameters)		
		## mutate parents
		offspring <- mutationFunction(offspring,parameters)
		#optional local search
		if(!is.null(localSearchFunction) & localSearchRate>0){
			if(localSearchRate<1){
				subsetsize <- ceiling(length(offspring)*localSearchRate)
				offspringsubset <- sample(length(offspring),subsetsize)
			}else{
				offspringsubset <- 1:length(offspring)
			}	
			evaluated <- rep(FALSE,length(offspring))
			tempfitness <- NULL 
			for(i in offspringsubset){
				if(localSearchSettings$budget > (budget-count)) #local search budget exceeds remaining budget
					localSearchSettings$budget <- budget-count #set to remaining budget
				if(localSearchSettings$budget < 2) #local search budget too small
					break
				res <- localSearchFunction(x=offspring[i],fun=fun,control=localSearchSettings) 
				if(archive){ #archive local search results
					xhist <- append(xhist,res$x)  
					fithist <-  c(fithist, res$y)
				} 
				offspring[[i]] <- res$xbest 
				evaluated[i] <- TRUE
				tempfitness <- c(tempfitness,res$ybest)
				count <- count + res$count #add local search counted evaluations to evaluation counter of main loop.
			}
			if(any(evaluated)){ #add only the evaluated individuals to the population, the rest is added later. this is for efficiency reasons only.
				offspring <- offspring[-evaluated]
				offspringLocal <- offspring[evaluated]
				population <- c(population,  offspringLocal)
				if(!is.null(selfAdaption)){
					populationSelfAdapt <- rbind(populationSelfAdapt,offspringSelfAdapt[evaluated,])
					offspringSelfAdapt <- offspringSelfAdapt[-evaluated,]
				}
				# remember best
				newbest <- min(tempfitness,na.rm=TRUE)
				if(newbest < fitbest){
					fitbest <- newbest
					xbest <- offspringLocal[[which.min(tempfitness)]]
				}  			
				fitness <- c(fitness, tempfitness)
			}				
		}
		if(length(offspring)>0 & budget > count){
			## append offspring to population, but remove duplicates first. duplicates are replaced by random, unique solutions.		
			offspring <- removeDuplicatesOffspring(population,offspring, creationFunction,duplicated)		
			newfit <- rep(NA,length(offspring)) # the vector of new fitness values
			evaluated <- rep(FALSE,length(offspring)) # the vector of indicators, whether the respective offspring is already evaluated
			## if archive available, take fitness from archive
			if(archive){ 
				inArchiveIDs <- fmatch(offspring,xhist) #fastmatch package
				inArchive <- !is.na(inArchiveIDs)
				if(any(inArchive)){					
					newfit[inArchive] <- fithist[inArchiveIDs[inArchive]]					
					evaluated[inArchive] <- TRUE
				}
			}			
			## evaluate the rest with fun
			if(any(!evaluated)){
				## remove offspring which violate the budget 
				evaluateOffspring <- offspring[!evaluated]
				possibleEvaluations <- min(budget-count,length(evaluateOffspring)) #number of possible evaluations, given the budget
				evaluateOffspring <- evaluateOffspring[1:possibleEvaluations] 
				#the non-evaluated ones receive NAs.
				
				## evaluate
				if(vectorized)
					newfit[!evaluated][1:possibleEvaluations] <- fun(evaluateOffspring)
				else
					newfit[!evaluated][1:possibleEvaluations] <- unlist(lapply(evaluateOffspring,fun))		
			
				## save to archive
				if(archive){
					xhist <- append(xhist,evaluateOffspring)  
					fithist <-  c(fithist, newfit[!evaluated][1:possibleEvaluations])		
				}
				
				## mark evaluated			
				evaluated[!evaluated][1:possibleEvaluations] <- TRUE	
				
				#update count
				count <- count+ length(evaluateOffspring)
			}								
			## append results to population
			population <- c(population,  offspring[evaluated]) 		
			if(!is.null(selfAdaption)){ 
				populationSelfAdapt <- rbind(populationSelfAdapt,offspringSelfAdapt[evaluated,])
			}
			fitness <- c(fitness, newfit[evaluated]) 

			## remember best
			newbest <- min(newfit,na.rm=TRUE)
			if(newbest < fitbest){
				fitbest <- newbest
				xbest <- offspring[[which.min(newfit)]]
			}  			
		}
		if(length(population)>popsize){		
			#tournament selection 
			if(selection == "tournament"){ #tournament selection
				popindex <- tournamentSelection(fitness,tournamentSize,tournamentProbability,popsize)
			}else{ # truncation selection
				popindex <- order(fitness)[1:popsize]
			}
			population <- population[popindex]
			fitness <- fitness[popindex]
			if(!is.null(selfAdaption)){			
				populationSelfAdapt <- populationSelfAdapt[popindex,]
			}
		}	
		if(!is.null(stoppingCriterionFunction)) # calculation of additional stopping criteria
			run <- stoppingCriterionFunction(population,fitness)
		if(plotting){
      besthist <- c(besthist,fitbest)
			plot(besthist,type="l")
		}
    if(verbosity > 0){
    	print(paste("Generations: ",gen," Evaluations: ",count, "Best Fitness: ",fitbest))
    }
	}
	#stopping criteria information for user:
	msg <- "Termination message:"
	if(!run) #success
		msg=paste(msg,"Custom stopping criterion satisfied.")
	if(fitbest <= targetY) 
		msg=paste(msg,"Successfully achieved target fitness.")		
	else if(count >= budget) #budget exceeded
		msg=paste(msg,"Target function evaluation budget depleted.")		
	else if(gen >= generations) #generation limit exceeded
		msg=paste(msg,"Number of generations limit reached.")		
	
  if(archive)
    return(list(xbest=xbest,ybest=fitbest,x=xhist,y=fithist, count=count, message=msg, population=population, fitness=fitness, populationSelfAdapt= populationSelfAdapt))
  else
    return(list(xbest=xbest,ybest=fitbest,count=count, message=msg, population=population, fitness=fitness , populationSelfAdapt= populationSelfAdapt)) 
}


###################################################################################
#' Self-adaption of EA parameters.
#'
#' Learning / self-adaption of parameters of the evolutionary algorithm.
#'
#' @param params parameters to be self-adapted
#' @param inum boolean vector, which parameters are numeric
#' @param icat boolean vector, which parameters are discrete, factors
#' @param iint boolean vector, which parameters are integer
#' @param nnum number of numerical parameters
#' @param ncat number of discrete parameters
#' @param nint number of integer parameters
#' @param lower lower bounds (numeric, integer parameters only)
#' @param upper upper bounds (numeric, integer parameters only)
#' @param values values or levels of the discrete parameters
#'
#' @seealso \code{\link{optimEA}}
#' 
#' @export
#' @keywords internal
###################################################################################
selfAdapt <- function(params,inum,icat,iint,nnum,ncat,nint,lower,upper,values,tau,p){
	noff <- nrow(params)/2 #number of offspring, assumes 2 parents for each offspring
	#tau <- 1/sqrt(2) # global parameter for this. mutation rate of the numeric/integer paremters
	#p <- 0.5# global parameter for this. mutationr ate of the categorical parameters.
	ainum <- any(inum)
	aiint <- any(iint)
	aicat <- any(icat)

	## first: recombine
	if(ainum) ## numerical: intermediate xover
		params[1:noff,inum] <-  (as.numeric(params[1:noff,inum]) + as.numeric(params[(1:noff)+noff,inum])) / 2
	if(aiint)	## integer: round, intermediate
		params[1:noff,iint] <-  round((as.numeric(params[1:noff,iint]) + as.numeric(params[(1:noff)+noff,iint])) / 2)
	if(aicat){ ## discrete: dominant xover
		index <- sample(1:(noff*ncat),noff*ncat / 2) #select 50 % of the discrete individuals and parameters, to be taken from parent 2
		params[1:noff,icat][index] <-  params[(1:noff)+noff,icat][index]
	}
	params <- params[1:noff,]	

	## second: mutate the parameters
	if(ainum)
		params[,inum] <- as.numeric(params[,inum]) * matrix(exp(tau*rnorm(nnum*noff,0,1)),noff,nnum)
	if(aiint)	
		params[,iint] <- round(as.numeric(params[,inum]) * matrix(exp(tau*rnorm(nint*noff,0,1)),noff,nint))
	if(aicat){
		rand <- matrix(runif(ncat*noff),noff,ncat) #get random number
		ind <- rand < p # if larger than 0.5, change strategy parameter
		if(any(ind)){
			params[,icat][ind] <-sapply(values,FUN=sample,size=noff,replace=TRUE)[ind]
		}
	}
	
	## repair: fix to lower/upper
	if(ainum|aiint){
		params[,inum|iint] <-  pmax(lower,params[,inum|iint]) 
		params[,inum|iint] <-  pmin(upper,params[,inum|iint]) 
	}
	params
}

###################################################################################
#' Self-adaptive mutation operator
#'
#' This mutation function selects an operator and mutationRate (provided in parameters$mutationFunctions)
#' based on self-adaptive parameters chosen for each individual separately.
#'
#' @param population List of permutations
#' @param parameters list, contains the available single mutation functions (\code{mutationFunctions}), 
#' and a data.frame that collects the chosen function and mutation rate for each individual (\code{selfAdapt}).
#'
#' @seealso \code{\link{optimEA}}, \code{\link{recombinationSelfAdapt}}
#' 
#' @export
#'
#' @examples
#' seed=0
#' N=5
#' #distance
#' dF <- distancePermutationHamming
#' #mutation
#' mFs <- c(mutationPermutationSwap,mutationPermutationInterchange,
#'	mutationPermutationInsert,mutationPermutationReversal)
#' rFs <- c(recombinationPermutationCycleCrossover,recombinationPermutationOrderCrossover1,
#'	recombinationPermutationPositionBased,recombinationPermutationAlternatingPosition)
#' mF <- mutationSelfAdapt
#' selfAdaptiveParameters <- list(
#' 	mutationRate = list(type="numeric", lower=1/N,upper=1, default=1/N),
#' 	mutationOperator = list(type="discrete", values=1:4, default=expression(sample(4,1))), 
#'	#1: swap, 2: interchange, 3: insert, 4: reversal mutation
#'	recombinationOperator = list(type="discrete", values=1:4, default=expression(sample(4,1))) 
#'	#1: CycleX, 2: OrderX, 3: PositionX, 4: AlternatingPosition
#')
#' #recombination
#' rF <-  recombinationSelfAdapt	 
#' #creation
#' cF <- function()sample(N)
#' #objective function
#' lF <- landscapeGeneratorUNI(1:N,dF)
#' #start optimization
#' set.seed(seed)
#' res <- optimEA(,lF,list(parameters=list(mutationFunctions=mFs,recombinationFunctions=rFs),
#'		creationFunction=cF,mutationFunction=mF,recombinationFunction=rF,
#'		popsize=15,budget=100,targetY=0,verbosity=1,selfAdaption=selfAdaptiveParameters,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res$xbest 
###################################################################################
mutationSelfAdapt <- function(population,parameters){
	mutfuns <- parameters$mutationFunctions
	mutfunsi <- mutfuns[as.numeric(parameters$selfAdapt[,"mutationOperator"])]
	pars <- parameters$selfAdapt[,"mutationRate"]
	population <- sapply(population,list,simplify=FALSE) #make each parameter to list, because required by operators. #TODO: more efficient to have separate functions that are able to deal with single non-list input
	population <- mapply(function(f, x, p) unlist(f(x,list(mutationRate=p))), mutfunsi, population,pars,SIMPLIFY=FALSE)
	population
}

###################################################################################
#' Self-adaptive recombination operator
#'
#' This recombination function selects an operator (provided in parameters$recombinationFunctions)
#' based on self-adaptive parameters chosen for each individual separately.
#'
#' @param population List of permutations
#' @param parameters list, contains the available single mutation functions (\code{mutationFunctions}), 
#' and a data.frame that collects the chosen function and mutation rate for each individual (\code{selfAdapt}).
#'
#' @seealso \code{\link{optimEA}}, \code{\link{mutationSelfAdapt}}
#' 
#' @export
###################################################################################
recombinationSelfAdapt <- function(population,parameters){
	recfuns <- parameters$recombinationFunctions
	recfunsi <- recfuns[as.numeric(parameters$selfAdapt[,"recombinationOperator"])]
	n <- length(population)
	population <- sapply(population,list,simplify=FALSE) #make each parameter to list, because required by operators. #TODO: more efficient to have separate functions that are able to deal with single non-list input
	population <- mapply(function(f, x1, x2) unlist(f(append(x1,x2))), recfunsi, population[1:(n/2)],population[(1+n/2):n],SIMPLIFY=FALSE)
	population
}
