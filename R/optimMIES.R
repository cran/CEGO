# R version of a mixed integer evolution strategy by Martin Zaefferer, intended for use in the CEGO package.

###################################################################################
#' Mixed Integer Evolution Strategy (MIES)
#' 
#' An optimization algorithm from the family of Evolution Strategies, designed to
#' optimize mixed-integer problems: The search space is composed of continuous (real-valued) parameters,
#' ordinal integers and categorical parameters.
#' Please note that the categorical parameters need to be coded as integers 
#' (type should not be a factor or character).
#' It is an implementation (with a slight modification) of MIES as described by Li et al. (2013).
#' Note, that this algorithm always has a step size for each solution parameter, unlike Li et al.,
#' we did not include the option to change to a single step-size for all parameters.
#' Dominant recombination is used for solution parameters (the search space parameters), 
#' intermediate recombination for strategy parameters (i.e., step sizes).
#' Mutation: Self-adaptive, step sizes sigma are optimized alongside the solution parameters. 
#' Real-valued parameters are subject to variation based on independent normal distributed random variables. 
#' Ordinal integers are subject to variation based on the difference of geometric distributions. 
#' Categorical parameters are changed at random, with a self-adapted probability.
#' Note, that a more simple bound constraint method is used. Instead of the Transformation Ta,b(x) 
#' described by Li et al., optimMIES simply replaces any value that exceeds the bounds by respective boundary value.
#'
#' The control variables types, lower, upper and levels are especially important.
#'
#' @param x Optional start individual(s) as a list. If NULL (default), \code{creationFunction} (in \code{control} list) is used to create initial design. 
#' If \code{x} has less individuals than the population size, creationFunction will fill up the rest.
#' @param fun target function to be minimized.
#' @param control (list), with the options:
#' \describe{
#'   \item{\code{budget}}{The limit on number of target function evaluations (stopping criterion) (default: 1000).}
#'   \item{\code{popsize}}{Population size (default: 100).}
#'   \item{\code{generations}}{Number of generations (stopping criterion) (default: Inf).}
#'   \item{\code{targetY}}{Target function value (stopping criterion) (default: -Inf).}
#'   \item{\code{vectorized}}{Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE.}
#'   \item{\code{verbosity}}{Level of text output during run. Defaults to 0, no output.}
#'   \item{\code{plotting}}{Plot optimization progress during run (TRUE) or not (FALSE). Default is FALSE.}
#'   \item{\code{archive}}{Whether to keep all candidate solutions and their fitness in an archive (TRUE) or not (FALSE). Default is TRUE.}
#'   \item{\code{stoppingCriterionFunction}}{Custom additional stopping criterion. Function evaluated on the population, receiving all individuals (list) and their fitness (vector). If the result is FALSE, the algorithm stops.}
#'   \item{\code{types}}{A vector that specifies the data type of each variable: "numeric", "integer" or "factor".}
#'   \item{\code{lower}}{Lower bound of each variable. Factor variables can have the lower bound set to NA.}
#'   \item{\code{upper}}{Upper bound of each variable. Factor variables can have the upper bound set to NA.}
#'   \item{\code{levels}}{List of levels for each variable (only relevant for categorical variables). 
#' Should be a vector of numerical values, usually integers, but not necessarily a sequence.
#' HAS to be given if any factors/categoricals are present. Else, set to NA.}
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
#' @references Rui Li, Michael T. M. Emmerich, Jeroen Eggermont, Thomas Baeck, Martin Schuetz, Jouke Dijkstra, and Johan H. C. Reiber. 2013. Mixed integer evolution strategies for parameter optimization. Evol. Comput. 21, 1 (March 2013), 29-64.
#'
#' @examples
#' set.seed(1)
#' controlList <- list(lower=c(-5,-5,1,1,NA,NA),upper=c(10,5,10,10,NA,NA),
#' 		types=c("numeric","numeric","integer","integer","factor","factor"),
#'		levels=list(NA,NA,NA,NA,c(1,3,5),1:4),
#' 		vectorized = FALSE)
#' objFun <- function(x){		
#' 		x[[3]] <- round(x[[3]])
#' 		x[[4]] <- round(x[[4]])
#' 		y <- sum(as.numeric(x[1:4])^2) 
#' 		if(x[[5]]==1 & x[[6]]==4)
#' 			y <- exp(y)
#' 		else
#' 			y <- y^2
#' 		if(x[[5]]==3)
#' 			y<-y-1	
#' 		if(x[[5]]==5)
#' 			y<-y-2	
#' 		if(x[[6]]==1)
#' 			y<-y*2
#' 		if(x[[6]]==2)
#' 			y<-y * 1.54
#' 		if(x[[6]]==3)
#' 			y<- y +2
#' 		if(x[[6]]==4)
#' 			y<- y * 0.5	
#' 		if(x[[5]]==1)
#' 			y<- y * 9	
#' 		y	
#' 	}
#' res <- optimMIES(,objFun,controlList)
#' res$xbest
#' res$ybest
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimRS}}, \code{\link{optimEA}}, \code{\link{optim2Opt}}, \code{\link{optimMaxMinDist}} 
#' 
#' @export
###################################################################################
optimMIES <- function(x=NULL,fun,control=list()){ #TODO: providing x seems to be buggy
	#default controls: #todo: document
	con<-list(budget = 1000 #default controls:
			 , popsize = 100 # mu
			 , lambda = 2 #
			 , generations = Inf
			 , targetY = -Inf
			 , vectorized=FALSE
			 , isotropic = FALSE #one step size for all parameters, or separate for each dimension?
			 , strategy = "plus" # mu + lambda strategy ("plus") or mu,lambda strategy ("comma")
#			 , lower = 0, #lower bounds for all parameters. HAS to be given.
#			 , upper = 1, #upper bounds for all parameters. HAS to be given.
#			 , types = NA, # vector with data types for each dimension. HAS to be given.
#			 , levels = NA, # list of levels for each variable (only relevant for categorical variables). HAS to be given if any factors/categoricals are present.
       , archive = FALSE
			 , stoppingCriterionFunction = NULL
			 , verbosity = 0 
			 , plotting = FALSE 
			 );
	con[names(control)] <- control;
	control<-con;

  archive <- control$archive
	budget <- control$budget
	vectorized <- control$vectorized
	popsize <- control$popsize
	generations <- control$generations
	targetY <- control$targetY
	stoppingCriterionFunction <- control$stoppingCriterionFunction
	verbosity <- control$verbosity
	plotting <- control$plotting
	
	lambda <- control$lambda
	lower <- control$lower
	upper <- control$upper
	types <- control$types
	levels <- control$levels
	isotropic <- control$isotropic
	sigma0 <- control$sigma0
	strategy <- control$strategy
	isotropic
	types
	levels
	lower
	upper
	
	if(is.null(types))
			stop("Please provide the vector types, which gives the data type (numeric, integer, factor) for each variable.")
		
	
	if(is.null(lower)|is.null(upper))
		stop("Please provide lower and upper bounds for the search space in optimMIES, via control list.")
	
	npar <- length(lower)  #number of parameters
	nreal <- sum(types=="numeric")  #number of real parameters
	nint <-  sum(types=="integer") #number of integer parameters
	ncat <-  sum(types=="factor")  #number of categorical parameters
	ireal <- which(types=="numeric")
	iint <- which(types=="integer")
	icat <- which(types=="factor")
	tauReal <- 1/sqrt(2*nreal) #global learning rate
	tauRealDash <- 1/sqrt(2*sqrt(nreal)) #local learning rate
	tauInt <- 1/sqrt(2*nint) #global learning rate
	tauIntDash <- 1/sqrt(2*sqrt(nint)) #local learning rate
	tauCat <- 1/sqrt(2*ncat) #global learning rate
	tauCatDash <- 1/sqrt(2*sqrt(ncat)) #local learning rate

	ranges <- upper-lower

	if(is.null(sigma0)){
		sigma0 <- numeric(npar)
		sigma0[ireal] <- ranges[ireal] * 0.1
		sigma0[iint] <- ranges[iint] * 0.33
		sigma0[icat] <- 0.1
	}
	
	if(npar != length(upper))
		stop("Lower and upper bounds for the search space in optimMIES have to be of the same length. Check the control list.")

		
		
	if(length(sigma0)!=npar)
			stop("Length of initial step size vector (sigma0) should be one, or same length as lower/upper bounds.")

	
	creationFunction <- function(){
	  x <- numeric(npar)
		if(length(ireal)>0)
			x[ireal] <- runif(nreal,lower[ireal],upper[ireal])
		if(length(iint)>0)
			x[iint] <- floor(runif(nint,lower[iint],upper[iint]+ 1-.Machine$double.eps ))
		if(length(icat)>0)
			x[icat] <- sapply(levels[icat], sample,1,simplify=T)
		## append strategy parameters
		x <- c(x,sigma0)
		x
	}
	#creationFunction()
	
	recombinationFunction <- function(parent1,parent2){ 
		## dominant recombination for solution parameters
		inds <- sample(1:npar,ceiling(npar/2))
		parent1[inds] <- parent2[inds]
		## intermediate recombination for strategy parameters
		parent1[npar+1] <- (parent1[npar+1] + parent2[npar+1]) / 2
		parent1
	}
	
	#Tlu <- function(x,a,b){ #todo!
	#	y <- (x-a) / (b-a)
	#	ydash <- y
	#	inds <- (floor(y) %% 2) == 0
	#	ydash[inds] <-  abs(y[inds]-floor(y[inds]))
	#	ydash[!inds] <-  1 - abs(y[!inds]-floor(y[!inds]))
	#	xdash <- a + (b - a) * ydash
	#	xdash
	#}
	#Tlu(1:4,c(2,2,2,2),c(3,3,3,3))
	#Tlu(-1:4,c(1,1,1,1,-1,-1),c(2,3,4,2,10,10))
	
	sigids <- (npar+1):(npar*2)
	
	#TODO: auslagern
	mutationFunctionReal <- function(individual){ 
		Nc <- rnorm(1,0,1)
		sigma <- individual[sigids][ireal]
		sigmaDash <- sigma * exp(tauReal * Nc + tauRealDash * rnorm(nreal,0,1))
		newval <-  as.numeric(individual[ireal]) + sigmaDash * rnorm(nreal,0,1)
		#individual[ireal] <- Tlu(newval,lower,upper)
		newval <- pmin(upper[ireal],newval) #fix solution parameters to bounds #TODO: Tab(x) implementation
		individual[ireal] <- pmax(lower[ireal],newval) 				
		individual[sigids][ireal] <- sigmaDash
		individual
	}	

	
	mutationFunctionInt <- function(individual){ 
		Nc <- rnorm(1,0,1)
		sigma <- individual[sigids][iint]
		sigmaDash <- pmax(1,sigma * exp(tauInt * Nc + tauIntDash * rnorm(nint,0,1)))
		#u1 <- runif(nint)
		#u2 <- runif(nint)
		p <- sigmaDash/nint
		p <- 1-p/(1+sqrt(1+p^2))
		G1 <- rgeom(nint,prob=p) #TODO: may be NA for very small p (~1e-10)
		G2 <- rgeom(nint,prob=p)
		newval <- as.numeric(individual[iint]) + G1-G2
		#individual[iint] <- Tlu(newval,lower[iint],upper[iint]) #fix solution parameters to bounds with transformation
		newval <- pmin(upper[iint],newval) #fix solution parameters to bounds
		individual[iint] <- pmax(lower[iint],newval) 
		individual[sigids][iint] <- sigmaDash
		individual
	}	

	
	mutationFunctionCat <- function(individual){  
		Nc <- rnorm(1,0,1)
		sigma <- individual[sigids][icat]
		sigmaDash <- 1/(1+((1-sigma)/sigma * exp(-tauCat * Nc - tauCatDash * rnorm(ncat,0,1))))
		#sigmaDash <- Tlu(sigmaDash,rep(1/(3*ncat),ncat),rep(0.5,ncat))
		sigmaDash <-  pmin(0.5,sigmaDash)
		sigmaDash <-  pmax(1/(3*ncat),sigmaDash)
		u <- runif(ncat)
		inds <- u < sigmaDash
		newval <- sapply(levels[icat], sample,1)
		individual[icat][inds] <- newval[inds]
		individual[sigids][icat] <- sigmaDash
		individual
	}
	#x1 <- creationFunction()
	#x1
	#mutationFunctionCat(x1)	
	
	## Create initial population
	population <- designRandom(x,creationFunction,popsize)
	
	if(vectorized) 
		fitness <- fun(sapply(population,'[',-sigids,simplify=F)) #note: this also first cuts off strategy parameters
	else
		fitness <- unlist(lapply(sapply(population,'[',-sigids,simplify=F),fun))#note: this also first cuts off strategy parameters
		
	count <- popsize	
	gen <- 1
	
  fitbest <- min(fitness,na.rm=TRUE)
  xbest <- population[[which.min(fitness)]][-sigids]
  if(archive){
    fithist <- fitness
    xhist <- population
  }
	besthist <- fitbest # initialization for plotting
	run <- TRUE
	while((count < budget) & (gen < generations) & (fitbest > targetY) & (run)){
		gen <- gen+1
		#recombine
		c1 <- sample(popsize,lambda,replace=TRUE) #first parents
		c2 <- sample(popsize,lambda,replace=TRUE) #second parents
		offspring <- mapply(FUN=recombinationFunction,population[c1],population[c2],SIMPLIFY=FALSE)
		#mutate real, integer, factors:		
		offspring <- sapply(offspring,mutationFunctionReal,simplify=FALSE)#todo: if any
		offspring <- sapply(offspring,mutationFunctionInt,simplify=FALSE)#todo: if any
		offspring <- sapply(offspring,mutationFunctionCat,simplify=FALSE)#todo: if any
		
		
		if(length(offspring)>0 & budget > count){
			## remove offspring which violate the budget 
			offspring <- offspring[1:min(budget-count,length(offspring))]
			#evaluate
			if(vectorized)
				newfit <- fun(sapply(offspring,'[',-sigids,simplify=F))#note: this also first cuts off strategy parameters
			else
				newfit <- unlist(lapply(sapply(offspring,'[',-sigids,simplify=F),fun))#note: this also first cuts off strategy parameters
			####newfit <- fun(unname(split(offspring[,1:npar],1:nrow(offspring))))#note: this also first cuts off strategy parameters
			
			
			#update count
			count=count+ length(newfit)    
			# keep archive
			if(archive){
        xhist <- append(xhist,offspring)  
        fithist <-  c(fithist, newfit)
      } 
			# remember best
			newbest <- min(newfit,na.rm=TRUE)
			if(newbest < fitbest){
				fitbest <- newbest
				xbest <- offspring[[which.min(newfit)]][-sigids]
			}  				
			if(strategy == "plus"){
				population <- c(population,  offspring)
				fitness <- c(fitness, newfit) 
			}else if(strategy == "comma"){
				population <- offspring
				fitness <- newfit 
			}else{
				stop("Invalid strategy string for MIES, please use plus or comma.")		
			}
		}
		if(length(population)>popsize){		
			#tournament selection 
			#if(selection == "tournament"){ #tournament selection
			#	popindex <- tournamentSelection(fitness,tournamentSize,tournamentProbability,popsize)
			#}else{ # truncation selection
				popindex <- order(fitness)[1:popsize] #MIES seems to use truncation selection only.
			#}
			population <- population[popindex]
			fitness <- fitness[popindex]			
		}	
		if(!is.null(stoppingCriterionFunction)) # calculation of additional stopping criteria
			run <- stoppingCriterionFunction(population,fitness)
		if(plotting){
      besthist <- c(besthist,fitbest)
			plot(besthist,type="l")
		}
    if(verbosity > 0){
    	print(paste("Generations: ",gen," Evaluations: ",count, "Best Fitness: ", min(fitness,na.rm=TRUE)))
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
    return(list(xbest=xbest,ybest=fitbest,x=xhist,y=fithist, count=count, message=msg, population=population, fitness=fitness))
  else
    return(list(xbest=xbest,ybest=fitbest,count=count, message=msg, population=population, fitness=fitness)) 
}
