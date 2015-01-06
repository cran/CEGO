#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

##############################################################
#' Negative Logarithm of Expected Improvement
#'
#' This function calculates the Expected Improvement" of candidate solutions,
#' based on predicted means, standard deviations (uncertainty) and
#' the best known objective function value so far.
#'
#' @param mean predicted mean values 
#' @param sd predicted standard deviation
#' @param min minimum of all observations so far
#'
#' @return Returns the negative logarithm of the Expected Improvement.
#' @export
##############################################################
infillExpectedImprovement <- function(mean,sd,min){
	im <- min - mean
	imsd <- im/sd
	termOne <- im*pnorm(imsd)
	termTwo <- sd*dnorm(imsd) # == sd/sqrt(2*pi)*exp(-0.5*(imsd)^2)
	ei <- -log10(termOne+termTwo+(.Machine$double.xmin)) 
	ei[is.na(ei)] <- Inf
	ei
}

###################################################################################
#' Initialize Random Design
#' 
#' Create a random initial population or experimental design, given a specifed creation function,
#' as well as a optional set of user-specified design members and a maximum design size.
#' Also removes duplicates from the design/population.
#'
#' @param par Optional list of user specified solutions to be added to the design/population, defaults to NULL
#' @param cf Creation function, creates random new individuals
#' @param size size of the design/population
#'
#' @return Returns list with experimental design (or population) without duplicates
#'
#' @keywords internal
#' @export
###################################################################################
initializeDesign <- function(par=NULL,cF,size){
	## initialization
	if(is.null(par)){
    x <- list()
    k=0
  }else{ #given start population
		x <- par
    k=length(x)
  }
		
	if(k>size){
		x <- x[1:size]
	}else if(k<size){
		## CREATE initial population
		x <- c(x, replicate(size-k , cF(),simplify=FALSE))
	}#else if k==size do nothing.
		
	## REPLACE duplicates from initial population with unique individuals
	x <- removeDuplicates(x, cF)
}

###################################################################################
#' Remove Duplicates
#' 
#' Remove duplicates in \code{x}, replace with non-duplicated individuals according to \code{cf}.
#'
#' @param x List of individuals
#' @param cf Creation function, creates random new individuals
#'
#' @return Returns \code{x} without duplicates
#'
#' @keywords internal
#' @export
###################################################################################
removeDuplicates <- function(x,cf){
	while(any(duplicated(x))){ 
		duplicates <- which(duplicated(x))
		for(i in 1:length(duplicates))
			x[[duplicates[i]]]=cf()
	}
	x
}	

###################################################################################
#' Remove Duplicates from Offsprings
#' 
#' Remove duplicates in \code{c(xhist,off)}, replace with non-duplicated individuals according to \code{cf}.
#'
#' @param xhist List of previous individuals
#' @param off List of offspring individuals
#' @param cf Creation function, creates random new individuals
#'
#' @return Returns \code{off} without duplicates
#'
#' @keywords internal
#' @export
###################################################################################
removeDuplicatesOffspring <- function(xhist,off,cf){
	x <- c(xhist,off)
	while(any(duplicated(x))){ 
		duplicates <- which(duplicated(x))
		for(i in 1:length(duplicates))
			x[[duplicates[i]]]=cf()
	}
	x[(length(xhist)+1):length(x)]
}	

###################################################################################
#' Tournament Selection
#' 
#' Simple Tournament Selection implementation.
#'
#' @param fitness Fitness values of individuals
#' @param tournamentSize Tournament Size
#' @param tournamentProbability Tournament Probability
#' @param selectFromN Number of tournament winners
#'
#' @return index of tournament winners
#'
#' @seealso \code{\link{combinatorialKriging}}
#'
#' @keywords internal
#' @export
###################################################################################
tournamentSelection <- function(fitness, tournamentSize, tournamentProbability, selectFromN){ 
	index <- rep(0,selectFromN)
	N <- length(fitness)
	tournamentSize <- min(tournamentSize, N) #can not select more than in population for each tournament.
	tmp <- seq(0,tournamentSize-1)
	pvec <- tournamentProbability*(1-tournamentProbability)^tmp #probabilities for individual selection
	cump <- cumsum(pvec) #cumulative probability vector
	cump[tournamentSize] <- 1 #make sure that sum is one.
	tf <- function(x){
		individuals <- sample(N,tournamentSize,FALSE,NULL)#select TSIZE individuals for a tournament, randomly.		
		fitnessrank <- order(fitness[individuals])
		rnd <- runif(1)
		i=which(cump>rnd)[1]
		individuals[fitnessrank[i]]
	}
	unlist(lapply(integer(selectFromN),tf))
}

