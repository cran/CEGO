###################################################################################
#' String Generator Function
#' 
#' Returns a function that generates random strings of length N, with given letters.
#' Can be used to generate individual solutions for permutation problems, e.g., Travelling Salesperson Problem
#'
#' @param N length of the permutations returned
#' @param lts letters allowed in the string
#'
#' @return returns a function, without any arguments
#'
#' @examples
#' fun <- solutionFunctionGeneratorString(10,c("A","C","G","T"))
#' fun()
#' fun()
#' fun()
#' 
#' @export
###################################################################################
solutionFunctionGeneratorString <- function(N,lts=c("A","C","G","T")){
	N #lazy evaluation fix, faster than force()
	function()paste(sample(lts,N,replace = TRUE),collapse="")
}

###################################################################################
#' Mutation for Strings
#' 
#' Given a population of strings, this function mutates all 
#' individuals by randomly changing an element of the string.
#'
#' @param population List of permutations
#' @param parameters list of parameters, with \code{parameters$mutationRate} and \code{parameters$lts}.
#' \code{parameters$mutationRate} should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of interchanges
#' performed, relative to the permutation length (N). 0 means none. 1 means N interchanges.
#' The default is 1/N. \code{parameters$lts} are the possible letters in the string.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationStringRandomChange <- function(population, parameters=list()){
  N <- nchar(population[[1]])	
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N
  if(is.null(parameters$lts)) parameters$lts <- c("A","C","G","T")
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	lts <- parameters$lts
	if(mutations==0)
		return(population)
	samples <- mutations * popsize
	index <- sample.int(N,samples,TRUE,NULL) 
	newLetter <- sample(lts,samples,TRUE,NULL)
	newpop <- list()	
	for(i in 1:popsize){
		individual <- population[[i]]
		if(mutations == 1){
		  ind <- index[i]
			substr(individual,ind,ind) <- newLetter[i]
		}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				ind <- index[jj]
				substr(individual,ind,ind) <- newLetter[jj]
			}	
		}		
		newpop <- c(newpop, list(individual))
	}
	newpop
}

###################################################################################
#' Single Point Crossover for Strings
#' 
#' Given a population of strings, this function recombines each
#' individual with another random individual.
#' Note, that \code{\link{optimEA}} will not pass the whole population
#' to recombination functions, but only the chosen parents.
#'
#' @param population List of strings
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationStringSinglePointCrossover <- function(population, parameters){
  N <- nchar(population[[1]])	
	popsize <- length(population)/2
	newpop <- list()	
	index <- sample.int(N-1,popsize,TRUE,NULL) #index for switch between parent 1 and 2 for each operation
	for(i in 1:popsize){
		j <- popsize + i ## second parent
		ind <- index[i]
		parent1 <- population[[i]]
		parent2 <- population[[j]]
		child <- paste(substr(parent1,1,ind),substr(parent2,ind+1,N),sep="")
		newpop <- c(newpop, child)
	}	
	newpop
}