#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Permutation Generator Function
#' 
#' Returns a function that generates random permutations of length N.
#' Can be used to generate individual solutions for permutation problems, e.g., Travelling Salesperson Problem
#'
#' @param N length of the permutations returned
#'
#' @return returns a function, without any arguments
#'
#' @export
###################################################################################
solutionFunctionGeneratorPermutation <- function(N){
	N #lazy evaluation fix, faster than force()
	function()sample(1:N,replace=FALSE)
}

###################################################################################
#' Interchange Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly interchanging two arbitrary elements of the permutation.
#'
#' @param population List of permutations
#' @param parameters not used.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationInterchange <- function(population, parameters){
	N=length(population[[1]])
	popsize= length(population)
	newpop <- list()
	index1 <- sample(N,popsize,TRUE,NULL) #todo not allowed on cran afaik
	index2 <- sample(N,popsize,TRUE,NULL)
	for(i in 1:popsize){				
		individual <- population[[i]]
		val1= individual[index1[i]]
		individual[index1[i]]= individual[index2[i]]
		individual[index2[i]]= val1
		newpop <- c(newpop, list(individual))
	}	
	newpop
}


###################################################################################
#' Swap Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly interchanging two adjacent elements of the permutation.
#'
#' @param population List of permutations
#' @param parameters not used
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationSwap <- function(population, parameters){
	N=length(population[[1]])
	popsize= length(population)
	newpop <- list()
	index1 <- sample(N-1,popsize,TRUE,NULL) #todo not allowed on cran afaik
	index2 <- index1 +1
	for(i in 1:popsize){				
		individual <- population[[i]]
		val1= individual[index1[i]]
		individual[index1[i]]= individual[index2[i]]
		individual[index2[i]]= val1
		newpop <- c(newpop, list(individual))
	}
	newpop
}


###################################################################################
#' Cycle Crossover (CX) for Permutations
#' 
#' Given a population of permutations, this function recombines each
#' individual with another random individual.
#' Note, that \code{\link{optimEA}} will not pass the whole population
#' to recombination functions, but only the chosen parents.
#'
#' @param population List of permutations
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationPermutationCycleCrossover <- function(population, parameters){
	N=length(population[[1]])
	popsize= length(population)
	newpop <- list()	
	for(i in 1:popsize){
		j <- (1:popsize)[-i][sample(popsize-1,1,FALSE,NULL)] #draw second parent
		parent1 <- population[[i]]
		parent2 <- population[[j]]
		e1 <- parent1[1]
		e2 <- parent2[1]		
		parent2[1] <- e1
		while(e1 != e2){
			e1 <- e2
			rplc <- which(parent1==e1)
			e2 <- parent2[rplc]
			parent2[rplc] <- e1
		}		
		newpop <- c(newpop, list(parent2))
	}	
	newpop
}