#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Binary String Generator Function
#' 
#' Returns a function that generates random bit-strings of length N.
#' Can be used to create individuals of NK-Landscapes or other problems with binary representation.
#'
#' @param N length of the bit-strings
#'
#' @return returns a function, without any arguments
#'
#' @export
###################################################################################
solutionFunctionGeneratorBinary <- function(N){
	N #lazy evaluation fix, faster than force()
	function()sample(c(0,1),N,replace=TRUE)
}


###################################################################################
#' Cycle Mutation for Bit-strings
#' 
#' Given a population of bit-strings, this function mutates all 
#' individuals by cyclical shifting the string to the right or left.
#'
#' @param population List of bit-strings
#' @param parameters list of parameters: parameters$mutationRate => mutation rate, specifying number of bits flipped. Should be in range between zero and one
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationBinaryCycle <- function(population, parameters){
  mutationRate <- parameters$mutationRate
	N<-length(population[[1]])
	popsize<- length(population)
	newpop <- list()
	cmp <- max(min(round(mutationRate*N),N-1),1)
	direction <- sample(c(TRUE,FALSE),popsize,replace=T)
	for(i in 1:popsize){		
		individual <- population[[i]]
		if(direction[i])
			individual <- c(individual[-(1:cmp)],individual[1:cmp])
		else
			individual <- c(individual[N+1-(cmp:1)],individual[-(N+1-(cmp:1))])
		newpop <- c(newpop, list(individual))
	}	
	newpop
}

###################################################################################
#' Block Inversion Mutation for Bit-strings
#' 
#' Given a population of bit-strings, this function mutates all 
#' individuals by inverting a whole block, randomly selected.
#'
#' @param population List of bit-strings
#' @param parameters list of parameters: parameters$mutationRate => mutation rate, specifying number of bits flipped. Should be in range between zero and one
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationBinaryBlockInversion <- function(population, parameters){
  mutationRate <- parameters$mutationRate
	N<-length(population[[1]])
	popsize<- length(population)
	newpop <- list()
	cmp <- max(min(round(mutationRate*N),N-1),1)
	startIndex <- N-cmp
	for(i in 1:popsize){		
		index1 <- sample(startIndex,1,FALSE,NULL)
		index2 <- index1+cmp
		sel <- index1:index2
		individual <- population[[i]]
		individual[sel] <- as.numeric(!individual[sel])
		newpop <- c(newpop, list(individual))
	}	
	newpop
}


###################################################################################
#' Bit-flip Mutation for Bit-strings
#' 
#' Given a population of bit-strings, this function mutates all 
#' individuals by randomly inverting one or more bits in each individual. 
#'
#' @param population List of bit-strings
#' @param parameters list of parameters: parameters$mutationRate => mutation rate, specifying number of bits flipped. Should be in range between zero and one
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationBinaryBitFlip <- function(population, parameters){
  mutationRate <- parameters$mutationRate
	N<-length(population[[1]])
	popsize<- length(population)
	newpop <- list()
	cmp <- max(min(round(mutationRate*N),N),1)
	for(i in 1:popsize){		
		index<-sample(N,cmp,FALSE,NULL)
		individual <- population[[i]]
		individual[index]=as.numeric(!individual[index])
		newpop <- c(newpop, list(individual))
	}	
	newpop
}

###################################################################################
#' Single Bit-flip Mutation for Bit-strings
#' 
#' Given a population of bit-strings, this function mutates all 
#' individuals by randomly inverting one bit in each individual. 
#' Due to the fixed mutation rate, this is computationally faster.
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationBinarySingleBitFlip <- function(population, parameters){
	N<-length(population[[1]])
	popsize<- length(population)
	newpop <- list()
	index <- sample(N,popsize,TRUE,NULL)
	for(i in 1:popsize){		
		individual <- population[[i]]
		individual[index[i]]=as.numeric(!individual[index[i]])
		newpop <- c(newpop, list(individual))
	}	
	newpop
}


###################################################################################
#' Uniform Crossover for Bit Strings
#' 
#' Given a population of bit-strings, this function recombines each
#' individual with another individual by randomly picking bits from each parent. 
#' Note, that \code{\link{optimEA}} will not pass the whole population
#' to recombination functions, but only the chosen parents.
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationBinaryUniform <- function(population, parameters){
	N <- length(population[[1]])
	popsize <- length(population)/2  ## assumes nParents == 2
	newpop <- list()	
	for(i in 1:popsize){
		index<-sample(N,N*0.5,FALSE,NULL)
		j <- popsize + i ## second parent
		parent1 <- population[[i]]
		parent1[-index] <- population[[j]][-index] #contribution of second parent			
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}

###################################################################################
#' Single Point Crossover for Bit Strings
#' 
#' Given a population of bit-strings, this function recombines each
#' individual with another individual by randomly specifying a single position.
#' Information before that position is taken from the first parent,
#' the rest from the second.
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationBinary1Point <- function(population, parameters){
	N <- length(population[[1]])
	popsize <- length(population)/2  ## assumes nParents == 2
	newpop <- list()	
	for(i in 1:popsize){
		index <- sample(2:(N-1),1,FALSE,NULL)
		inds <- 1:index
		j <- popsize + i ## second parent
		parent1 <- population[[i]] #first parent
		parent1[inds] <- population[[j]][inds] #contribution of second parent			
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}

###################################################################################
#' Two Point Crossover for Bit Strings
#' 
#' Given a population of bit-strings, this function recombines each
#' individual with another individual by randomly specifying 2 positions.
#' Information in-between is taken from one parent, the rest from the other.
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationBinary2Point <- function(population, parameters){
	N <- length(population[[1]])
	popsize <- length(population)/2  ## assumes nParents == 2
	newpop <- list()	
	for(i in 1:popsize){
		index1 <- sample(N,1,FALSE,NULL)
		index2 <- sample(2:(N-1),1,FALSE,NULL)
		inds <- index1:index2
		j <- popsize + i ## second parent
		parent1 <- population[[i]] #first parent
		parent1[inds] <- population[[j]][inds] #contribution of second parent			
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}

###################################################################################
#' Arithmetic (AND) Crossover for Bit Strings
#' 
#' Given a population of bit-strings, this function recombines each
#' individual with another individual by computing \code{parent1 & parent2} (logical AND).
#'
#' @param population List of bit-strings
#' @param parameters not used
#'
#' @return population of recombined offspring
#'
#' @export
###################################################################################
recombinationBinaryAnd <- function(population, parameters){
	popsize <- length(population)/2  ## assumes nParents == 2
	newpop <- list()	
	for(i in 1:popsize){
		j <- popsize + i ## second parent
		parent1 <- population[[i]] * population[[j]] #logical and (in 0/1 encoding)
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}
