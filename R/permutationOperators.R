#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

#TODO:
# insert/shift mutation
# reversal mutation

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
#' @examples
#' fun <- solutionFunctionGeneratorPermutation(10)
#' fun()
#' fun()
#' fun()
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
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of interchanges
#' performed, relative to the permutation length (N). 0 means none. 1 means N interchanges.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationInterchange <- function(population, parameters=list()){
  N <- length(population[[1]])	
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations<=0)
		return(population)
	samples <- mutations * popsize
	index1 <- sample.int(N,samples,TRUE,NULL) 
	index2 <- sample.int(N,samples,TRUE,NULL)
	mutationPermutationInterchangeCore(population,popsize,mutations,index1,index2)
}
#mutationPermutationInterchange(list(1:5),list(mutationRate=1))

###################################################################################
#' Interchange of permutation elements
#' 
#' Support function for \code{\link{mutationPermutationInterchange}} and \code{\link{mutationPermutationSwap}}.
#'
#' @param population List of permutations
#' @param popsize population size
#' @param mutations number of mutated elements for each individual
#' @param index1 vector of first indices, one element for each interchange
#' @param index2 vector of second indices, one element for each interchange
#'
#' @return mutated population
#'
#' @keywords internal
###################################################################################
mutationPermutationInterchangeCore <- function(population,popsize,mutations,index1,index2){
	newpop <- list()
	for(i in 1:popsize){
		individual <- population[[i]]
		if(mutations == 1){
			val1 <- individual[index1[i]]
			individual[index1[i]] <- individual[index2[i]]
			individual[index2[i]] <- val1
		}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				i1 <- index1[jj]
				i2 <- index2[jj]
				val1= individual[i1]
				individual[i1]= individual[i2]
				individual[i2]= val1
			}	
		}		
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
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of swaps
#' performed, relative to the permutation length (N). 0 means none. 1 means N swaps.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationSwap <- function(population,parameters=list()){
	N <- length(population[[1]])
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations<=0)
		return(population)
	samples <- mutations * popsize
	index1 <- sample.int(N-1,samples,TRUE,NULL) 
	index2 <- index1 +1
	mutationPermutationInterchangeCore(population,popsize,mutations,index1,index2)
}
#mutationPermutationSwap(list(1:5),list(mutationRate=1))


###################################################################################
#' Reversal Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly selecting two indices, and reversing the respective sub-permutation.
#'
#' @param population List of permutations
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of reversals
#' performed, relative to the permutation length (N). 0 means none. 1 means N reversals.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationReversal <- function(population, parameters=list()){
	N <- length(population[[1]])
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N  
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations<=0)
		return(population)
	samples <- mutations * popsize
	newpop <- list()
	index1 <- sample.int(N,samples,TRUE,NULL) 
	index2 <- sample.int(N,samples,TRUE,NULL)
	for(i in 1:popsize){				
		individual <- population[[i]]
		if(mutations == 1){
			individual[index1[i]:index2[i]] <- individual[index2[i]:index1[i]]
		}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				i1 <- index1[jj]
				i2 <- index2[jj]
				individual[i1:i2] <- individual[i2:i1]
			}	
		}		
		newpop <- c(newpop, list(individual))
	}	
	newpop
}

###################################################################################
#' Insert Mutation for Permutations
#' 
#' Given a population of permutations, this function mutates all 
#' individuals by randomly selecting two indices.
#' The element at index1 is moved to positition index2, other elements
#' 
#'
#' @param population List of permutations
#' @param parameters list of parameters, currently only uses parameters$mutationRate, 
#' which should be between 0 and 1 (but can be larger than 1). The mutation rate determines the number of reversals
#' performed, relative to the permutation length (N). 0 means none. 1 means N reversals.
#' The default is 1/N.
#'
#' @return mutated population
#'
#' @export
###################################################################################
mutationPermutationInsert <- function(population, parameters=list()){
	N <- length(population[[1]])
  if(is.null(parameters$mutationRate)) parameters$mutationRate <- 1/N  
	mrate <- parameters$mutationRate
	popsize <- length(population)
	mutations <- ceiling(N * mrate)
	if(mutations<=0)
		return(population)
	samples <- mutations * popsize
	newpop <- list()
	index1 <- sample.int(N,samples,TRUE,NULL) 
	index2 <- sample.int(N,samples,TRUE,NULL)
	for(i in 1:popsize){				
		individual <- population[[i]]
		#if(mutations == 1){
		#	individual[index1[i]:index2[i]] <- individual[index2[i]:index1[i]]
		#}else{
			j <- ((i-1)*mutations+1) : (i*mutations)
			for(jj in j){
				i1 <- index1[jj]
				i2 <- index2[jj]
				#print(i1)
				#print(i2)
				ind1 <- individual[i1] #move out
				individual <- individual[-i1]
				if(i2==1) 
					individual <- c(ind1,individual) #and insert
				else if(i2==N)	
					individual <- c(individual,ind1)  #and insert
				else
					individual <- c(individual[1:(i2-1)],ind1,individual[i2:(N-1)]) #and insert
			}	
		#}		
		newpop <- c(newpop, list(individual))
	}	
	newpop
}
#mutationPermutationInsert(list(1:5))

###################################################################################
#' Cycle Crossover (CX) for Permutations
#' 
#' Given a population of permutations, this function recombines each
#' individual with another individual.
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
	popsize <- length(population)/2  ## assumes nParents == 2
	newpop <- list()	
	for(i in 1:popsize){
		j <- popsize + i ## second parent
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

###################################################################################
#' Order Crossover 1 (OX1) for Permutations
#' 
#' Given a population of permutations, this function recombines each
#' individual with another individual.
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
recombinationPermutationOrderCrossover1 <- function(population, parameters){
	popsize <- length(population)/2  ## assumes nParents == 2
	N <- length(population[[1]]) #number of elements
	newpop <- list()	
	for(i in 1:popsize){
		j <- popsize + i ## second parent
		parent1 <- population[[i]]
		parent2 <- population[[j]]
		idx <- sort(sample(N,2)	)## select part from first parent
		pnew <- setdiff(parent2,parent1[idx[1]:idx[2]])## identify parts from second parent which are not in the selected part
		parent1[-(idx[1]:idx[2])] <- pnew #insert
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}

###################################################################################
#' Position Based Crossover (POS) for Permutations
#' 
#' Given a population of permutations, this function recombines each
#' individual with another individual.
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
recombinationPermutationPositionBased <- function(population, parameters){
	popsize <- length(population)/2  ## assumes nParents == 2
	N <- length(population[[1]]) #number of elements
	newpop <- list()	
	for(i in 1:popsize){
		j <- popsize + i ## second parent
		parent1 <- population[[i]]
		parent2 <- population[[j]]
		idx <- sample(N,N/2,replace=FALSE)
		parent1[-idx] <- setdiff(parent2,parent1[idx])
		newpop <- c(newpop, list(parent1))
	}	
	newpop
}

###################################################################################
#' Alternating Position Crossover (AP) for Permutations
#' 
#' Given a population of permutations, this function recombines each
#' individual with another individual.
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
recombinationPermutationAlternatingPosition <- function(population, parameters){
	popsize <- length(population)/2  ## assumes nParents == 2
	N <- length(population[[1]]) #number of elements
	newpop <- list()	
	for(i in 1:popsize){
		j <- popsize + i ## second parent
		parent1 <- population[[i]]
		parent2 <- population[[j]]
		pnew <- NULL
		for(i in 1:N){
			if(i%%2==1)
				pnew <- c(pnew,setdiff(parent1,pnew)[1])		
			if(i%%2==0)
				pnew <- c(pnew,setdiff(parent2,pnew)[1])		
		}
		newpop <- c(newpop, list(pnew))
	}	
	newpop
}
