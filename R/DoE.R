###################################################################################
#' Random Design
#' 
#' Create a random initial population or experimental design, given a specifed creation function,
#' as well as a optional set of user-specified design members and a maximum design size.
#' Also removes duplicates from the design/population.
#'
#' @param x Optional list of user specified solutions to be added to the design, defaults to NULL
#' @param cf Creation function, creates random new individuals
#' @param size size of the design
#' @param control not used
#'
#' @return Returns list with experimental design without duplicates
#'
#' @seealso \code{\link{optimRS}}, \code{\link{designMaxMinDist}}
#' @keywords internal
#' @examples
#' # Create a design of 10 permutations, each with 5 elements
#' design <- designRandom(NULL,function()sample(5),10)
#' # Create a design of 20 real valued 2d vectors
#' design <- designRandom(NULL,function()runif(2),20)
#' @export
###################################################################################
designRandom <- function(x=NULL,cf,size,control=list()){
	## initialization
	if(is.null(x)){
    x <- list()
    k=0
  }else{ #given start population
    k=length(x)
  }
		
	if(k>size){
		x <- x[1:size]
	}else if(k<size){
		## CREATE initial population
		x <- c(x, replicate(size-k , cf(),simplify=FALSE))
	}#else if k==size do nothing.
		
	## REPLACE duplicates from initial population with unique individuals
	x <- removeDuplicates(x, cf)
}

###################################################################################
#' Max-Min-Distance Design
#' 
#' Build a design of experiments in a sequential manner: First candidate solution is created at random.
#' Afterwards, candidates are added sequentially, maximizing the minimum distances to the existing candidates.
#' Each max-min problem is resolved by random sampling.
#' The aim is to get a rather diverse design.
#'
#' @param x Optional list of user specified solutions to be added to the design/population, defaults to NULL
#' @param cf Creation function, creates random new individuals
#' @param size size of the design
#' @param control list of controls. \code{control$distanceFunction} requires a distance function to compare two candidates created by cf. 
#' \code{control$budget} is the number of candidates for the random sampling, defaults to 100.
#'
#' @return Returns list with experimental design without duplicates
#'
#' @seealso \code{\link{optimMaxMinDist}}, \code{\link{designRandom}}
#' @keywords internal
#' @examples
#' # Create a design of 10 permutations, each with n=5 elements, 
#' # and with 50 candidates for each sample.
#' # Note, that in this specific case the number of candidates 
#' # should be no larger than factorial(n).
#' # The default (hamming distance) is used.
#' design <- designMaxMinDist(NULL,function()sample(5),10,
#' 		control=list(budget=50))
#' # Create a design of 20 real valued 2d vectors, 
#' # with 100 candidates for each sample
#' # using euclidean distance.
#' design <- designMaxMinDist(NULL,function()runif(2),20,
#'		control=list(budget=100,
#'		distanceFunction=function(x,y)sqrt(sum((x-y)^2))))
#' # plot the resulting design
#' plot(matrix(unlist(design),,2,byrow=TRUE))
#' @export
###################################################################################
designMaxMinDist <- function(x=NULL,cf,size,control=list()){
	con<-list(budget=100, 
            distanceFunction=distancePermutationHamming
			 )
	con[names(control)] <- control
	control<-con

	## initialization
	if(is.null(x)){
    x <- list(cf())
    k=1
  }else{ #given start population
    k=length(x)
  }		
    
	if(k>size){
		x <- x[1:size]
	}else if(k<size){
		## CREATE initial population
    for(ki in (k+1):size){
			fun <- function(xnew) -min(distanceVector(xnew,x,control$distanceFunction))
			res <- optimRS(,fun,control=list(creationFunction=cf,budget=control$budget))
      x <- c(x,list(res$xbest))
    }
	}#else if k==size do nothing.
	## REPLACE duplicates from initial population with unique individuals
	x <- removeDuplicates(x, cf)
}
