
###################################################################################
#' Two-Opt
#' 
#' Implementation of a Two-Opt local search.
#'
#' @param x start solution of the local search
#' @param fun function that determines cost or length of a route/permutation
#' @param control (list), with the options:
#' \describe{
#' \item{\code{archive}}{ Whether to keep all candidate solutions and their fitness in an archive (TRUE) or not (FALSE). Default is TRUE.}
#' \item{\code{budget}}{ The limit on number of target function evaluations (stopping criterion) (default: 100)}
#' \item{\code{creationFunction}}{ Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6}
#' \item{\code{vectorized}}{ Boolean. Defines whether target function is vectorized (takes a list of solutions 
#' as argument) or not (takes single solution as argument). Default: FALSE}
#' }
#'
#' @return a list with:
#' \describe{
#' \item{\code{xbest}}{ best solution found}
#' \item{\code{ybest}}{ fitness of the best solution}
#' \item{\code{count}}{ number of performed target function evaluations }
#' }
#' 
#' @examples
#' seed=0
#' #distance
#' dF <- distancePermutationHamming
#' #creation
#' cF <- function()sample(5)
#' #objective function
#' lF <- landscapeGeneratorUNI(1:5,dF)
#' #start optimization
#' set.seed(seed)
#' res <- optim2Opt(,lF,list(creationFunction=cF,budget=100,
#'    vectorized=TRUE)) ##target function is "vectorized", expects list of solutions as input
#' res
#'
#' @references Wikipedia contributors. "2-opt." Wikipedia, The Free Encyclopedia. Wikipedia, The Free Encyclopedia, 13 Jun. 2014. Web. 21 Oct. 2014. 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimEA}}, \code{\link{optimRS}}, \code{\link{optimMaxMinDist}} 
#' 
#' @export
###################################################################################
optim2Opt <- function(x=NULL,fun,control=list()){ 
	con<-list(budget=100
					, vectorized=FALSE
					, creationFunction = solutionFunctionGeneratorPermutation(6)
          , archive =TRUE
			 )
	con[names(control)] <- control
	control <- con

  archive <- control$archive
  budget <- control$budget
  vectorized <- control$vectorized
	creationFunction <- control$creationFunction	
	
  if(is.null(x)){
    route=creationFunction()
  }else{ #given start population
    route=x[[1]]
  }  

	improvement=TRUE
	bestRoute=route
  if (vectorized)
    bestDist = fun(list(route))
  else
    bestDist = fun(route)
  if(archive){
    fithist <- bestDist
    xhist <- route
  } 
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
        if(archive){
          xhist <- append(xhist,newRoute)  
          fithist <-  c(fithist, newDist)
        }     
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
  if(archive)
    return(list(xbest=bestRoute,ybest=bestDist,x=xhist,y=fithist, count=count))
  else
    return(list(xbest=bestRoute,ybest=bestDist,count=count))
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

