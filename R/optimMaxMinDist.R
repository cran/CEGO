
###################################################################################
#' Max-Min-Distance Optimizer
#' 
#' One-shot optimizer: Create a design with maximum sum of distances, and evaluate.
#' Best candidate is returned.
#'
#' @param x Optional set of solution(s) as a list, which are added to the randomly generated solutions and are also evaluated with the target function.
#' @param fun target function to be minimized
#' @param control (list), with the options:
#' \describe{
#' \item{\code{budget}}{ The limit on number of target function evaluations (stopping criterion) (default: 100).}
#' \item{\code{vectorized}}{ Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE.}
#' \item{\code{creationFunction}}{ Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6.}
#' \item{\code{designBudget}}{ budget of the design function \code{\link{designMaxMinDist}}, which is the number of randomly created candidates in each iteration.}
#' }
#'
#' @return a list: 	
#' \describe{
#' \item{\code{xbest}}{ best solution found}
#' \item{\code{ybest}}{ fitness of the best solution}
#' \item{\code{x}}{ history of all evaluated solutions}
#' \item{\code{y}}{ corresponding target function values f(x)}
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
#' res <- optimMaxMinDist(,lF,list(creationFunction=cF,budget=20,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res$xbest 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimEA}}, \code{\link{optimRS}}, \code{\link{optim2Opt}}
#' 
#' @export
###################################################################################
optimMaxMinDist <- function(x=NULL,fun,control=list()){
	con<-list( budget=100
						,vectorized=FALSE
						,creationFunction = solutionFunctionGeneratorPermutation(6)
						,distanceFunction = distancePermutationHamming
            ,designBudget=100
			 )
	con[names(control)] <- control
	control<-con
	
  budget <- control$budget
	vectorized <- control$vectorized
 	creationFunction <- control$creationFunction	
	
	## Create random solutions without duplicates, filled up with x
	x <- designMaxMinDist(x,creationFunction,budget,list(budget=control$designBudget,distanceFunction=control$distanceFunction))
  
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

