
###################################################################################
#' Combinatorial Random Search
#' 
#' Random Search for mixed or combinatorial optimization. Solutions are generated completely at random.
#'
#' @param x Optional set of solution(s) as a list, which are added to the randomly generated solutions and are also evaluated with the target function.
#' @param fun target function to be minimized
#' @param control (list), with the options:
#' \describe{
#' \item{\code{budget}}{ The limit on number of target function evaluations (stopping criterion) (default: 100)}
#' \item{\code{vectorized}}{ Boolean. Defines whether target function is vectorized (takes a list of solutions as argument) or not (takes single solution as argument). Default: FALSE}
#' \item{\code{creationFunction}}{ Function to create individuals/solutions in search space. Default is a function that creates random permutations of length 6}
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
#' res <- optimRS(,lF,list(creationFunction=cF,budget=100,
#'		vectorized=TRUE)) ##target function is "vectorized", expects list as input
#' res$xbest 
#'
#' @seealso \code{\link{optimCEGO}}, \code{\link{optimEA}}, \code{\link{optim2Opt}}, \code{\link{optimMaxMinDist}} 
#' 
#' @export
###################################################################################
optimRS <- function(x=NULL,fun,control=list()){
	con<-list( budget=100
						,vectorized=FALSE
						,creationFunction = solutionFunctionGeneratorPermutation(6)
			 )
	con[names(control)] <- control
	control<-con
	
	## Create random solutions without duplicates, filled up with x
	x <- designRandom(x,control$creationFunction,control$budget)
  
  #evaluate
  if(control$vectorized) 
		y <- fun(x)
	else
		y <- unlist(lapply(x,fun))
  
  #best value found:  
	j <- which.min(y)
  
  #return
	list(xbest=x[[j]],ybest=y[j],x=x,y=y, count=control$budget)
}

