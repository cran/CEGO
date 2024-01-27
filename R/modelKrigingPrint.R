
###################################################################################
#' Print Function: modelKriging
#'
#' Print information about a Kriging fit, as produced by \code{\link{modelKriging}}.
#'
#' @rdname print
#' @method print modelKriging
#  @S3method print modelKriging
#' @param x	fit returned by \code{\link{modelKriging}}.
#' @param ... additional parameters	
#' @export
#' @keywords internal
###################################################################################
print.modelKriging <- function(x,...){
	cat("------------------------\n")
	cat("Kriging model fit of class modelKriging\n")
	cat("\n")
	if(!is.null(x$theta)){
		cat("Estimated parameters of the correlation function (kernel):\n")
		cat(paste(round(x$theta,4),sep=" "))
		cat("\n \n")
	}
	if(!is.null(x$distanceParameters)){
		cat("Estimated parameters of the distance function:\n")
		cat(paste(round(x$distanceParameters,4),sep=" "))
		cat("\n \n")
	}
	if(x$combineDistances & is.list(x$distanceFunction)){
		cat("Several distance functions are combined in this model.\n")
		cat("The following weights are used to combine them:\n")
		cat(paste(round(x$distanceWeights,4),sep=" ",collaps=" "))
		cat("\n \n")
	}		
	if(x$useLambda){
		cat("Estimated regularization constant (nugget) lambda:\n")	
		cat(x$lambda)
		cat("\n \n")
	}else{
		cat("No regularization constant (nugget) was used.\n")	
		cat("\n")
	}
	cat("Number of Likelihood evaluations during MLE:\n")	
	cat(x$nevals)
	cat("\n")
	cat("Minimal NegLnLikelihood:\n")	
	cat(x$like)
	cat("\n")
	cat("------------------------\n")
}

