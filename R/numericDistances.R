###################################################################################
#' Hamming Distance for Vectors
#' 
#' The number of unequal elements of two vectors (which may be of unequal length), divided by the number of elements (of the larger vector).
#'
#' @param x first vector
#' @param y second vector
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1
#'
#' @examples
#' #e.g., used for distance between bit strings
#' x <- c(0,1,0,1,0)
#' y <- c(1,1,0,0,1)
#' distanceNumericHamming(x,y)
#' p <- replicate(10,sample(c(0,1),5,replace=TRUE),simplify=FALSE)
#' distanceMatrix(p,distanceNumericHamming)
#'
#' @export
###################################################################################
distanceNumericHamming <- function(x, y){
	sum(x != y)/max(length(x),length(y))
}

###################################################################################
#' Levenshtein Distance for Numeric Vectors
#' 
#' Levenshtein distance for two numeric vectors, e.g., bit vectors.
#'
#' @param x first vector (numeric)
#' @param y second vector (numeric)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1
#'
#' @examples
#' #e.g., used for distance between bit strings
#' x <- c(0,1,0,1,0)
#' y <- c(1,1,0,0,1)
#' distanceNumericLevenshtein(x,y)
#' p <- replicate(10,sample(c(0,1),5,replace=TRUE),simplify=FALSE)
#' distanceMatrix(p,distanceNumericLevenshtein)
#'
#' @export
###################################################################################
distanceNumericLevenshtein <- function(x, y){
	#.Call("numericDistanceLevenshtein", as.numeric(x),as.numeric(y), PACKAGE="CEGO")/max(length(x),length(y))
	.Call(C_numericDistanceLevenshtein, as.numeric(x),as.numeric(y))/max(length(x),length(y))
}

###################################################################################
#' Longest Common Substring for Numeric Vectors
#' 
#' Longest common substring distance for two numeric vectors, e.g., bit vectors.
#'
#' @param x first vector (numeric)
#' @param y second vector (numeric)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1
#'
#' @examples
#' #e.g., used for distance between bit strings
#' x <- c(0,1,0,1,0)
#' y <- c(1,1,0,0,1)
#' distanceNumericLCStr(x,y)
#' p <- replicate(10,sample(c(0,1),5,replace=TRUE),simplify=FALSE)
#' distanceMatrix(p,distanceNumericLCStr)
#'
#' @export
###################################################################################
distanceNumericLCStr <- function(x, y){
	#.Call("numericDistanceLongestCommonSubstring",as.numeric(x),as.numeric(y), PACKAGE="CEGO")/max(length(x),length(y))
	.Call(C_numericDistanceLongestCommonSubstring,as.numeric(x),as.numeric(y))/max(length(x),length(y))
}