
###################################################################################
#' Levenshtein Distance forsSequences of numbers
#' 
#' Levenshtein distance for two sequences of numbers
#'
#' @param x first vector (numeric vector)
#' @param y second vector (numeric vector)
#'
#' @return numeric distance value \deqn{d(x,y)}
#'
#' @examples
#' #e.g., used for distance between integer sequence
#' x <- c(0,1,10,2,4)
#' y <- c(10,1,0,4,-4)
#' distanceSequenceLevenshtein(x,y)
#' p <- replicate(10,sample(1:5,3,replace=TRUE),simplify=FALSE)
#' distanceMatrix(p,distanceSequenceLevenshtein)
#'
#' @export
###################################################################################
distanceSequenceLevenshtein <- function(x, y){
  .Call(C_numericDistanceLevenshtein, as.numeric(x),as.numeric(y))
}