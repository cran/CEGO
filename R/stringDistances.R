#   Copyright (c) 2014-2016 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Hamming Distance for Strings
#' 
#' Number of unequal letters in two strings.
#'
#' @param x first string (class: character)
#' @param y second string (class: character)
#'
#' @return numeric distance value \deqn{d(x,y)}
#'
#' @examples
#' distanceStringHamming("ABCD","AACC")
#'
#' @export
#' 
###################################################################################
distanceStringHamming <- function(x, y){
	#.Call("stringDistanceHamming", as.character(x), as.character(y), PACKAGE="CEGO")
	.Call(C_stringDistanceHamming, as.character(x), as.character(y))
}

###################################################################################
#' Levenshtein Distance for Strings
#' 
#' Number of insertions, deletions and substitutions to transform one string into another
#'
#' @param x first string (class: character)
#' @param y second string (class: character)
#'
#' @return numeric distance value \deqn{d(x,y)}
#'
#' @examples
#' distanceStringLevenshtein("ABCD","AACC")
#'
#' @export
#' 
###################################################################################
distanceStringLevenshtein <- function(x, y){
	#.Call("stringDistanceLevenshtein", as.character(x), as.character(y), PACKAGE="CEGO")
	.Call(C_stringDistanceLevenshtein, as.character(x), as.character(y))
}

###################################################################################
#' Longest Common Substring distance
#' 
#' Distance between strings, based on the longest common substring.
#'
#' @param x first string (class: character)
#' @param y second string (class: character)
#'
#' @return numeric distance value \deqn{d(x,y)}
#'
#' @examples
#' distanceStringLCStr("ABCD","AACC")
#'
#' @export
#' 
###################################################################################
distanceStringLCStr <- function(x, y){
	#.Call("stringDistanceLongestCommonSubstring", as.character(x), as.character(y), PACKAGE="CEGO")
	.Call(C_stringDistanceLongestCommonSubstring, as.character(x), as.character(y))
}
