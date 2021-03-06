% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binaryOperators.R
\name{mutationBinarySingleBitFlip}
\alias{mutationBinarySingleBitFlip}
\title{Single Bit-flip Mutation for Bit-strings}
\usage{
mutationBinarySingleBitFlip(population, parameters)
}
\arguments{
\item{population}{List of bit-strings}

\item{parameters}{not used}
}
\value{
mutated population
}
\description{
Given a population of bit-strings, this function mutates all 
individuals by randomly inverting one bit in each individual. 
Due to the fixed mutation rate, this is computationally faster.
}
