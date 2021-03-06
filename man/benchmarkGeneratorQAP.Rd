% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutationBenchmarkFunctions.R
\name{benchmarkGeneratorQAP}
\alias{benchmarkGeneratorQAP}
\title{Create Quadratic Assignment Problem (QAP) Benchmark}
\usage{
benchmarkGeneratorQAP(a, b)
}
\arguments{
\item{a}{distance matrix}

\item{b}{flow matrix}
}
\value{
the function of type cost=f(permutation)
}
\description{
Creates a benchmark function for the Quadratic Assignment Problem.
}
\examples{
set.seed(1)
n=5
#ceate a flow matrix
A <- matrix(0,n,n) 
for(i in 1:n){
	for(j in i:n){
		if(i!=j){
			A[i,j] <- sample(100,1)
			A[j,i] <- A[i,j]
 	}
	}
}
#create a distance matrix
locations <- matrix(runif(n*2)*10,,2)
B <- as.matrix(dist(locations))
#create QAP objective function 
fun <- benchmarkGeneratorQAP(A,B)
#evaluate
fun(1:n)
fun(n:1)

}
\seealso{
\code{\link{benchmarkGeneratorFSP}}, \code{\link{benchmarkGeneratorTSP}}, \code{\link{benchmarkGeneratorWT}}
}
