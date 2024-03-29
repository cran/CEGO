################################################################################
##	This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

################################################################################
#' Combinatorial Efficient Global Optimization
#'
#' Model building, surrogate model
#' based optimization and Efficient Global Optimization in combinatorial
#' or mixed search spaces. This includes methods for distance calculation,
#' modeling and handling of indefinite kernels/distances.
#'
#' \tabular{ll}{
#' Package: \tab CEGO\cr
#' Type: \tab Package\cr
#' Version: \tab 2.4.3\cr
#' Date: \tab 2024-01-27\cr
#' License: \tab GPL (>= 3)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
# @name CEGO-package
# @aliases CEGO
# @docType package
#' @title Combinatorial Efficient Global Optimization in R
#' @author Martin Zaefferer \email{mzaefferer@@gmail.com} 
#' @references Zaefferer, Martin; Stork, Joerg; Friese, Martina; Fischbach, Andreas; Naujoks, Boris; Bartz-Beielstein, Thomas. (2014). Efficient global optimization for combinatorial problems. In Proceedings of the 2014 conference on Genetic and evolutionary computation (GECCO '14). ACM, New York, NY, USA, 871-878. DOI=10.1145/2576768.2598282
#' @references Zaefferer, Martin; Stork, Joerg; Bartz-Beielstein, Thomas. (2014). Distance Measures for Permutations in Combinatorial Efficient Global Optimization. In Parallel Problem Solving from Nature - PPSN XIII (p. 373-383). Springer International Publishing.
#' @references Zaefferer, Martin and Bartz-Beielstein, Thomas (2016). Efficient Global Optimization with Indefinite Kernels. Parallel Problem Solving from Nature-PPSN XIV. Accepted, in press. Springer. 
#' @keywords package
#' @seealso Interface of main function: \code{\link{optimCEGO}} 
#' @import MASS
#' @import graphics
#' @import stats
#' @import DEoptim
#' @importFrom quadprog solve.QP
#' @importFrom Matrix nearPD
#' @importFrom methods formalArgs
#' @importFrom anticlust balanced_clustering
#' @useDynLib CEGO, .registration = TRUE, .fixes = "C_"
#' 
#' @section Acknowledgments:
#' This work has been partially supported by the Federal Ministry of Education
#' and Research (BMBF) under the grants CIMO (FKZ 17002X11) and
#' MCIOP (FKZ 17N0311).
#'
"_PACKAGE" #ends description
################################################################################