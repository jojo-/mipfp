# File mipfp/R/utils.R
# by Johan Barthelemy and Thomas Suesse
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

# ------------------------------------------------------------------------------
# This file provides functions to convert array to vectors and conversely as
# well as a function to remove the linearly dependent columns from a given
# matrix
# ------------------------------------------------------------------------------

Array2Vector <- function(arr) {
  # Transform a N-dimensional array a to vector, where last index moves fastest.
  #
  # Author: T. Suesse
  #
  # Args:
  #   arr: The array to be transformed into a vector.
  #
  # Returns: A vector containing the input array data.
  
  # checking if an input array is specified
  if (is.null(arr) == TRUE) {
    stop('Error: no array arr specified!')
  }
  
  dim.array <- dim(arr)
  if (is.null(dim.array) == FALSE ) {
    arr <- aperm(arr, seq(length(dim.array), 1, by = -1))    
  }
  return(c(arr))
  
}

Vector2Array <- function(vect, dim.out) {
  # Transform a vector to an array with given dimensions, where last index moves 
  # fastest.
  #
  # Author: T. Suesse
  #
  # Args:
  #   vector: The vector to be transformed into an array.
  #   dim.out: The dimension of the generated array.
  #
  # Returns: An array filled with the input vector data.
  
  # checking if an input array is specified
  if (is.null(vect) == TRUE) {
    stop('Error: no vector vect specified!')
  }
  
  # checking if an input array is specified
  if (is.null(dim.out) == TRUE) {
    stop('Error: no dimension dim.out specified for the output array!')
  }
  
  l.dim.out<-length(dim.out)
  arr <- array(vect, dim.out[seq(l.dim.out,1,by=-1)])
  arr <- aperm(arr, seq(l.dim.out, 1, by = -1)) 
  return(arr)
  
}

GetLinInd <- function(mat, tol = 1e-10) {
  # Removing the linearly dependant columns of matrix to obtain a matrix of full
  # rank using QR decomposition. See Matrix Computations from Golub and Van Loan
  # (2012) for a detailed description of the procedure.
  #
  # Author: J. Barthelemy
  #
  # Args:
  #   mat: The matrix possibly containing linearly dependant columns.
  #   tol: Rank estimation tolerance.
  #
  # Returns: The extracted columns of mat.
  
  # checking if an input matrix is specified
  if (is.null(mat) == TRUE) {
    stop('Error: mat is missing!')
  }
  
  # checking if tol is positive
  if (tol < 0.0) {
    stop('Error: tol must be positive!')
  }
  
  # QR decomposition
  mat.li <- as.matrix(mat)
  mat.qr = qr(mat.li)
  
  # if input matrix is of full rank, nothing to do
  if( mat.qr$rank == dim(mat.li)[2] ) {
    idx = seq(1,dim(mat.li)[2])
    result = list("mat.li" = mat.li, "idx" = idx)
    return(result) 
  }
  
  if( is.vector(qr.R(mat.qr)) == FALSE ) {
    diagr = abs(diag(qr.R(mat.qr)))
  } else {
    diagr = qr.R(mat.qr)[1]
  }
  
  # rank computation
  r = tail(which(diagr >= tol * diagr[1]), n = 1)
  
  # selection the r linearly independant columns
  idx <- sort(mat.qr$pivot[1:r])
  mat.li <- mat.li[,idx]  
  
  # returning the matrix and the selected row indices
  result = list("mat.li" = mat.li, "idx" = idx)
  return(result)
  
}

GetConfInt <- function(list.est, alpha = 0.05) {
  # Computing the confidence interval for the estimates produced either
  # by the Ipfp() or ObtainModelEstimates() functions.
  #
  # Author: J. Barthelemy
  #
  # Args:
  #   list.est: The list produced by either Ipfp() or ObtainModelEstimates().
  #   alpha: The confidence level.
  #
  # Returns: a list containing the confidence interval of the estimates and
  #          their probabilities
  
  # checking that a list.est is provided
  if (is.null(list.est) == TRUE)  {
    stop('Error: a list containing the estimate and their standard deviation
         is missing!')
  }
  
  # checking that alpha is in [0,1]
  if (alpha < 0.0 | alpha > 1.0) {
    stop('Error: the confidence level alpha should be in [0,1]!')
  }
  
  # checking that list.est contains the required elements
  if (is.null(list.est$x.hat) == TRUE | is.null(list.est$x.hat.se) == TRUE) {
    stop('Error: list.est does not have the component(s) x.hat and/or x.hat.se')
  }
  if (is.null(list.est$p.hat) == TRUE | is.null(list.est$p.hat.se) == TRUE) {
    stop('Error: list.est does not have the component(s) p.hat and/or p.hat.se')
  }
  
  # checking if the lenghts of the required inputs are consistent
  if (length(list.est$x.hat) != length(list.est$x.hat.se)) {
    stop('Error: lengths of x.hat and x.hat.se components are not equal!')
  }  
  if (length(list.est$p.hat) != length(list.est$p.hat.se)) {
    stop('Error: lengths of p.hat and p.hat.se components are not equal!')
  }
  
  # computing the confidence interval
  #n <- 1 / sqrt(sum(list.est$x.hat))
  l <- qnorm(1 - alpha * 0.5)  
  l < 1.96
  # ... lower bound (counts)
  ci.lo <- Array2Vector(list.est$x.hat) - l * list.est$x.hat.se
  ci.lo <- Vector2Array(ci.lo, dim(list.est$x.hat))
  dimnames(ci.lo) <- dimnames(list.est$x.hat)
  # ... upper bound (counts)
  ci.up <-Array2Vector(list.est$x.hat) + l * list.est$x.hat.se
  ci.up <- Vector2Array(ci.up, dim(list.est$x.hat))
  dimnames(ci.up) <- dimnames(list.est$x.hat)
  # ... lower bound (probabilities)  
  ci.lo.p <-Array2Vector(list.est$p.hat) - l * list.est$p.hat.se
  ci.lo.p <- Vector2Array(ci.lo.p, dim(list.est$x.hat))
  dimnames(ci.lo.p) <- dimnames(list.est$x.hat)
  # ... upper bound (probabilities)
  ci.up.p <-Array2Vector(list.est$p.hat) + l * list.est$p.hat.se
  ci.up.p <- Vector2Array(ci.up.p, dim(list.est$x.hat))
  dimnames(ci.up.p) <- dimnames(list.est$x.hat)
  
  # returning the result
  result <- list("lower.x" = ci.lo, "upper.x" = ci.up,
                 "lower.p" = ci.lo.p, "upper.p" = ci.up.p)
  return(result)
  
}
