# File mipfp/R/ipfpMultiDim.R
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
# This file provides the functions Ipfp and Ipfp.Covar that respectively
# implements the iterative proportional fitting procedure and a Delta method
# that computes the covariance matrix of the estimate produced by Ipfp.
# ------------------------------------------------------------------------------

library(cmm)

Ipfp <- function(seed, target.list, target.data, print = FALSE, iter = 1000, 
                 tol = 1e-10, na.target = FALSE, compute.cov = FALSE, ...) {
  # Update an array using the iterative proportional fitting procedure.
  #
  # Author: J. Barthelemy
  #  
  # Args:
  #   seed: The initial multi-dimensional array to be updated. Each cell must
  #         be non-negative.
  #   target.list: A list of the target margins provided in target.data. Each
  #                component of the list is an array whose cells indicates
  #                 which dimension the corresponding margin relates to.
  #   target.data: A list containing the data of the target margins. Each
  #                component of the list is an array storing a margin.
  #                The list order must follow the one defined in target.list. 
  #                Note that the cells of the arrays must be non-negative, but
  #                can contains NA values.
  #   print: Verbose parameter: if TRUE prints the current iteration number
  #          and the value of the stopping criterion.
  #   iter: The maximum number of iteration allowed; must be greater than 0.
  #   tol: If the maximum absolute difference between two iteration is lower
  #        than the value specified by tol, then ipfp has reached convergence
  #        (stopping criterion); must be greater than 0.
  #   na.target: If set to TRUE, allows the targets to have NA cells. In that
  #              case the margins consistency is not checked.
  #   compute.cov: If set to TRUE, then the function also return the covariance
  #                matrices of the updated cells and cells proportion.
  #   ...: Additional arguments that can be passed to the IpfpCovar function
  #        if compute.cov = TRUE. See IpfpCovar documentation.
  #
  # Returns: A list generated whose elements are
  #   xi.hat: An array of the same dimension of seed whose margins match the
  #           ones specified in target.list.
  #   stp.crit: The final value of the stopping criterion.
  #   evol.stp.crit: Evolution of the stopping criterion over the iterations.
  #   conv: A boolean indicating whether the algorithm converged to a solution.
  #   check.margins: A list returning, for each margin, the absolute maximum 
  #                 deviation between the desired and generated margin.
  
  # checking if a seed is provided
  if (is.null(seed) == TRUE) {
    stop('Error: no seed specified!')
  }
  
  # checking if target are provided
  if (is.null(target.data) == TRUE | is.null(target.data) == TRUE) {
    stop('Error: target.data and/or target.data not specified!')
  }
  
  # checking if NA in target cells if na.target is set to FALSE
  if (is.na(min(sapply(target.data, min))) == TRUE & na.target == FALSE)  {
    stop('Error: NA values present in the margins - use na.target = TRUE!')
  }
  
  # checking if NA in seed
  if (is.na(min(seed)) == TRUE) {
    stop('Error: NA values present in the seed!')
  }
  
  # checking non negativity condition for the seed and the target
  if (min(sapply(target.data, min), na.rm = na.target) < 0 | min(seed) < 0) {
    stop('Error: Target and Seed cells must be non-negative!')    
  }  
  
  # checking the strict positiviy of tol and iter
  if (iter < 1 | tol <= 0.0) {
    stop('Error: tol and iter must be strictly positive!')
  }
  
  # checking if NA allowed and requesting the covariance matrices
  if (na.target == TRUE & compute.cov == TRUE) {
    warning('Missing values allowed in the target margins.
             Computation of the covariance matrices set to FALSE!')
    compute.cov <- FALSE
  }
  
  # checking the margins consistency if no missing values in the targets
  check.margins <- TRUE  
  if (na.target == FALSE) {
    if (length(target.data) > 1) {
      for (m in 2:length(target.data)) {      
        if (abs(sum(target.data[[m-1]]) - sum(target.data[[m]])) > 1e-10) {
          check.margins <- FALSE
          warning('Target not consistents - shifting to probabilities!
                  Check input data!\n')
          break
        }      
      }
    }
  } else {
    if (print == TRUE) {
      cat('NOTE: Missing values present in target cells. ')
      cat('Margins consistency not checked!\n')  
    }        
  }
  
  # if margins are not consistent, shifting from frequencies to probabilities
  if (check.margins == FALSE) {
    seed <- seed / sum(seed)
    for (m in 1:length(target.data)) {
      target.data[[m]] <- target.data[[m]] / sum(target.data[[m]])
    }
  }
  
  if (print == TRUE & check.margins == TRUE & na.target == FALSE) {
    cat('Margins consistency checked!\n')
  }
    
  # initial value is the seed
  result <- seed  
  converged <- FALSE
  tmp.evol.stp.crit <- vector(mode="numeric", length = iter)
  
  # ipfp iterations
  for (i in 1:iter) {
    
    if (print) {
      cat('... ITER', i, '\n')
    } 
    
    # saving previous iteration result (for testing convergence)
    result.temp <- result
            
    # loop over the constraints
    for (j in 1:length(target.list)) {
      # ... extracting current margins
      temp.sum <- apply(result, target.list[[j]], sum)
      # ... computation of the update factor, taking care of 0 and NA cells   
      update.factor <- ifelse(target.data[[j]] == 0 | temp.sum == 0, 0,
                              target.data[[j]] / temp.sum)
      if (na.target == TRUE) {
        update.factor[is.na(update.factor)] <- 1;
      }
      # ... apply the update factor
      result <- sweep(result,target.list[[j]], update.factor, FUN = "*")
    }
    
    # stopping criterion
    stp.crit <- max(abs(result - result.temp))
    tmp.evol.stp.crit[i] <- stp.crit
    if (stp.crit < tol) {
      converged <- TRUE
      if (print) {
        cat('Convergence reached after', i, 'iterations!\n')
      } 
      break
    }
    
    if (print) {
      cat('       stoping criterion:', stp.crit, '\n')
    }
    
  }
    
  # checking the convergence
  if (converged == FALSE) {
    warning('IPFP did not converged after ', iter, ' iteration(s)! 
            This migh be due to 0 cells in the seed, maximum number 
            of iteration too low or tolerance too small\n')
  }        
  
  # computing final max difference between generated and target margins
  diff.margins <- vector(mode = "numeric", length = length(target.list))
  if (na.target == FALSE) {
    for (j in 1:length(target.list)) {
      diff.margins[j] = max(abs(target.data[[j]] 
                                - apply(result, target.list[[j]], sum))) 
    }
  }
  
  # storing the evolution of the stopping criterion
  evol.stp.crit <- tmp.evol.stp.crit[1:i]
  
  # computing the proportions
  result.prop <- result / sum(result)
  
  # gathering the results in a list
  results.list <- list("x.hat" = result, "p.hat" = result.prop, 
                      "stp.crit" = stp.crit, "conv" = converged, 
                      "check.margins" = diff.margins, 
                      "evol.stp.crit" = evol.stp.crit);
  
  # adding covariance if requested
  if (compute.cov == TRUE) {
    results.list$p.hat.cov <- IpfpCov(result, seed, target.list, ...)
    results.list$x.hat.cov <- results.list$p.hat.cov * sum(result)^2
    results.list$p.hat.se <- sqrt(diag(results.list$p.hat.cov))
    results.list$x.hat.se <- sqrt(diag(results.list$x.hat.cov))
  }
  
  # returning the result
  return(results.list)
  
}

IpfpCov <- function(estimate, seed, target.list, replace.zeros = 1e-10) {
  # Compute the covariance matrix of the estimators produced by Ipfp.
  #
  # This function determines the covariance matrix of the estimated proportions
  # using the Delta method given in the paper "Models for Contingency Tables 
  # With Known Margins When Target and Sampled Populations Differ" written by 
  # Little and Wu (1991).
  #
  # Author: J. Barthelemy
  #
  # Args:
  #   estimate: The array of estimate produced by the Ipfp function.
  #   seed: The initial multi-dimensional array updated by Ipfp. 
  #   target.list: A list of the target margins used by the Ipfp function. Each
  #                component of the list is an array whose cells indicates
  #                which dimension the corresponding margin relates to.
  #   replace.zeros: If 0-cells are to be found in either the seed or the
  #                  estimate arrays, then their values are replaced with this
  #                  value.
  #
  # Returns: A covariance matrix of the estimated probabilities (last index move
  #          fastest).
    
  n <- sum(seed)  
  seed.prob <- Array2Vector(seed / sum(seed))
  estimate.prob <- Array2Vector(estimate / sum(estimate))
  
  # checking if 0-cells values and replace them with a small value  
  seed.prob <- ifelse(seed.prob == 0, replace.zeros, seed.prob)    
  estimate.prob <- ifelse(estimate.prob == 0, replace.zeros, estimate.prob)
  
  # computation of the diagonal matrix filled with the inverses of seed and
  # estimated probabilities
  D.seed <- diag(1 / seed.prob)
  D.estimate <- diag(1 / estimate.prob)
  
  # computation of A such that A' * vector(estimate) = vector(target.data)
  # ... one line filled with ones
  A.transp <- matrix(1, nrow = 1, ncol = length(estimate.prob))  
  # ... constrainst (removing the first one since it is redundant information)
  for (j in 1:length(target.list)) {
    marg.mat <- cmm::MarginalMatrix(var = 1:length(dim(seed)), 
                                    marg = target.list[[j]], 
                                    dim = dim(seed))[-1,]
    A.transp <- rbind(marg.mat, A.transp, deparse.level = FALSE)
  }  
  A <- t(A.transp)
  
  # removing the linearly dependant columns from A (redundant constrainst)
  A <- GetLinInd(A)$mat.li
    
  # computation of the orthogonal complement of A (using QR decomposition)
  K <- qr.Q(qr(A), complete = TRUE)[, (dim(A)[2] + 1):dim(A)[1]]

  # computation of the variance  
  estimate.var <- (1 / n) * K %*% solve((t(K) %*% D.estimate %*% K)) %*%
                  t(K) %*% D.seed %*% K %*%
                  solve(t(K) %*% D.estimate %*% K) %*% t(K)   

  # returning the result
  return(estimate.var)
  
}
