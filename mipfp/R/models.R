# File mipfp/R/models.R
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
# This file provides a function to estimate a contingency table using several
# model-based approaches.
# ------------------------------------------------------------------------------

library(Rsolnp)
library(cmm)
library(numDeriv)

ObtainModelEstimates <- function(seed, target.list, target.data, method = "ml", 
                                 replace.zeros = 1e-10, compute.cov = FALSE,
                                 ...) {
  # Estimates N-way tables using max likelihood, min chi2 or least squares.
  # 
  # This function provides several alternative estimating methods to the IPFP 
  # when estimating multiway table subject to known constrains/totals: Maximum 
  # likelihood method (ML), minimum chi-squared (CHI2) and weighted least 
  # squares (LSQ). Note that the resulting estimators are probabilities.
  # 
  # The covariance matrix of the estimators (as defined in the paper "Models for 
  # contingency tables with known margins when target and sampled populations 
  # differ" by Little and Wu (1991) are also provided. Also in the case of the
  # ML method, the covariance matrix defined by Lang in "Multinomial-Poisson 
  # homogeneous models for contingency tables" is also returned.
  #
  # Author: T. Suesse
  #  
  # Args:
  #   seed: The initial multi-dimensional array to be updated. Each cell must 
  #         be non-negative.
  #   target.list: A list of the target margins provided in target.data. Each 
  #                component of the list is an array whose cells indicates which
  #                dimension the corresponding margin relates to.
  #   target.data: A list containing the data of the target margins. Each 
  #                component of the list is an array storing a margin. The list 
  #                order must follow the one defined in target.list. Note that 
  #                the cells of the arrays must be non-negative.
  #   method: by default "ml" (Maximum likelihood), other options "chi2" (
  #           minimum chi-squared) and "lsq" (least squares).
  #   replace.zeros: constant that is added to zero cell counts, as the 
  #                  procedures require strictly positive cell counts.
  #   ...: Additional parameters that can be passed to control the optimisation
  #        process.
  #     
  # Returns: A list whose elements are  
  #   pi.hat: Array of the estimated table probabilities.
  #   xi.hat: Array of the estimated table frequencies.  
  #   pi.se: Vector of estimates' standard errors for pi.hat using Delta method.
  #   pi.cov: Asymptotic covariance matrix for pi.hat using Delta method.
  #   check.margins: for each list element of target.data, check.margins shows 
  #                  the maximum absolute value of A * pi.hat - margins.vector.
  #                  The elements should approximate zero, otherwise margins are
  #                  not met.
  #   solnp.res: For optimisation it uses the R package Rsolnp and solnp is the 
  #              corresponding object returned by Rsolnp.
  #   lang$pi.cov: Asymptotic covariance matrix for pi.hat using Lang's method.
  #   lang$pi.se: Vector of estimates' standard errors for pi.hat using Lang's
  #               estimator.  
  #   lang$G2: Log-likelihood statistic for testing that the target constraints
  #            are met.
  #   lang$W2: Wald statistic for testing that the target constraints are met.
  #   lang$X2: Person chi-squared statistic for testing that the target 
  #            constraints are met.
  #   lang$df: Degrees of freedom for G2, W2 and X2 statistics.  
  
  # checking if a seed is provided
  if (is.null(seed) == TRUE) {
    stop('Error: no seed specified!')
  }
  
  # checking if target are provided
  if (is.null(target.data) == TRUE | is.null(target.list) == TRUE) {
    stop('Error: target.data and/or target.list not specified!')
  }
  
  # checking non negativity condition for the seed and the target
  if (min(sapply(target.data, min)) < 0.0 | min(seed) < 0.0) {
    stop('Error: Target and Seed cells must be non-negative!')    
  }
  
  # checking if NA in target cells
  if (is.na(min(sapply(target.data, min))) == TRUE)  {
    stop('Error: NA values present in the margins')
  }
  
  # checking if NA in seed
  if (is.na(min(seed)) == TRUE) {
    stop('Error: NA values present in the seed!')
  }
  
  # checking the strict positivy of replace.zeros
  if ( replace.zeros <= 0.0 ) {
    stop('Error: replace.zeros must be strictly positive!')
  }
  
  # checking the margins consistency if no missing values in the targets
  check.margins <- TRUE
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

  # if margins are not consistent, shifting from frequencies to probabilities
  if (check.margins == FALSE) {
    seed <- seed / sum(seed)
    for (m in 1:length(target.data)) {
      target.data[[m]] <- target.data[[m]] / sum(target.data[[m]])
    }
  }  
  
  # create vector version the seed compatible with the function 'MarginalMatrix' 
  # from cmm that requires that the last index moves fatest
  seed.vector <- Array2Vector(seed)
    
  # dimensions of the problem
  n.sets.margins <- length(target.list)    
  K <- dim(seed)
  K.length <- length(K)
  K.prod <- prod(K)
  n <- sum(seed.vector)
  target.n = sum(target.data[[1]])
  
  # scaling input to probabilities and removing 0 cells
  seed.vector[seed.vector == 0] <- replace.zeros * n
  seed.prop.vector <- seed.vector / n  
        
  # generation of the constraints matrix A such that A * pi.hat' = target.data',
  # where pi.hat' and target.data' are obtained by removing the first element of
  # each target dimension in order for A to be full rank
  length.margins <- rep(0, K.length)
  margins.vector <- NULL
  A <- matrix(0, nrow = 0, ncol = K.prod)
  
  # ... generating the marginal matrix row by row
  for (k in 1:n.sets.margins) {
    
    # convert k-th margin to an array
    temp.margins <- Array2Vector(target.data[[k]])
    temp.margins <- temp.margins[!is.na(temp.margins)]
    
    # remove first condition as it is redundant information
    temp.margins <- temp.margins[-1] / sum(temp.margins)
    margins.vector <- c(margins.vector, temp.margins)
    
    # construct the current row of the marginal matrix
    A.temp <- cmm::MarginalMatrix(1:K.length, list(target.list[[k]]), K)
    
    # remove first entry as it is always redundant
    A <- rbind(A, A.temp[-1, ]) 
    length.margins[k] <- length(temp.margins)
    
  }
  
  # removing the linearly dependant rows
  temp <- GetLinInd(t(A))
  A <- t(temp$mat.li)    
  dim.A <- dim(A)      
  
  # and updating the vector of target margins accordingly
  margins.vector <- margins.vector[temp$idx]  
  
  # adding a final row of 1 to A to generate A.new
  A.new <- rbind(A, 1) 
    
  # defining some functions used by the optimisation
  # ... log-likelihood function
  fun.loglik <- function(p) {
    return(- sum(seed.vector * log(p)))
  }
  
  # ... chi-square function
  fun.chisq <- function(p) {
    return(sum((p - seed.prop.vector)^2 / p))
  }
  
  # ... least-square function
  fun.lsq <- function(p) {
    return(sum((p - seed.prop.vector)^2 / seed.prop.vector))
  }
  
  # ... equality constraints
  eqfun1 <- function(p) {    
    return(A.new %*% p - c(margins.vector, 1))
  }
  
  # switching to the desired user method
  # ... converting method to a integer index
  method.num <- switch(method, ml = 1, chi2 = 2, lsq = 3, 4)
  if (method.num == 4) {
    warning("'method' must be 'ml', 'chi2' or 'lsq', switching to 'ml'!")
    method.num <- 1
  }
  # ... index determine the method: 1 = ML, 2 = CHI2 and 3 = LSQ
  switch(method.num,
    fun <- fun.loglik,
    fun <- fun.chisq,
    fun <- fun.lsq
  )
  
  # calling of "solnp" of package "Rsolnp"
  rsolnp <- try(Rsolnp::solnp(pars=seed.prop.vector, fun = fun,
                              LB = rep(0 + replace.zeros^1.5, K.prod), 
                              UB = rep(1 - replace.zeros^1.5, K.prod), 
                              eqfun = eqfun1, eqB = c(rep(0, dim(A.new)[1])),                              
                              control = list(trace = 0, ...)))
  
  # assessing solnp's convergence, returning an error if no convergence
  if (inherits(rsolnp, "try-error") | rsolnp$convergence > 0 ){
    warning("No reliable solutions found by solnp!
            Check the delta and tol parameters.
            replace.zeros might be too low!\n")   
  }
      
  # saving solution: probabilities and frequencies
  pi.hat <- rsolnp$pars
  pi.hat.array <- Vector2Array(pi.hat, dim.out = K)  
  xi.hat.array <- pi.hat.array * target.n;    
    
  # computing final max difference between generated and target margins
  check.margins <- vector(mode = "numeric", length = n.sets.margins)  
  for (j in 1:n.sets.margins) {
    check.margins[j] = max(abs(target.data[[j]] 
                          - apply(xi.hat.array, target.list[[j]], sum))) 
  }      
  
  # gathering the results
  results <- list("x.hat" = xi.hat.array, "p.hat" = pi.hat.array,                  
                  "check.margins" = check.margins, "solnp.res" = rsolnp)  
  
  
  # compute variance estimators using the Delta method (Little and Wu, 1991)
  if (compute.cov == TRUE) {
    # ... computing D1 and D2 matrices according to the chosen method
    switch(method.num, { 
      # ML 
      D1.inv <- diag( 1 / (pi.hat^2 / seed.prop.vector))
      D2.inv <- D1.inv
    }, {
      # CHI2
      D1.inv <- diag(1 / (pi.hat^4 / seed.prop.vector^3))
      D2.inv <- diag(1 / (pi.hat^4 / seed.prop.vector^3))
    }, {
      # LSQ
      D1.inv <- diag(1 / seed.prop.vector)
      D2.inv <- diag(1 / (seed.prop.vector^3 / pi.hat^2))
    })  
    
    # ... obtain orthogonal complement of A.new using QR decomposition
    A.new.t <- t(A.new)     
    U <- qr.Q(qr(A.new.t), 
              complete = TRUE)[,(dim(A.new.t)[2] + 1):dim(A.new.t)[1]]
    
    if (is.null(dim(U)) == TRUE) {
      U <- t(U)
    } 
    
    # ... computing the variance
    pi.VCov <- (1 / n) * U %*% solve(t(U) %*% D1.inv %*% U) %*% 
               t(U) %*% D2.inv %*% U %*% solve(t(U) %*% D1.inv %*% U) %*% t(U)
    xi.VCov <- pi.VCov * sum(xi.hat.array)^2
    
    # estimates' standart error
    pi.se <- sqrt(diag(pi.VCov)) 
    xi.se <- sqrt(diag(xi.VCov))
    
    # updating the results
    results$p.hat.cov <- pi.VCov
    results$x.hat.cov <- xi.VCov
    results$p.hat.se <- pi.se
    results$x.hat.se <- xi.se    
    
  }
  
  # calculate Lang's covariance matrix and various statistic test if method = ML
  if (method.num == 1) {          
  
    # constraint function h(pi) = A * pi - margins
    h.fct <- function(m) {            
      return(A %*% (m / sum(m)) - margins.vector)
    }
    
    # compute statistics for testing if constraints are met and appending to
    # the results
    H.seed <- t(numDeriv::jacobian(h.fct, seed.vector))
    h.Y <- h.fct(seed.vector)
    # ... log-likelihood ration
    results$G2 <- 2 * sum(seed.vector * log(seed.prop.vector / pi.hat))
    # ... Wald test statistic
    results$W2 <- as.double(t(h.Y) %*% solve(t(H.seed) %*% 
                            diag(seed.vector) %*% H.seed) %*% h.Y)
    # ... Pearson Chi-square
    results$X2 <- as.double(t(seed.vector - n * pi.hat) %*%
                              diag(1 / (n * pi.hat)) %*%
                              (seed.vector - n * pi.hat))
    # ... associated degree of freedoms (dim of constraint function h)
    results$df <- dim.A[1]
    
    # Lang's covariance if requested
    if (compute.cov == TRUE) {
      H.pi <- t(numDeriv::jacobian(h.fct, pi.hat))
      D.pi <- diag(pi.hat)
      pi.VCov.Lang <- 1 / n * (D.pi - pi.hat %*% t(pi.hat) - D.pi %*% H.pi %*% 
                               solve(t(H.pi) %*% D.pi %*% H.pi) %*% t(H.pi)
                               %*% D.pi)
      xi.VCov.Lang <- pi.VCov.Lang * sum(xi.hat.array)^2
      
      # extracting standart errors
      pi.se.Lang <- sqrt(diag(pi.VCov.Lang))
      xi.se.Lang <- sqrt(diag(xi.VCov.Lang))
      
      # ... appending the Lang's covariance estimation to the results      
      results$lang <- list("p.hat.se" = pi.se.Lang, "p.hat.cov" = pi.VCov.Lang,
                           "x.hat.se" = xi.se.Lang, "x.hat.cov" = xi.VCov.Lang)
      
    }    
               
  }  

  # returning results
  return(results)

}
