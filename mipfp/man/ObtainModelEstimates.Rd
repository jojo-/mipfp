\name{ObtainModelEstimates}
\alias{ObtainModelEstimates}
\title{Estimating a contingency table using model-based approaches
}
\description{This function provides several alternative estimating methods to 
the IPFP when estimating a multiway table subject to known constrains/totals: 
maximum likelihood method (ML), minimum chi-squared (CHI2) and weighted least 
squares (WLSQ). Note that the resulting estimators are probabilities.
   
The covariance matrix of the estimated proportions (as defined by Little and Wu, 
1991) are also provided. Also in the case of the ML method, the covariance 
matrix defined by Lang (2004) is also returned.

}
\usage{
ObtainModelEstimates(seed, target.list, target.data, method="ml", 
                     replace.zeros = 1e-10, compute.cov = FALSE, ...)
}
\arguments{
 \item{seed}{
    The initial multi-dimensional array to be updated. Each cell must
    be non-negative.
}
  \item{target.list}{
    A list of the target margins provided in \code{target.data}. Each component
    of the list is an array whose cells indicates which dimension the
    corresponding margin relates to.
}
  \item{target.data}{
    A list containing the data of the target margins. Each
    component of the list is an array storing a margin.
    The list order must follow the one defined in \code{target.list}.
    Note that the cells of the arrays must be non-negative.
}
  \item{method}{Determine the model to be used for estimating the contingency
    table. By default the method is \code{ml} (maximum likelihood); other 
    options available are \code{chi2} (minimum chi-squared) and \code{lsq} 
    (least squares).
}
  \item{replace.zeros}{
    Constant that is added to zero cell found in the seed,
    as procedures require strictly positive cells. Default value is 
    \eqn{1e^{-10}}{1e-10}.
}
  \item{compute.cov}{
    If set to TRUE, then the function also return the (asymptotic) covariance 
    matrices of the updated cells and cells proportion as well as the 
    associated standard errors. Default = FALSE.
}
  \item{...}{Additional parameters that can be passed to control the 
    optimisation process (see \link[Rsolnp:solnp]{solnp} from the package 
    \link[Rsolnp]{Rsolnp}).
}
}
\value{
A list containing the final estimated table as well as the covariance matrix of
the estimated proportion and other convergence informations.
  \item{x.hat}{
    Array of the estimated table frequencies.
}
  \item{p.hat}{
    Array of the estimated table probabilities.
}
  \item{x.hat.cov}{
    Asymptotic covariance matrix for \code{x.hat} using the Delta method (if 
    \code{compute.cov = TRUE}).
}
  \item{x.hat.se}{
    Vector of standard errors for \code{x.hat} using the Delta method (if 
    \code{compute.cov = TRUE}).
}
  \item{p.hat.cov}{
    Asymptotic covariance matrix for \code{p.hat} using the Delta method (if 
    \code{compute.cov = TRUE}).
}
  \item{p.hat.se}{
    Vector of standard errors for \code{p.hat} using the Delta method (if 
    \code{compute.cov = TRUE}).
}
  \item{check.margins}{
    For each list element of \code{target.data}, \code{check.margins} shows the 
    maximum absolute deviation between the element and the corresponding 
    estimated margin. Note that the deviations should approximate zero, 
    otherwise the target margins are not met.
}
  \item{solnp.res}{
    The estimation process uses the \code{solnp} optimisation function from 
    the R package Rsolnp and \code{solnp.res} is the corresponding object 
    returned by the solver.
}

Additionnaly, it the estimation method is \code{ml}, then the function also 
returns the following elements.
  \item{G2}{
    Log-likelihood statistic for testing that the target constraints are met.
}
  \item{W2}{
    Wald statistic for testing that the target constraints are met.
}
  \item{X2}{
    Pearson chi-squared statistic for testing that the target constraints are 
    met.
}
  \item{df}{
    Degrees of freedom relevant for the \code{G2}, \code{W2} and \code{X2} 
    statistics.
}
  \item{lang$x.hat.cov}{
    Asymptotic covariance matrix for \code{x.hat} using Lang's method (if
    \code{compute.cov = TRUE}).
}
  \item{lang$x.hat.se}{
    Vector of standard errors for \code{x.hat} using Lang's method (if 
    \code{compute.cov = TRUE}).
}
  \item{lang$p.hat.cov}{
    Asymptotic covariance matrix for \code{p.hat} using Lang's method (if 
    \code{compute.cov = TRUE}).
}
  \item{lang$p.hat.se}{
    Vector of standard errors for \code{p.hat} using Lang's method (if 
    \code{compute.cov = TRUE}).
}
}
\references{
 Little, R. J., Wu, M. M. (1991)
 Models for contingency tables with known margins when target and sampled 
 populations differ.
 \emph{Journal of the American Statistical Association} 86 (413): 87-95.
 
 Lang, J.B. (2004) 
 Multinomial-Poisson homogeneous models for contingency tables. 
 \emph{Annals of Statistics} 32(1): 340-383.
}
\author{
Thomas Suesse 

Maintainer: Johan Barthelemy <johan@uow.edu.au>.
}
\seealso{
 \code{\link[Rsolnp:solnp]{solnp}} function documentation of the package
 \code{\link[Rsolnp]{Rsolnp}} for the details of the 
 \code{solnp.res} object returned by the function.
}
\examples{
# set-up an initial 3-way table of dimension (2 x 2 x 2)
seed <- Vector2Array(c(80, 60, 20, 20, 40, 35, 35, 30), dim = c(c(2, 2, 2)))

# building target margins
margins12 <- c(2000, 1000, 1500, 1800)
margins12.array <- Vector2Array(margins12, dim=c(2, 2))
margins3 <- c(4000,2300)
margins3.array <- Vector2Array(margins3, dim = 2) 
target.list <- list(c(1, 2), 3)
target.data <- list(margins12.array, margins3.array)

# estimating the new contingency table using the ml method
results.ml <- ObtainModelEstimates(seed, target.list, target.data, 
                                   compute.cov = TRUE)
print(results.ml)

# estimating the new contingency table using the chi2 method
results.chi2 <- ObtainModelEstimates(seed, target.list, target.data, 
                                     method = "chi2", compute.cov = TRUE)
print(results.chi2)

# estimating the new contingency table using the lsq method
results.lsq <- ObtainModelEstimates(seed, target.list, target.data,
                                    method = "lsq", compute.cov = TRUE)
print(results.lsq)
}
\concept{multiway contingency table estimation}
\concept{multivariate contingency table estimation}
\keyword{models}
\keyword{multivariate}
\keyword{array}