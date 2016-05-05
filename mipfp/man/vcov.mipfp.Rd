\name{vcov.mipfp}
\alias{vcov.mipfp}
\title{
Calculate variance-covariance matrix for mipfp objects
}
\description{
This function determines the (asymptotic) covariance matrix of the estimates 
in an \code{mipfp} object using either the Delta formula designed by Little 
and Wu (1991) or Lang's formula (2004).
}
\usage{
\method{vcov}{mipfp}(object, method.cov = "delta", seed = NULL, 
     target.data = NULL, target.list = NULL, replace.zeros = 1e-10, ...)
}
\arguments{
  \item{object}{
    An object of class mipfp.
}
  \item{method.cov}{
    Select the method to use for the computation of the covariance.
    The available methods are \code{delta} and \code{lang}.
}
  \item{seed}{
    The initial multi-dimensional array used to create \code{object}
    (optional). 
}
  \item{target.data}{
    A list containing the data of the target margins used to create 
    \code{object}. Each component of the list is an array storing a margin. 
    The list order must follow the one defined in \code{target.list}
    (optional).
}
  \item{target.list}{
    A list of the target margins used to create \code{object} function. Each
    component of the list is an array whose cells indicates which dimension the
    corresponding margin relates to (optional).
}
  \item{replace.zeros}{
    If 0-cells are to be found, then their values are replaced with this
    value.
}
  \item{...}{
    Not used.  
}
}
\details{
The asymptotic covariance matrix of the estimates probabilities using Delta's 
formula has the form (Little and Wu, 1991)
\deqn{K(K^T D1^{-1} K)^{-1} K^T D2^{-1} K (K^T D1^{-1} K)^{-1} K^T}{
      K * inv(t(K) * inv(D1) * K) * t(K) * inv(D2) * K * inv(t(K) * inv(D1)
      * K) * t(K)}
where
\itemize{
  \item{\code{K} is the orthogonal complement of the marginal matrix, i.e. the 
    matrix \eqn{A}{A} required to obtain the marginal frequencies \eqn{m}{m};
}
  \item{\code{D1} and \code{D2} are two diagonal matrices whose components
    depends on the estimation process used to generate \code{object}.
}
}
If the estimation process has been done using
\itemize{
  \item{\code{ipfp}}{
    then \eqn{diag(D1) = \hat{p}}{diag(D1) = p.hat} and 
    \eqn{diag(D2) = p_*}{diag(D2) = p.seed};
}
  \item{\code{ml}}{
    then \eqn{diag(D1) = \frac{\hat{p}^2}{p_*}}{diag(D1) = p.hat^2 / p.seed} and 
    \eqn{diag(D2) = diag(D1)}{diag(D2) = diag(D1)};
}
  \item{\code{chi2}}{
    then \eqn{diag(D1) = \frac{\hat{p}^4}{p_*^3}}{diag(D1) = p.hat^4 / p.seed^3}
    and \eqn{diag(D2) = diag(D1)}{diag(D2) = diag(D1)};
}
  \item{\code{lsq}}{
    then \eqn{diag(D1) = p_*}{diag(D1) = p.seed} and 
    \eqn{diag(D2) = \frac{p_*3}{\hat{p}^2}}{diag(D2) = p.seed^3 / p.hat^2};
}
}
where \eqn{\hat{p}}{p.hat} is the vector of estimated probabilities and
\eqn{p_*}{p.seed} is the vector of the seed probabilities.

Using Lang's formula (2004), the covariance matrix becomes
\deqn{\frac{1}{N} 
      \left( D - \hat{p}\hat{p}^T - D H(H^T D H)^{-1} H^T D \right)}{%
      1/N (D - p.hat * t(p.hat) - D * H * inv(t(H) * D * H) * t(H) * D) }
where
\itemize{
  \item{\code{D}}{
    is a diagonal matrix of the estimated probabilities \eqn{\hat{p}}{p.hat};
  }
  \item{\code{H}}{
    denotes the Jacobian evaluated in \eqn{\hat{p}}{p.hat} of the function
    \eqn{h(p) = A^T p - m}{h(p) = t(A) * p - m}.
  }
}
}
\value{
A list with the following components:
  \item{x.hat.cov}{
    A covariance matrix of the estimated counts (last index move fastest)
    computed using the method specified in \code{cov.method}.
}
  \item{p.hat.cov}{
    A covariance matrix of the estimated probabilities (last index 
    move fastest) computed using the method specified in \code{cov.method}.
}
  \item{x.hat.se}{
    The standard deviation of the estimated counts (last index move fastest)
    computed using the method specified in \code{cov.method}.
}
  \item{p.hat.se}{
    The standard deviation of the estimated probabilities (last index move 
    fastest) computed using the method specified in \code{cov.method}.
}
  \item{df}{
    Degrees of freedom of the estimates.
}
  \item{method.cov}{
    The method used to compute the covariance matrix.
}
}
\references{
 Lang, J.B. (2004) 
 Multinomial-Poisson homogeneous models for contingency tables. 
 \emph{Annals of Statistics} 32(1): 340-383.

 Little, R. J., Wu, M. M. (1991)
 Models for contingency tables with known margins when target and seed 
 populations differ.
 \emph{Journal of the American Statistical Association} 86 (413): 87-95.
}
\author{
Johan Barthelemy.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\seealso{
\code{\link[mipfp:Estimate]{Estimate}} function to create an object of class
\code{mipfp} and to update an initial multidimensional array with respect to 
given constraints.
}
\examples{
# true contingency (2-way) table
true.table <- array(c(43, 44, 9, 4), dim = c(2, 2))
# generation of sample, i.e. the seed to be updated
seed <- ceiling(true.table / 10)
# desired targets (margins)
target.row <- apply(true.table, 2, sum)
target.col <- apply(true.table, 1, sum)
# storing the margins in a list
target.data <- list(target.col, target.row)
# list of dimensions of each marginal constrain
target.list <- list(1, 2)
# calling the Estimate function
res <- Estimate(seed, target.list, target.data)
# printing the variance-covariance matrix
print(vcov(res))
}
\concept{multiway contingency table estimation}
\concept{multivariate contingency table estimation}
\concept{variance-covariance matrix}
\concept{covariance}
\keyword{models}
\keyword{multivariate}
\keyword{array}