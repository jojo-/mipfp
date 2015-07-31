\name{IpfpCov}
\alias{IpfpCov}
\title{
Covariance matrix of the estimators produced by Ipfp
}
\description{
This function determines the (asymptotic) covariance matrix of the estimates 
produced by the iterative proportional fitting procedure using the formula 
designed by Little and Wu (1991).
}
\usage{
IpfpCov(estimate, seed, target.list, replace.zeros = 1e-10)
}
\arguments{
  \item{estimate}{
    The array of estimates produced by the \code{\link[mipfp:Ipfp]{Ipfp}}
    function.
}
  \item{seed}{
    The intial array (seed) that was updated by the 
    \code{\link[mipfp:Ipfp]{Ipfp}} function.
}
  \item{target.list}{
    A list of dimensions of the marginal target constrains. Each component of 
    the list is an array whose cells indicate which dimension the corresponding 
    margin relates to.
}
  \item{replace.zeros}{
    If a cell of the \code{estimate} or the \code{seed} has a value equals to 
    0, then it is replaced with this value. Default is 1e-10.
}
}
\details{
 The asymptotic covariance matrix of the estimates produced by the iterative
 proportional fitting procedure has the form (Little and Wu, 1991)
 \deqn{K(K^T D1^{-1} K)^{-1} K^T D2^{-1} K (K^T D1^{-1} K)^{-1} K^T}{
 K * inv(t(K) * inv(D1) * K) * t(K) * inv(D2) * K * inv(t(K) * inv(D1)
 * K) * t(K)}
 where
 \itemize{
  \item{\code{K} is the orthogonal complement of the marginal matrix, i.e. the 
    matrix required to obtain the marginal frequencies;}
  \item{\code{D1} is a diagonal matrix of the estimates probabilities}
  \item{\code{D2} is a diagonal matrix of the seed probabilities}
 }
}
\value{
  A matrix of dimension \code{length(estimate) x length(estimate)}
  of the asymptotic variance of the estimates produced by \code{Ipfp}.
}
\references{
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
\code{\link[mipfp:Ipfp]{Ipfp}} function to update an initial multidimensional
array with respect to given constraints.
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
# calling the Ipfp function
res <- Ipfp(seed, target.list, target.data)

# computation of the covariance matrix of the produced estimated probabilities
res.cov <- IpfpCov(res$x.hat, seed, target.list)

# 0.95 level confidence interval of the estimates
n <- sum(res$x.hat)
# ... lower bound
ci.lb <- Array2Vector(res$x.hat) - 1.96 * sqrt(n * diag(res.cov))
# ... upperbound
ci.ub <- Array2Vector(res$x.hat) + 1.96 * sqrt(n * diag(res.cov))
}
\concept{multiway contingency table estimation}
\concept{multivariate contingency table estimation}
\keyword{models}
\keyword{multivariate}
\keyword{array}