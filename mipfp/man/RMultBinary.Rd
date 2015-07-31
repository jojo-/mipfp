\name{RMultBinary}
\alias{RMultBinary}
\alias{Bernoulli}
\alias{RMultBernoulli}
\title{Simulating a multivariate Bernoulli distribution
}
\description{
This function generates a sample from a multinomial distribution of \eqn{K}
dependent binary (Bernoulli) variables
\eqn{(X_1, X_2, ..., X_K)}{(X_1, X_2, ..., X_K)} defined by an array 
(of 2^K cells) detailing the joint-probabilities.
}
\usage{
RMultBinary(n = 1, mult.bin.dist, target.values = NULL)
}
\arguments{
  \item{n}{
    Desired sample size. Default = 1.
}
  \item{mult.bin.dist}{
    A list describing the multivariate binary distribution. It can be generated
    by the \code{\link[mipfp]{ObtainMultBinaryDist}} function.
    The list contains at least the element \code{joint.proba}, an array 
    detailing the joint-probabilities of the \eqn{K} 
    binary variables. The array has \eqn{K} dimensions of size 2, referring to 
    the 2 possible outcomes of the considered variable. Hence, the total number 
    of elements is \eqn{2^K}{2^K}.
    Additionnaly the list can also provides the element \code{var.label}, a list
    containing the names of the \eqn{K} variables.    
}
  \item{target.values}{
    A list describing the possibles outcomes of each binary variable, for 
    instance \{1, 2\}. Default = \{0, 1\}.
}
}
\value{
A list whose elements are detailed herehunder.
  \item{binary.sequences}{
  The generated \eqn{K \times n}{K x n} random sequence.
}
  \item{possible.binary.sequences}{
    The possible binary sequences, i.e. the domain.
}
  \item{chosen.random.index}{
    The index of the random draws in the domain.
}
}
\references{
  Lee, A.J. (1993).
  Generating Random Binary Deviates Having Fixed Marginal Distributions and 
  Specified Degrees of Association.
  \emph{The American Statistician}  47 (3): 209-215.
 
  Qaqish, B. F., Zink, R. C., and Preisser, J. S. (2012).
  Orthogonalized residuals for estimation of marginally specified association
  parameters in multivariate binary data.
  \emph{Scandinavian Journal of Statistics} 39, 515-527.
}
\author{
Thomas Suesse 

Maintainer: Johan Barthelemy <johan@uow.edu.au>.
}
\seealso{
  \code{\link[mipfp]{ObtainMultBinaryDist}} for estimating the 
  joint-distribution required by this function.
}
\examples{
# from Qaqish et al. (2012)
or <- matrix(c(Inf, 0.281, 2.214, 2.214,
               0.281, Inf, 2.214, 2.214,
               2.214, 2.214, Inf, 2.185,
               2.214, 2.214, 2.185, Inf), nrow = 4, ncol = 4, byrow = TRUE)
rownames(or) <- colnames(or) <- c("Parent1", "Parent2", "Sibling1", "Sibling2")

# hypothetical marginal probabilities
p <- c(0.2, 0.4, 0.6, 0.8)

# estimating the joint-distribution
p.joint <- ObtainMultBinaryDist(odds = or, marg.probs = p)

# simulating 100,000 draws from the obtained joint-distribution
y.sim <- RMultBinary(n = 1e5, mult.bin.dist = p.joint)$binary.sequences

# checking results
cat('dim y.sim =', dim(y.sim)[1], 'x', dim(y.sim)[2], '\n')
cat('Estimated marginal probs from simulated data\n')
apply(y.sim,2,mean)
cat('True probabilities\n')
print(p)
cat('Estimated correlation from simulated data\n')
cor(y.sim)
cat('True correlation\n')
Odds2Corr(or,p)$corr

# generating binary outcomes with outcome different than 0, 1
RMultBinary(n = 10, mult.bin.dist = p.joint, 
            target.values = list(c("A", "B"), c(0, 1), c(1, 2), c(100, 101)))
}
\concept{Bernoulli}
\concept{random draws}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}