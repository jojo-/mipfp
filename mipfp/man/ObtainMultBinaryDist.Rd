\name{ObtainMultBinaryDist}
\alias{ObtainMultBinaryDist}
\title{Generating a multivariate Bernoulli joint-distribution
}
\description{
This function applies the IPFP procedure to obtain a joint distribution of 
\eqn{K} multivariate binary (Bernoulli) variables \eqn{X_1}, ..., \eqn{X_K}.

It requires as input the odds ratio or alternatively the correlation as a 
measure of association between all the binary variables and a vector of marginal 
probabilities. 

This function is useful when one wants to simulate and draw from a 
multivariate binary distribution when only first order (marginal probabilities)
and second order moments (correlation or odds ratio) are available.
}
\usage{
ObtainMultBinaryDist(odds = NULL, corr = NULL, marg.probs, ...)
}
\arguments{
  \item{odds}{  
   A \eqn{K \times K}{K x K} matrix where the \eqn{i}-th row and the \eqn{j}-th 
   column represents the Odds ratio between variables \eqn{i} and \eqn{j}. Must 
   be provided if \code{corr} is not.
}
  \item{corr}{
    A \eqn{K \times K}{K x K} matrix where the \eqn{i}-th row and the \eqn{j}-th
    column represents the correlation between variables \eqn{i} and \eqn{j}. 
    Must be provided if \code{odds} is not.
}
  \item{marg.probs}{
    A vector with \eqn{K} elements of marginal probabilities where the 
    \eqn{i}-th entry refers to \eqn{P(X_i = 1)}{P(X_i = 1)}.
}
  \item{...}{
    Additional arguments that can be passed to the \code{Ipfp} function such as 
    \code{tol}, \code{iter}, \code{print} and \code{compute.cov}.
}

}

\value{
A list whose elements are mainly determined by the \code{Ipfp} function.
  \item{joint.proba}{
    The resulting multivariate joint-probabilities (from \code{Ipfp}).
}
  \item{stp.crit}{
    The final value of the \code{Ipfp} stopping criterion.
}
  \item{conv}{
    Boolean indicating whether the \code{Ipfp} algorithm converged to a 
    solution.
}
  \item{check.margins}{
    A list returning, for each margin, the absolute maximum deviation between 
    the desired and generated margin. Ideally the elements should approximate
    0 (from \code{Ipfp}).
}
  \item{label}{
    The names of the variables.
}
}
\references{
  Lee, A.J. (1993).
  Generating Random Binary Deviates Having Fixed Marginal Distributions and 
  Specified Degrees of Association
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
\note{
It is important to note that either the odds ratio defined in \code{odds} or
the correlations described in \code{corr} must be provided.
}
\seealso{
\code{\link[mipfp:Ipfp]{Ipfp}} for the function used to estimate the 
distribution; \code{\link[mipfp:RMultBinary]{RMultBinary}} to simulate the
estimated joint-distribution; \code{\link[mipfp:Corr2Odds]{Corr2Odds}} and
\code{\link[mipfp:Odds2Corr]{Odds2Corr}} to convert odds ratio to correlation
and conversely.
}
\examples{
# initial odds ratios from Qaqish et al. (2012)
or <- matrix(c(Inf, 0.281, 2.214, 2.214, 
               0.281, Inf, 2.214, 2.214,
               2.214, 2.214, Inf, 2.185,
               2.214, 2.214, 2.185, Inf), nrow = 4, ncol = 4, byrow = TRUE)
rownames(or) <- colnames(or) <- c("Parent1", "Parent2", "Sibling1", "Sibling2")

# hypothetical marginal probabilities
p <- c(0.2, 0.4, 0.6, 0.8)

# estimating the joint-distribution
p.joint <- ObtainMultBinaryDist(odds = or, corr = NULL, marg.probs = p)
print(p.joint$joint.proba)

# obtain identical solution when providing correlation
corr <- Odds2Corr(odds = or, marg.probs = p)$corr
p.joint.alt <- ObtainMultBinaryDist(corr = corr, marg.probs = p)

# checking if the results are truly identicals
diff <- sum(abs(p.joint.alt$joint.proba - p.joint$joint.proba))
cat('Sum of the absolute deviations: ', diff, '\n')
}
\concept{Bernoulli}
\keyword{distribution}
\keyword{multivariate}
