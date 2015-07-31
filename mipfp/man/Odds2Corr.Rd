\name{Odds2Corr}
\alias{Odds2Corr}
\title{Converting odds ratio to correlation
}
\description{

For \eqn{K} binary (Bernoulli) random variables 
\eqn{X_1}{X_1}, ..., \eqn{X_K}{X_K}, this function transforms the odds ratios 
measure of association \eqn{O_{ij}}{O_ij} between every pair
\eqn{(X_i, X_j)}{(X_i, X_j)} to the correlation \eqn{C_{ij}}{C_ij} where
\deqn{C_{ij} = \frac{cov(X_i, X_j)}{\sqrt{var(X_i) * var(X_j)}}}{
C_ij = cov(X_i, X_j) / sqrt(var(X_i) * var(X_j))}
and 
\deqn{
O_{ij} = \frac{P(X_i = 1, X_j = 1) * P(X_i = 0, X_j = 0)}
            {P(X_i = 1, X_j = 0) * P(X_i = 0, X_j = 1)}.
}{
O_ij = P(X_i = 1, X_j = 1) * P(X_i = 0, X_j = 0) / 
       P(X_i = 1, X_j = 0) * P(X_i = 0, X_j = 1).}
}
\usage{
Odds2Corr(odds, marg.probs)
}
\arguments{
  \item{odds}{      
    A \eqn{K \times K}{K x K} matrix where the \eqn{i}-th row and the \eqn{j}-th
    column represents the odds ratio \eqn{O_{ij}}{O_ij} between variables
    \eqn{i} and \eqn{j}.  
}
  \item{marg.probs}{
    A vector with \eqn{K} elements of marginal probabilities where the 
    \eqn{i}-th entry refers to \eqn{P(X_i = 1)}{P(X_i = 1)}.
}
}
\value{
The function return a list with the correlations and the pairwise probabilities.
  \item{corr}{
    A matrix of the same dimension as \code{odds} containing the correlations
}
  \item{pair.proba}{
    A matrix of the same dimension as \code{odds} containing the pairwise 
    probabilities.
}
}
\references{
  Lee, A.J. (1993).
  Generating Random Binary Deviates Having Fixed Marginal Distributions and 
  Specified Degrees of Association
  \emph{The American Statistician} 47 (3): 209-215.
 
  Qaqish, B. F., Zink, R. C., and Preisser, J. S. (2012).
  Orthogonalized residuals for estimation of marginally specified association
  parameters in multivariate binary data.
  \emph{Scandinavian Journal of Statistics} 39, 515-527.
}
\author{
Thomas Suesse.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\seealso{
 \code{\link[mipfp]{Corr2Odds}} for converting correlation to odds ratio.
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

# converting odds ratio to correlation
corr <- Odds2Corr(odds = or, marg.probs = p)
print(corr)
}
\concept{correlation}
\concept{Bernoulli}
\concept{odds ratio}
\keyword{multivariate}
\keyword{array}