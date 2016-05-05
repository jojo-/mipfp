\name{Corr2Odds}
\alias{Corr2Odds}
\title{Converting correlation to odds ratio
}
\description{
For \eqn{K} binary (Bernoulli) random variables 
\eqn{X_1}{X_1}, ..., \eqn{X_K}{X_K}, this function transforms the correlation
measure of association \eqn{C_{ij}}{C_ij} between every pair
\eqn{(X_i, X_j)}{(X_i, X_j)} to the odds ratio \eqn{O_{ij}}{O_ij} where
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
Corr2Odds(corr, marg.probs)
}
\arguments{
  \item{corr}{      
    A \eqn{K \times K}{K x K} matrix where the \eqn{i}-th row and the \eqn{j}-th
    column represents the correlation \eqn{C_{ij}}{C_ij} between variables
    \eqn{i} and \eqn{j}.  
}
  \item{marg.probs}{
    A vector with \eqn{K} elements of marginal probabilities where the 
    \eqn{i}-th entry refers to \eqn{P(X_i = 1)}{P(X_i = 1)}.
}
}
\value{
The function return a list with the correlations and the pairwise probabilities.
  \item{odds}{
    A matrix of the same dimension as \code{corr} containing the correlations
}
  \item{pair.proba}{
    A matrix of the same dimension as \code{corr} containing the pairwise 
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
\code{\link[mipfp:Corr2Odds]{Corr2Odds}} for converting correlation to odds 
ratio.
}
\examples{
# correlation matrix from Qaqish et al. (2012)
cr <- matrix(c( 1.000, -0.215, 0.144, 0.107,
               -0.215,  1.000, 0.184, 0.144,
                0.144,  0.184, 1.000, 0.156,
                0.107,  0.144, 0.156, 1.000), nrow = 4, ncol = 4, byrow = TRUE)
rownames(cr) <- colnames(cr) <- c("Parent1", "Parent2", "Sibling1", "Sibling2")

# hypothetical marginal probabilities
p <- c(0.2, 0.4, 0.6, 0.8)

# converting correlation to odds ratio and getting pairwise probabilities
or <- Corr2Odds(corr = cr, marg.probs = p)
print(or)
}
\concept{correlation}
\concept{odds ratio}
\concept{Bernoulli}
\concept{odds ratio}
\keyword{multivariate}
\keyword{array}