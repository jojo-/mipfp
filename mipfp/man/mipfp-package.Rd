\name{mipfp-package}
\alias{mipfp-package}
\alias{mipfp}
\docType{package}
\title{
Multidimensional Iterative Proportional Fitting and Alternative Models
}
\description{
An implementation of several methods for updating an initial \eqn{N}-dimensional
array (called a seed) with respect to given target marginal distributions. 
Those targets can also be multi-dimensional. The procedures are also able to 
estimate a (multi-dimensional) contingency table (encoded as an array) 
matching a given set of (multi-dimensional) margins. In that case, each cell 
of the seed must simply be set to 1.

The package provides the iterative proportional fitting procedure (IPFP), also 
known as the RAS algorithm in economics and matrix raking or matrix scaling in 
computer science. Additionnaly several alternative estimating methods to 
the IPFP are also included, namely the maximum likelihood (ML), minimum 
chi-squared (CHI2) and weighted least squares (WLSQ) model-based approaches.

The package also includes an application of the IPFP to simulate and estimate 
the parameters of multivariate Bernoulli distributions.

Finally a function extracting the linearly independant columns from a matrix,
hence returning a matrix of full rank is provided.
}
\details{
\tabular{ll}{
Package: \tab mipfp\cr
Type: \tab Package\cr
Version: \tab 2.2.1\cr
Date: \tab 2015-09-18\cr
Depends: \tab cmm, numDeriv, Rsolnp, R(>= 2.10.0)\cr
License: \tab GPL-2\cr
}
This package provides an implementation of several fitting procedures for 
updating a \eqn{N}-dimensional array with respect to given target marginal 
distributions. Those targets can also multi-dimensional. The available methods
are listed herehunder.

\itemize{
  \item The function \link[mipfp]{Ipfp} provides the iterative proportionnal 
    fitting Procedure. The variance-covariance matrix for the updated cells 
    probabilities is determined by the Delta method (Little and Wu, 1991) can
    alsobe obtained using the function \link[mipfp]{IpfpCov}.
  \item Maximum likelihood, minimum Chi-square and weighted least squares
    approaches are availables in the function
    \link[mipfp]{ObtainModelEstimates}. Note that this function provides the 
    updated table as well as the variance-covariance matrix of the estimated 
    cell probabilities using either the Lang's (2004) or the Delta method 
    (Little and Wu, 1991).
}

The package also includes an application of the IPFP to simulate and estimate 
the parameters of multivariate Bernoulli distributions, respectively in the
functions \link[mipfp]{RMultBinary} and \link[mipfp]{ObtainMultBinaryDist}.
In addition, the functions \link[mipfp]{Corr2Odds}, \link[mipfp]{Odds2Corr}, 
\link[mipfp]{Corr2PairProbs}, \link[mipfp]{Odds2PairProbs} are in turn 
responsible for converting correlation to odds ratio, odds ratio to correlation, 
correlation to pairwise probability and odds ratio to pairwise probability.

Finally, auxillary functions are also provided. \link[mipfp]{Array2Vector} and 
\link[mipfp]{Vector2Array} transforms an array to a vector and vice-versa.
The function \link[mipfp]{GetLinInd} extracting the linearly independant columns 
from a matrix (using QR decomposition) and returning a matrix of full rank is 
also provided.
}
\author{
Johan Barthelemy and Thomas Suesse with contributions from Mohammad Namazi-Rad.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\references{
 Bacharach, M. (1965).
 Estimating Nonnegative Matrices from Marginal Data.
 \emph{International Economic Review} (Blackwell Publishing) 6 (3): 294-310.
 
 Bishop, Y. M. M., Fienberg, S. E., Holland, P. W. (1975).
 \emph{Discrete Multivariate Analysis: Theory and Practice.}
 MIT Press. ISBN 978-0-262-02113-5.
 
 Deming, W. E., Stephan, F. F. (1940).
 On a Least Squares Adjustment of a Sampled Frequency Table When the Expected
 Marginal Totals are Known.
 \emph{Annals of Mathematical Statistics} 11 (4): 427-444.
 
 Fienberg, S. E. (1970).
 An Iterative Procedure for Estimation in Contingency Tables.
 \emph{Annals of Mathematical Statistics} 41 (3): 907-917.

 Golub, G. H., Van Loan C. F. (2012)
 \emph{Matrix Computations. Third Edition.}
 Johns Hopkins University Press.

 Lang, J.B. (2004) 
 Multinomial-Poisson homogeneous models for contingency tables. 
 \emph{Annals of Statistics} 32(1): 340-383.

 Lee, A.J. (1993).
 Generating Random Binary Deviates Having Fixed Marginal Distributions and 
 Specified Degrees of Association
 \emph{The American Statistician}  47 (3): 209-215.
 
 Little, R. J., Wu, M. M. (1991)
 Models for contingency tables with known margins when target and sampled 
 populations differ.
 \emph{Journal of the American Statistical Association} 86 (413): 87-95.
 
 Qaqish, B. F., Zink, R. C., and Preisser, J. S. (2012). 
 Orthogonalized residuals for estimation of marginally specified association
 parameters in multivariate binary data.
 \emph{Scandinavian Journal of Statistics} 39, 515-527.
 
 Stephan, F. F. (1942).
 Iterative method of adjusting frequency tables when expected margins are known.
 \emph{Annals of Mathematical Statistics} 13 (2): 166-178.
}
\seealso{
 \code{ipfp} for a package implementing the ipfp to solve problems of the form
 \eqn{Ax = b}{Ax = b}.
}
\examples{
# generation of an intial 2-ways table to be updated
seed <- array(1, dim=c(2, 2))
# desired targets (margins)
target.row <- c(87, 13)
target.col <- c(52, 48)
# storing the margins in a list
target.data <- list(target.col, target.row)
# list of dimensions of each marginal constrain
target.list <- list(1, 2)
# calling the fitting methods
r.ipfp <- Ipfp(seed, target.list, target.data)
r.ml <- ObtainModelEstimates(seed, target.list, target.data, method = "ml")
r.chi2 <- ObtainModelEstimates(seed, target.list, target.data, method = "chi2")
r.lsq <- ObtainModelEstimates(seed, target.list, target.data, method = "lsq")
}
\concept{iterative proportional fitting procedure}
\concept{ipfp}
\concept{multiway contingency table update}
\concept{multivariate contingency table update}
\concept{RAS algorithm}
\concept{matrix raking}
\concept{matrix scaling}
\concept{Bernoulli}
\keyword{models}
\keyword{multivariate}
\keyword{array}
\keyword{distribution}
\keyword{package}
