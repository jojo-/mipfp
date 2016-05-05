\name{gof.estimates}
\alias{gof.estimates}
\alias{gof.estimates.default}
\alias{gof.estimates.mipfp}
\title{
Wald, Log-likelihood ratio and Person Chi-square statistics for mipfp object
}
\description{
This method computes three statistics to perform a test wheter the seed
agrees with the target data. The statistics are the Wilk's log-likelihood ratio
statistic, the Wald statistic and the Person Chi-square statistic.

The method also returns the associated degrees of freedom.
}
\usage{
\method{gof.estimates}{mipfp}(object, seed = NULL, target.data = NULL, 
              target.list = NULL, replace.zeros = 1e-10, ...)
}
\arguments{
  \item{object}{
    The object of class \code{mipfp} containing.
}
  \item{seed}{
    The seed used to compute the estimates (optional). If not provided, the
    method tries to determine the \code{seed} automatically.
}
  \item{target.data}{
    A list containing the data of the target margins. Each component of the 
    list is an array storing a margin. The list order must follow the one 
    defined in \code{target.list}. Note that the cells of the arrays must be 
    non-negative (and can even be NA if \code{method = ipfp}) (optional). 
    If not provided, the method tries to dermine \code{target.data}
    automatically.
}
  \item{target.list}{
    A list of the target margins provided in \code{target.data}. Each
    component of the list is an array whose cells indicates which dimension the 
    corresponding margin relates to (optional). If not provided, the method
    tries to determine \code{target.list} automatically.
}
  \item{replace.zeros}{
    If 0-cells are to be found, then they are replaced with this value.
}
  \item{\dots}{
    Not used.
}
}
\details{
The test is formally expressed as:
\deqn{H_0 ~ : ~ h(\pi) = 0 \quad vs \quad H_1~:~h(\pi) \neq 0}{%
H0 : h(p.seed) = 0 vs  H1 : h(p.seed) != 0}
where \eqn{\pi}{p.seed} is the vector of the seed probabilities and 
\eqn{h(x) = A^T x - m}{h(x) = t(A) * x - m} with \eqn{A} and \eqn{m} being
respectively the marginal matrix and the margins vector of the estimation
problem.

The three statistics are then defined as:
\itemize{
  \item{Wilk's log-likelihoold ratio}{
    \deqn{G^2 = 2 \sum x_i \ln \frac{\pi_i}{\hat{\pi}_i}}{%
    G2 = 2 * \sum ( x.seed * ln (p.seed / p.hat) )}
  }
  \item{Wald's statistic}{
    \deqn{W^2 = h(x)^T ( H^T_x D_x H_x)^{-1} h(x)}{%
    W2 = t(h(x)) * inv(t(H) * diag(x.seed) * H) * h(x.seed)}
  }
  \item{Pearson Chi-square}{
    \deqn{\chi^2 = (x - n \hat{\pi})^T D^{-1}_{n\hat{\pi}} (x - n \hat{\pi})}{%
    X2 = t(x.seed - n * p.hat) * inv(diag(n * p.hat)) * (x.seed - n * p.hat)}
  }
}
where \eqn{x}{x.seed} is the vectorization of the \code{seed}, 
\eqn{n = \sum x_i}{n = sum(x.seed)}, \eqn{D_v}{diag(v)} is a diagonal matrix 
derived from the vector \eqn{v}{v} and \eqn{H}{H} denotes the Jacobian evaluated 
in \eqn{\hat{\pi}}{p.hat} (the vector of the estimated probabilities) of the 
function \eqn{h(x)}{h(x)}.

The degrees of freedom for these statistics corresponds to the number of
components in \eqn{m}{m}.
}
\value{
A list whose elements are detailed below.
  \item{G2}{
    The Log-likelihood statistic.
}
  \item{W2}{
    The Wald statistic.
}
  \item{X2}{
    The Pearson chi-squared statistic.
}
  \item{stats.df}{
    The degrees of freedom for the \code{G2}, \code{W2} and 
    \code{X2} statistics.
}
}
\references{
 Lang, J.B. (2004) 
 Multinomial-Poisson homogeneous models for contingency tables. 
 \emph{Annals of Statistics} 32(1): 340-383.
}
\author{
Johan Barthelemy

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\seealso{
\code{\link[mipfp:Estimate]{Estimate}} function to create an object of class
\code{mipfp} and to update an initial multidimensional array with respect to 
given constraints.
\code{\link[mipfp:summary.mipfp]{summary.mipfp}} can also retrieve the
statistics and their associated p-values.
}
\examples{
# loading the data
data(spnamur, package = "mipfp")
# subsetting the data frame, keeping only the first 3 variables
spnamur.sub <- subset(spnamur, select = Household.type:Prof.status)
# true table
true.table <- table(spnamur.sub)
# extracting the margins
tgt.v1        <- apply(true.table, 1, sum)
tgt.v1.v2     <- apply(true.table, c(1,2), sum)
tgt.v2.v3     <- apply(true.table, c(2,3), sum)
tgt.list.dims <- list(1, c(1,2), c(2,3))
tgt.data      <- list(tgt.v1, tgt.v1.v2, tgt.v2.v3)
# creating the seed, a 10 pct sample of spnamur
seed.df <- spnamur.sub[sample(nrow(spnamur), round(0.10*nrow(spnamur))), ]
seed.table <- table(seed.df)
# applying one fitting method (ipfp)
r.ipfp <- Estimate(seed=seed.table, target.list=tgt.list.dims, 
                   target.data = tgt.data)
# printing the G2, X2 and W2 statistics
print(gof.estimates(r.ipfp))
# alternative way (pretty printing, with p-values)
print(summary(r.ipfp)$stats.gof)
}
\keyword{multivariate}
\keyword{htest}

