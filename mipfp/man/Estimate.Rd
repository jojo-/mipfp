\name{Estimate}
\alias{Estimate}
\alias{print.mipfp}
\title{
Update an N-way table given target margins
}
\description{
This function provides several estimating methods to up multiway table 
(referred as the seed) subject to known constrains/totals: Iterative 
proportional fitting procedure (ipfp), maximum likelihood method (ml), minimum 
chi-squared (chi2) and weighted least squares (lsq). Note that the targets can 
also be multi-dimensional.
}
\usage{
Estimate(seed, target.list, target.data, method = "ipfp", keep.input = FALSE,
         ...)
}
\arguments{
  \item{seed}{
    The initial multi-dimensional array to be updated. Each cell must
    be non-negative if \code{method} is \code{ipfp} or strictly positive
    when \code{method} is \code{ml}, \code{lsq} or \code{chi2}.
}
  \item{target.list}{
    A list of dimensions of the marginal target constrains in
    \code{target.data}. Each component of the list is an array whose cells 
    indicate which dimension the corresponding margin relates to.
}
  \item{target.data}{
    A list containing the data of the target marginal tables. Each
    component of the list is an array storing a margin.
    The list order must follow the ordering defined in \code{target.list}.
    Note that the cells of the arrays must be non-negative.
}
  \item{method}{
    An optional character string indicating which method is to be used to
    update the \code{seed}. This must be on of the strings "ipfp", "ml",
    "chi2", or "lsq". Default is "ipfp".
}
  \item{keep.input}{
    A Boolean indicating if \code{seed}, \code{target.data} and 
    \code{target.list} when set to \code{TRUE}.
}
  \item{\dots}{
    Additionals argument that can be passed to the functions 
    \link[mipfp:Ipfp]{Ipfp} and 
    \link[mipfp:ObtainModelEstimates]{ObtainModelEstimates}. See their 
    respective documentation for more details.
}
}
\value{
An object of class \code{mipfp} is a list containing at least the 
following components:
  \item{x.hat}{
    An array with the same dimension of \code{seed} whose margins match those
    specified in \code{target.list}.
}
  \item{p.hat}{
    An array with the same dimension of \code{x.hat} containing the updated
    cell probabilities, i.e. \code{x.hat / sum(x.hat)}.
}
  \item{error.margins}{
    A list returning, for each margin, the absolute maximum deviation between 
    the desired and generated margin.
}
  \item{conv}{
    A boolean indicating whether the algorithm converged to a solution.
}
 \item{evol.stp.crit}{
    The evolution of the stopping criterion over the iterations (if selected
    \code{method} is \code{"ipfp"})).
}
  \item{solnp.res}{
    The estimation process uses the \code{solnp} optimisation function from 
    the R package Rsolnp and \code{solnp.res} is the corresponding object 
    returned by the solver (if selected \code{method} is not \code{"ipfp"}).
}
  \item{method}{
    The selected method for estimation.
}
  \item{call}{
    The matched call.  
}

The will be also added if \code{keep.input} has been set to \code{TRUE}: 
\code{seed}, \code{target.data}, \code{target.list}.

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
  
 Little, R. J., Wu, M. M. (1991)
 Models for contingency tables with known margins when target and sampled 
 populations differ.
 \emph{Journal of the American Statistical Association} 86 (413): 87-95.
 
 Lang, J.B. (2004) 
 Multinomial-Poisson homogeneous models for contingency tables. 
 \emph{Annals of Statistics} 32(1): 340-383.
  
 Stephan, F. F. (1942).
 Iterative method of adjusting frequency tables when expected margins are known.
 \emph{Annals of Mathematical Statistics} 13 (2): 166-178.
}
\author{
Johan Barthelemy.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\note{
It is important to note that if the margins given in \code{target.list} are
not consistent (i.e. the sums of their cells are not equals), the input data
is then normalised by considering probabilities instead of frequencies:
\itemize{
  \item the cells of the seed are divided by \code{sum(seed)};
  \item the cells of each margin \code{i} of the list \code{target.data} are 
        divided by \code{sum(target.data[[i]])}.
}
}
\seealso{
See the functions \code{\link[mipfp:Ipfp]{Ipfp}} and 
\code{\link[mipfp:ObtainModelEstimates]{ObtainModelEstimates}}
for more details about the estimation process.

\code{\link[mipfp:summary.mipfp]{summary.mipfp}} for summaries,  
\code{\link[mipfp:vcov.mipfp]{vcov.mipfp}} for the (asymptotic) covariance of 
the estimates and \code{\link[mipfp:gof.estimates.mipfp]{gof.estimates.mipfp}}
for testing if the seed agrees with the targets.

The genereric functions \code{\link[mipfp:print.mipfp]{print}} and
\code{\link[mipfp:coef.mipfp]{coef}}.
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
print(r.ipfp)
}
\keyword{models}
\keyword{multivariate}
\keyword{array}
\concept{iterative proportional fitting procedure}
\concept{ipfp}
\concept{multiway contingency table update}
\concept{multiway contingency table estimation}
\concept{multivariate contingency table estimation}
\concept{RAS algorithm}
\concept{matrix raking}
\concept{matrix scaling}

