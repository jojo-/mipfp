\name{Ipfp}
\alias{Ipfp}
\title{
Multidimensional Iterative Proportional Fitting
}
\description{
This function implements the iterative proportional fitting (IPFP) procedure. 
This procedure updates an initial N-dimensional array (referred as the seed) 
with respect to given target marginal distributions. Those targets can also 
be multi-dimensional. This procedure is also able to 
estimate a (multi-dimensional) contingency table (encoded as an array) matching 
a given set of (multi-dimensional) margins. In that case, each cell of the seed 
must simply be set to 1.

The IPFP is also known as the RAS algorithm in economics and matrix raking or 
matrix scaling in computer science.
}
\usage{
Ipfp(seed, target.list, target.data, print = FALSE, iter = 1000, tol = 1e-10,
     tol.margins = 1e-10, na.target = FALSE)
}
\arguments{
  \item{seed}{
    The initial multi-dimensional array to be updated. Each cell must
    be non-negative.
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
  \item{print}{
    Verbose parameter: if TRUE prints the current iteration number
    and the associated value of the stopping criterion. Default is FALSE.
}
  \item{iter}{
    Stopping criterion. The maximum number of iteration allowed; 
    must be greater than 0. Default is 1000.
}
  \item{tol}{
    Stopping criterion. If the maximum absolute difference between two iteration
    is lower than the value specified by \code{tol}, then ipfp has reached 
    convergence; must be greater than 0. Default is \eqn{1e^{-10}}{1e-10}.
}
  \item{tol.margins}{
    Tolerance for the margins consistency. Default is \eqn{1e^{-10}}{1e-10}.
}
  \item{na.target}{
    If set to TRUE, allows the targets to have NA cells. Note that in that 
    particular case the margins consistency is not checked.
}
}
\value{
A list containing the final updated array as well as other convergence
informations.
  \item{x.hat}{
    An array with the same dimension of \code{seed} whose margins match those
    specified in \code{target.list}.
}
  \item{p.hat}{
    An array with the same dimension of \code{x.hat} containing the updated
    cell probabilities, i.e. \code{x.hat / sum(x.hat)}.
}
  \item{evol.stp.crit}{
    The evolution of the stopping criterion over the iterations.
}
  \item{conv}{
    A boolean indicating whether the algorithm converged to a solution.
}
  \item{error.margins}{
    A list returning, for each margin, the absolute maximum deviation between 
    the desired and generated margin.
}
  \item{method}{
    The selected method for estimation (here it will always be \code{ipfpf}).  
}
  \item{call}{
    The matched call.  
}
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
The documentation of \code{\link[mipfp:IpfpCov]{IpfpCov}} provide details on the
the covariance matrices determination.
 
\code{\link[mipfp:ObtainModelEstimates]{ObtainModelEstimates}} for alternatives 
to the IPFP.
}
\examples{
# Example 1: 2-way table (V1,V2) of dim=(2,2)
# generating an intial 2-way table to be updated
seed.2d <- array(1,dim=c(2,2))
# desired targets (margins) : V1 and V2
target.row <- c(50,50)
target.col <- c(30,70)
# storing the margins in a list
tgt.data.2d <- list(target.col, target.row)
# list of dimensions of each marginal constrain
tgt.list.2d <- list(1,2)
# calling the Ipfp function
res.2d <- Ipfp(seed.2d, tgt.list.2d, tgt.data.2d)

# Example 2: 3-way table (V1,V2,V3) of dim=(2,4,2)
# seed
seed.3d <- array(1,c(2,4,2))
seed.3d[1,1,1] <- 4
seed.3d[1,3,1] <- 10
seed.3d[1,4,2] <- 6
# desired targets (margins) : V1 and (V2,V3)
target.V1 <- c(50, 16)
target.V2.V3 <- array(4, dim=c(4,2))
target.V2.V3[1,1] <- 10
target.V2.V3[3,1] <- 22
target.V2.V3[4,2] <- 14
# list of dimensions of each marginal constrain
tgt.data.3d <- list(target.V1, target.V2.V3)
# storing the description of target data in a list
tgt.list.3d <- list( 1, c(2,3) )
# calling the Ipfp function
res.3d <- Ipfp(seed.3d, tgt.list.3d, tgt.data.3d, iter=50, print=TRUE, tol=1e-5)

# Example 3: 2-way table (V1,V2) of dim=(2,3) with missing values in the targets
# generating an intial 2-way table to be updated
seed.2d.na <- array(1,dim=c(2,3))
# desired targets (margins) : V1 and V2
target.row.na <- c(40,60)
target.col.na <- c(NA,10,NA)
# storing the margins in a list
tgt.data.2d.na <- list(target.row.na, target.col.na)
# storing the description of target data in a list
tgt.list.2d.na <- list(1,2)
# calling the Ipfp function
res.2d.na <- Ipfp(seed.2d.na, tgt.list.2d.na, tgt.data.2d.na, na.target=TRUE)
}
\concept{iterative proportional fitting procedure}
\concept{ipfp}
\concept{multiway contingency table update}
\concept{multivariate contingency table estimation}
\concept{RAS algorithm}
\concept{matrix raking}
\concept{matrix scaling}
\keyword{models}
\keyword{multivariate}
\keyword{array}