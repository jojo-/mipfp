\name{ComputeA}
\alias{ComputeA}
\title{
Computes the marginal matrix A and margins vector m of an estimation problem
}
\description{
Given a set of marginal target constraints and the dimension of the array
\eqn{X}{X} to wich the targets relate to, this function computes the matrix 
\eqn{A}{A} of full rank and vector \eqn{m}{m} such that
\deqn{A^T \pi = (m, 1)^T}{t(A) * p = t(c(m, 1)) }
where vector \eqn{m}{m} contains all components but one of every target and
\eqn{\pi}{p} is a vector of the (unknow) components of \eqn{X}.
}
\usage{
ComputeA(dim.arr, target.list, target.data)
}
\arguments{
  \item{dim.arr}{
    The dimension of the array \eqn{X} to which the margins are applied.
}
  \item{target.list}{
    A list of the target margins provided in target.data. Each 
    component of the list is an array whose cells indicates which
    dimension the corresponding margin relates to.
}
  \item{target.data}{
    A list containing the data of the target margins. Each 
    component of the list is an array storing a margin. The list 
    order must follow the one defined in target.list. Note that 
    the cells of the arrays must be non-negative.
}
}
\value{
A list whose elements are defined below.
  \item{marginal.matrix}{
    The marginal matrix.
}
  \item{margins}{
    A vector containing the margins associated with A.
}
  \item{df}{
    The degree of freedom of the problem.
}
}
\author{
Johan Barthelemy

Maintainer: Johan Barthelemy <johan@uow.edu.au>.
}
\seealso{
\code{\link[cmm:MarginalMatrix]{MarginalMatrix}}.
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
# computing the associated marginal matrix and margins vector
res.marg <- ComputeA(dim(seed.table), tgt.list.dims, tgt.data)
print(res.marg)
}
\keyword{array}
\keyword{algebra}
