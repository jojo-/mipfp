\name{Array2Vector}
\alias{Array2Vector}
\title{
Transforming an array to a vector
}
\description{
Transform a N-dimensional array a to vector. The transformation is done
assuming that the last index of the array moves fastest. For instance, an array
\eqn{a} of dimensions (2,2,2) will produce the vector
\eqn{v = (a_{111}, a_{112}, a_{113}, a_{121}, a_{122},\ldots,a_{333})}{
v = ( a[1,1,1], a[1,1,2], a[1,1,3], a[1,2,1], a[1,2,2],\ldots,a[3,3,3])}.
}
\usage{
Array2Vector(arr)
}
\arguments{
  \item{arr}{
    The array to be transformed.
}
}
\value{
A vector filled with the data of the input array \code{arr}.
}
\author{
Thomas Suesse.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\seealso{
The inverse transformation is performed with the function
\code{\link[mipfp:Vector2Array]{Vector2Array}}
}
\examples{
# generating an array of dimension (3,3,3)
a <- array(seq(1:27),dim=c(3,3,3))
# transforming it into a vector
v <- Array2Vector(a)
}

\keyword{vector}
\keyword{array}
\keyword{manip}
\keyword{utilities}
