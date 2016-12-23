\name{flat}
\alias{flat}
\alias{flat.default}
\alias{flat.table}
\alias{flat.array}
\alias{flat.matrix}

\title{
Flatten a table, array or matrix
}
\description{
This function takes a multidimensional object and flattens it for a pretty 
printing. The row names are the concatenation of the original dimension names 
while the only column stores the initial data of the object.
}
\usage{
\method{flat}{array}(x, sep = ".", label = "value", l.names = 0, ...)

\method{flat}{table}(x, sep = ".", label = "value", l.names = 0, ...)

\method{flat}{matrix}(x, sep = ".", label = "value", l.names = 0, ...)
}
\arguments{
  \item{x}{ 
    An array, table or matrix.
}
  \item{sep}{
    The separator used to concatenate the dimension names.
}
  \item{label}{
    The name of the column storing the data.
}
  \item{l.names}{
    If set to a value greater than 0, then the dimnames will be shorten to
    a length of \code{l.names} characters.
}
  \item{\dots}{
    Not used.
}
}
\value{
An array containing a flattened version of \code{x}.
}
\author{
Johan Barthelemy.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\note{
The function is inspired from the function \code{wrap.array} from the package 
\code{R.utils} written by Henrik Bengtsson.
}
\seealso{
The function \code{wrap.array} from the \code{R.utils} package
(\url{https://cran.r-project.org/package=R.utils}).
}
\examples{
# loading the data and saving in a 3D-table
data(spnamur, package = "mipfp")
spnamur.sub <- subset(spnamur, select = Household.type:Prof.status)
tab <- table(spnamur.sub)

# flattening the table
tab.flat <- flat(tab)
print(tab.flat)
}
\keyword{methods}
\keyword{programming}
