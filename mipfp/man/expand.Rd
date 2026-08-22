\name{expand}
\alias{expand}
\alias{expand.default}
\alias{expand.table}

\title{
Expand a Table in a Data Frame
}
\description{
This function takes a multi-dimensionnal contingency table and expands it to a 
data frame containing individual records.
}
\usage{
expand(x, ...)

\method{expand}{table}(x, ...)
}
\arguments{
  \item{x}{
    An object of type \code{table} storing a N-dimensial contingency table.
}
  \item{\dots}{
    Further arguments passed to or from other methods.
}
}
\value{
A data frame of the individual records derived from \code{x}.
}
\references{
  Cookbook for R - 
  \url{http://www.cookbook-r.com/Manipulating_data/Converting_between_data_frames_and_contingency_tables/}
}
\author{
Johan Barthelemy.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\note{
The function is inspired from the "Cookbook for R".

It should also be noted that the cells od \code{x} are rounded before being
expanded in a data frame.
}
\seealso{
\code{\link[base]{expand.grid}} and \code{\link[base]{as.data.frame}}.
}
\examples{
# loading data
data(spnamur, package = "mipfp")
# subsetting the data frame, keeping only the first 3 variables
spnamur.sub <- subset(spnamur, select = Household.type:Prof.status)
# create a contingency table
t <- table(spnamur.sub)
# expand the table to a data frame
t.df <- expand(t)
}
\keyword{manip}
\keyword{datagen}
\keyword{multivariate}
