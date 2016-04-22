\name{error.margins}
\alias{error.margins}
\alias{error.margins.default}
\alias{error.margins.mipfp}
\title{
Extracts the deviation between every target and generated margin
}
\description{
This method returns the maximum deviation between each generated and 
desired margins of the input argument. It corresponds to the absolute
maximum deviation between each target margin used to generate the estimates
in the mipfp object and the generated one.
}
\usage{
\method{error.margins}{mipfp}(object, ...)
}
\arguments{
  \item{object}{
    An object of class \code{mipfp}.
}
  \item{\dots}{
    Further arguments passed to or from other methods. See 
    \code{\link[mipfp:CompareMaxDev]{CompareMaxDev}}.
}
}
\value{
An array containing the absolute maximum deviations for each margin.
}
\author{
Johan Barthelemy

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\note{
It is an alias for \code{\link[mipfp:CompareMaxDev]{CompareMaxDev}} when
only one object is passed to the function and the verbose parameter is set to
\code{FALSE}.
}
\seealso{
The estimation function \code{\link[mipfp:Estimate]{Estimate}}.

This function relies on 
\code{\link[mipfp:CompareMaxDev]{CompareMaxDev}}.
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
# creating the seed, a 10% sample of spnamur
seed.df <- spnamur.sub[sample(nrow(spnamur), round(0.10*nrow(spnamur))), ]
seed.table <- table(seed.df)
# applying a fitting method
r.ipfp <- Estimate(seed=seed.table, target.list=tgt.list.dims, 
                   target.data = tgt.data,  method = "ipfp")
# print the maximum absolute deviation between targets and generated margins
print(error.margins(r.ipfp))
}
\keyword{univar}

