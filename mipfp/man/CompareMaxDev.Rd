\name{CompareMaxDev}
\alias{CompareMaxDev}
\title{
Comparing deviations of mipfp objects
}
\description{
This function compares either the margins errors from different \code{mipfp} 
objects or the absolute maximum deviation between a given table and the 
estimates in the mipfp objects.
}
\usage{
CompareMaxDev(list.mipfp = list(), true.table = NULL, echo = FALSE)
}
\arguments{
  \item{list.mipfp}{
    The list produced by the function \code{\link[mipfp:Estimate]{Estimate}}.
}
  \item{true.table}{
    When provided, the estimates contained in the \code{mipfp} objects in the 
    list \code{list.mipfp} are compared against this table. It is an optional
    argument.
}
  \item{echo}{
    Verbose parameter. If \code{TRUE}, the function prints what is being 
    compared. Default is \code{FALSE}.
}
}
\value{
A table with as many rows as the number of mipfp objects in \code{list.mipfp}. 
Each row details the margins errors or the maximum absolute deviation of 
one \code{mipfp} object.
}
\author{
Johan Barthelemy

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\seealso{
The estimation function \code{\link[mipfp:Estimate]{Estimate}}.

This function is used by 
\code{\link[mipfp:error.margins.mipfp]{error.margins.mipfp}}.
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
# applying the different fitting methods
r.ipfp <- Estimate(seed=seed.table, target.list=tgt.list.dims, 
                   target.data = tgt.data,  method = "ipfp")
r.ml   <- Estimate(seed = seed.table, target.list = tgt.list.dims, 
                   target.data = tgt.data, method = "ml")
r.chi2 <- Estimate(seed = seed.table, target.list = tgt.list.dims, 
                   target.data = tgt.data, method = "chi2")
r.lsq  <- Estimate(seed = seed.table, target.list = tgt.list.dims, 
                   target.data = tgt.data, method = "lsq")
# print the maximum absolute deviation between targets and generated margins
CompareMaxDev(list(r.ipfp,r.ml,r.chi2,r.lsq), echo = TRUE)
# compute the maximum absolute deviation between the true and estimated tables
CompareMaxDev(list(r.ipfp,r.ml,r.chi2,r.lsq), echo = TRUE, 
              true.table = true.table)
}
\keyword{univar}

