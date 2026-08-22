\name{coef.mipfp}
\alias{coef.mipfp}

\title{
Extract the coefficients of the estimates from an object of class mipfp
}
\description{
This method extracts the coefficients of estimates of an mipfp object.
}
\usage{
\method{coef}{mipfp}(object, prop = FALSE, ...)
}
\arguments{
  \item{object}{
    An object of class \code{mipfp}
}
  \item{prop}{
    If this Boolean is set to \code{TRUE} then the method will return
    the estimated probabilities. Otherwise, it will return the estimated counts.
    Default is \code{False}.
}
  \item{...}{
    Not used.  
}
}
\value{
Coefficients of the estimates extracted from the \code{mipfp} object 
\code{object}.
}
\references{
 Chambers, J. M. and Hastie, T. J. (1992)
 \emph{Statistical Models in S.}
 Wadsworth & Brooks/Cole.
}
\author{
Johan Barthelemy.

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}

\seealso{
\code{\link[stats:coef]{coef}}.
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
# estimating a table using ipfp
r.ipfp <- Estimate(seed=seed.table, target.list=tgt.list.dims, 
                   target.data = tgt.data)
# exctracting and printing the coefficient of the estimates
print(coef(r.ipfp))
}
\keyword{multivariate}

