\name{Qaqish}
\alias{Qaqish}
\docType{data}
\title{
Qaqish
}
\description{
The data set provides the odds ratios and correlations as measures of 
associations of the binary outcome impaired pulmonary function for a family of 
four with two parents and two siblings.

These correlations and odds ratios are obtained from Qaqish et al. (2012) based
on a regression analysis of a common data set of parents and siblings with 
chronic obstructive pulmonary disease and their controls.
}
\usage{data(Qaqish)}
\format{
A list \code{Qaqish} containing 2 elements:
\enumerate{
  \item \code{cr} : the correlation matrix;
  \item \code{or} : the odd ratios matrix.
}
}
\source{
  Qaqish, B. F., Zink, R. C., and Preisser, J. S. (2012).
  Orthogonalized residuals for estimation of marginally specified association
  parameters in multivariate binary data.
  \emph{Scandinavian Journal of Statistics} 39, 515-527.
}
\examples{
data(Qaqish)
print(Qaqish$or)
print(Qaqish$cr)
}
\keyword{datasets}
