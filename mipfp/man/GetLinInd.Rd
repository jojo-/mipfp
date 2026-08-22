\name{GetLinInd}
\alias{GetLinInd}
\title{
Extracting the linearly independant columns from a matrix
}
\description{
Extracts the linearly dependant columns of matrix to obtain a matrix of full
rank using QR decomposition.
}
\usage{
GetLinInd(mat, tol = 1e-10)
}
\arguments{
  \item{mat}{
    The matrix possibly containing linearly dependant columns
}
  \item{tol}{
    Rank estimation tolerance. Default is \eqn{1e^{-10}}{1e-10}.
}
}
\value{
  A list containing the new matrix and the index of the selected colums.
  \item{mat.li}{
    A matrix made of the linearly independant columns of \code{mat}.
}
  \item{idx}{
    The index of the selected columns.    
}
}
\references{
 Golub, G. H., Van Loan C. F. (2012)
 \emph{Matrix Computations. Third Edition.}
 Johns Hopkins University Press.
}
\author{
Johan Barthelemy

Maintainer: Johan Barthelemy \email{johan@uow.edu.au}.
}
\seealso{
\code{\link[base:qr]{qr}}.
}
\examples{
# generation of a matrix with linearly dependant columns
A <- matrix(c(1, 2, 3,
              1, 2, 4,
              1, 2, 8), nrow = 3, ncol = 3, byrow = TRUE)
              
# extracting and printing the linearly independant columns
B <- GetLinInd(A)
print(B)
}
\keyword{array}
\keyword{algebra}
