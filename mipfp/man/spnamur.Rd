\name{spnamur}
\alias{spnamur}
\docType{data}
\title{
Synthetic population of Namur (Belgium)
}
\description{
This data drame contains a synthetic population of individuals for Belgian city
of Namur. The attributes details the gender, age class, socio-professional
status, education level and driving license ownership of every synthetic
individual.
}
\usage{data(spnamur)}

\format{
A data frame detailing the synthetic individuals whose columns are described
in the Table below.

\tabular{ll}{
  \bold{Attribute} \tab \bold{Values (levels)} \cr
  Household.type \tab C(couple); F (family with children); I (isolated); N (non
                      family) \cr
  Gender \tab F (female); H (male) \cr
  Prof.status  \tab A (active); E (student); I (inactive) \cr
  Education.level  \tab O (none); P (primary); S (high school); 
                        U (higher education) \cr
  Driving.license \tab O (no); P (yes) \cr
  Age.class	\tab 0 (0-5); 1 (6-17); 2 (18-39); 3 (40-59); 4 (60+) \cr  
}
}
\source{
  VirtualBelgium - \url{http://virtualbelgium.sourceforge.net}
}
\references{
  Barthelemy, J. and Toint, P.L. (2013)
  Synthetic population generation without a sample
  \emph{Transportation Science} 47 (2): 266-279
}
\examples{
data(spnamur)
# generating the contingency table of the synthetic population
table(spnamur)
}
\keyword{datasets}
