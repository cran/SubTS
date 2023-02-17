\name{getk1}
\alias{getk1}
\title{
Constant K_1
}
\description{
Evaluates the constant K_1, which is the normalizing constant for f_1.
}
\usage{
getk1(alpha, p)
}
\arguments{
  \item{alpha}{
Parameter >=0.
}
  \item{p}{
Parameter >1.
}
}
\details{
Evaluates

K_1 = int_1^infty exp(-x^p)*x^(-1-alpha) dx.

This is needed to simulate p-RDTS random variables.
}
\value{
Returns a positive real number.
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
 getk1(1.5,2.5)
}
