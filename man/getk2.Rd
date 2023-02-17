\name{getk2}
\alias{getk2}
\title{
Constant K_2
}
\description{
Evaluates the constant K_2, which is the normalizing constant for f_2.
}
\usage{
getk2(alpha, p)
}
\arguments{
  \item{alpha}{
Parameter in [0,1).
}
  \item{p}{
Parameter >1.
}
}
\details{
Evaluates

K_2 = int_0^1 ( exp(-x^p) - exp(-x) )*x^(-1-alpha) dx.

This is needed to simulate p-RDTS random variables.}
\value{
Returns a positive real number.
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
 getk2(0.5,2.5)
}
