\name{dF2}
\alias{dF2}
\title{Pdf for f_2}
\description{
Evaluates the pdf f_2(x) intruduced in Grabchak (2021).
}
\usage{
dF2(x, a, p)
}
\arguments{
  \item{x}{
Vector of real numbers.
}
  \item{a}{
Parameter in [0,1).
}
  \item{p}{
Parameter >1.
}
}
\details{
Evaluates the pdf

f_2(x) = (exp(-x^p) - exp(-x))*x^(-1-a)/K_2,   0<x<1

where K_2 is a normalizing constant. This distribution is needed to simulate p-RDTS random variables.}
\value{
Returns a vector of real numbers corresponding to the values of f_2(x).
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
  x = (0:10)/10
  dF2(x, .5, 1.5)
}

