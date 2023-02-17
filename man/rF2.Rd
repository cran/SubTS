\name{rF2}
\alias{rF2}
\title{Simulation from f_2}
\description{
Simulates from the pdf f_2(x) intruduced in Grabchak (2021).
}
\usage{
rF2(n, a, p)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{a}{
Parameter in [0,1).
}
  \item{p}{
Parameter >1.
}
}
\details{
Uses Algorithm 2 in Grabchak (2021) to simulate from the pdf

f_2(x) = (exp(-x^p) - exp(-x))*x^(-1-a)/K_2,   0<x<1

where K_2 is a normalizing constant. This is needed to simulate p-RDTS random variables.}
\value{
Returns a vector of n random numbers.
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rF2(10, .7, 2.5)
}
