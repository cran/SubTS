\name{rF1}
\alias{rF1}
\title{Simulation from f_1}
\description{
Simulates from the pdf f_1(x) intruduced in Grabchak (2021).
}
\usage{
rF1(n, a, p)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{a}{
Parameter >=0.
}
  \item{p}{
Parameter >1.
}
}
\details{
Uses Algorithm 1 in Grabchak (2021) to simulate from the pdf

f_1(x) = exp(-x^p)*x^(-1-a)/K_1,   x>1,

where K_1 is a normalizing constant. This is needed to simulate p-RDTS random variables.}
\value{
Returns a vector of n random numbers.
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rF1(10, .7, 2.5)
}
