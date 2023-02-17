\name{rGGa}
\alias{rGGa}
\title{
Simulates from the generalized gamma distribution}
\description{
Simulates from the generalized gamma distribution.
}
\usage{
rGGa(n, a, p, b)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{a}{
Parameter >0.
}
  \item{p}{
Parameter >0.
}
  \item{b}{
Parameter >0.}
}
\details{
Simulates from the generalized gamma distribution with density

g(x) = exp(-b*x^p)*x^(a-1)/K_3,   x>0,

where K_3 is a normalizing constant. The mathodology is explained in Section 4 of Grabchak (2021). This distribution is needed to simulate p-RDTS random variables with negative alpha values.}
\value{
Returns a vector of n random numbers.
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.

E.W. Stacy (1962) A generalization of the gamma distribution. Annals of Mathematical Statistics, 33(3):1187-1192.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rGGa(20, .5, 2, 2)
}

