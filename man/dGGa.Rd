\name{dGGa}
\alias{dGGa}
\title{Pdf of the generalized gamma distribution}
\description{
Evaluates the pdf of the generalized gamma distribution.
}
\usage{
dGGa(x, a, p, b)
}
\arguments{
  \item{x}{
Vector of real numbers.
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
Evaluates the pdf of the generalized gamma distribution with density

g(x) = exp(-b*x^p)*x^(a-1)/K_3,   x>0,

where K_3 is a normalizing constant. This distribution is needed to simulate p-RDTS random variables with negative alpha values.
}
\value{
Returns a vector of real numbers corresponding to the values of g(x).
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.

E.W. Stacy (1962) A generalization of the gamma distribution. Annals of Mathematical Statistics, 33(3):1187-1192.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
  x = (0:20)/10
  dGGa(x, 2.5, 1.5, 3.1)
}
