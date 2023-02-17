\name{dF1}
\alias{dF1}
\title{Pdf for f_1}
\description{
Evaluates the pdf f_1(x) intruduced in Grabchak (2021).}
\usage{
dF1(x, a, p)
}
\arguments{
  \item{x}{
Vector of real numbers.
}
  \item{a}{
Parameter >=0.
}
  \item{p}{
Parameter >1.
}
}
\details{
Evaluates the pdf

f_1(x) = exp(-x^p)*x^(-1-a)/K_1,   x>1

where K_1 is a normalizing constant. This is distribution is needed to simulate p-RDTS random variables.
}
\value{
Returns a vector of real numbers corresponding to the values of f_1(x).
}
\references{
M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
  x = (10:20)/10
  dF1(x, .5, 2)
}
