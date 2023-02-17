\name{rDickman}
\alias{rDickman}
\title{Simulation from the generalized Dickman distribution}
\description{Simulates from the generalized Dickman distribution using Algorithm 3.1 in Dassios, Qu, and Lim (2019).}
\usage{
rDickman(n, t, b = 1)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{t}{
Parameter > 0.
}
  \item{b}{
Parameter > 0.}
}
\details{
Simulates from the generalized Dickman distribution by using Algorithm 3.1 in Dassios, Qu, and Lim (2019). This distribution has Laplace transform

L(z) = exp( t int_0^b (e^(-xz)-1) x^(-1) dx),   z>0

and Levy measure

M(dx) = t x^(-1) 1(0<x<b) dx.

When b=1 and t=1, this is the Dickman distribution.
}
\value{
Returns a vector of n random numbers.
}
\references{
A. Dassios, Y. Qu, J.W. Lim (2019). Exact simulation of generalised Vervaat perpetuities. Journal of Applied Probability, 56(1):57-75.

M. Penrose and A. Wade (2004). Random minimal directed spanning trees and Dickman-type distributions. Advances in Applied Probability, 36(3):691-714.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rDickman(10, 1)
}

