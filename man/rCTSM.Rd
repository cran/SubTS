\name{rCTSM}
\alias{rCTSM}
\title{
Simulation of Fully Left Skewed CTS Distributions
}
\description{
Simulates from fully left skewed classical tempered stable (CTS) distributions in both the finite and infinite variation cases. When alpha=0 this is the negative of the gamma distribution.
}
\usage{
rCTSM(n, alpha, c, ell, epsilon = 0.5, p = 0.5)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{alpha}{
Parameter in [0,2).
}
  \item{c}{
Parameter >0.
}
  \item{ell}{
Tempering parameter >0.
}
  \item{epsilon}{
Tuning parameter in (0,1). Only for alpha>=1.
}
  \item{p}{
  Tuning parameter in (0,1). Only for alpha>=1.
  }
}
\details{
Simulates from a fully left skewed CTS distribution. The distribution has moment generating function

M(z) = exp( c int_(-infty)^0 (e^(xz)-1)e^(-x/ell) x^(-1-alpha) dx)

and Levy measure

M(dx) = c e^(-|x|/ell) |x|^(-1-alpha) 1(x<0)dx.
}
\value{
Returns a vector of n random numbers.
}
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.

M. Grabchak (2026). Exact Simulation from Tempered Stable Distributions with Infinite Variation (alpha>=1). <doi 10.48550/arXiv.2604.17732>.
}
\author{Michael Grabchak and Lijuan Cao}

\examples{
rCTSM(20, 1.5, 1, 1)
}
