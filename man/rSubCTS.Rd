\name{rSubCTS}
\alias{rSubCTS}
\title{Simulates of CTS subordinators}
\description{
Simulates from classical tempered stable (CTS) distributions. When alpha=0 this is the gamma distribution.
}
\usage{
rSubCTS(n, alpha, c, ell, method = NULL)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{alpha}{
Parameter in [0,1).
}
  \item{c}{
Parameter >0
}
  \item{ell}{
Tempering parameter >0
}
  \item{method}{
Parameter used by retstable in the copula package. When NULL restable selects the best method.}
}
\details{
Simulates a CTS subordinator. The distribution has Laplace transform

L(z) = exp( c int_0^infty (e^(-xz)-1)e^(-x/ell) x^(-1-alpha) dx),   z>0

and Levy measure

M(dx) = c e^(-x/ell) x^(-1-alpha) 1(x>0)dx.
}
\value{
Returns a vector of n random numbers.
}
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.
}
\author{Michael Grabchak and Lijuan Cao}
\note{
Uses the method retstable in the copula package.
}
\examples{
rSubCTS(20, .7, 1, 1)
}
