\name{dSubCTS}
\alias{dSubCTS}
\title{PDF of CTS subordinator}
\description{
Evaluates the pdf of the classical tempered stable (CTS) subordinator. When alpha=0 this is the pdf of the gamma distribution.}
\usage{
dSubCTS(x, alpha, c, ell)
}
\arguments{
  \item{x}{
Vector of real numbers.
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
}
\details{
Returns the pdf of a classical tempered stable subordinator. The distribution has Laplace transform

L(z) = exp( c int_0^infty (e^(-xz)-1)e^(-x/ell) x^(-1-alpha) dx),   z>0

and Levy measure

M(dx) = c e^(-x/ell) x^(-1-alpha) 1(x>0)dx.
}
\value{
Returns a vector of real numbers corresponding to the values of pdf.
}
\references{
M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.
}
\author{Michael Grabchak and Lijuan Cao}
\note{
Uses the method dtweedie in the Tweedie package.
}
\examples{
  x = (0:20)/10
  dSubCTS(x, .5, 1, 1.5)
}
