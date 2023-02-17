\name{rTrunS}
\alias{rTrunS}
\title{
Simulation from the truncated stable distribution}
\description{
Simulates from the truncated stable distribution.}
\usage{
rTrunS(n, t, alpha, b = 1, step = 1)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{t}{
Parameter > 0.}
  \item{alpha}{
Parameter in the open interval (0,1).}
  \item{b}{
Parameter > 0.}
  \item{step}{
Tuning parameter. The larger the step, the slower the rejection sampling, but the fewer the number of terms. See Hoefert (2011) or Section 4 in Grabchak (2019).}
}
\details{
Simulates from the truncated stable distribution using Algorithm 4.3 in Dassios, Qu, and Lim (2020). This distribution has Laplace transform

L(z) = exp( t * (alpha/Gamma(1-alpha)) * int_0^b (e^(-xz)-1) x^(-1-alpha) dx),   z>0

and Levy measure

M(dx) = t * (alpha/Gamma(1-alpha)) * x^(-1-alpha) 1(0<x<b) dx. 

Here Gamma() is the gamma function.
}
\value{
Returns a vector of n random numbers.
}
\references{
A. Dassios, Y. Qu, J.W. Lim (2020). Exact simulation of a truncated Levy subordinator. ACM Transactions on Modeling and Computer Simulation, 30(10), 17.

M. Grabchak (2019). Rejection sampling for tempered Levy processes. Statistics and Computing, 29(3):549-558

M. Hofert (2011). Sampling exponentially tilted stable distributions. ACM Transactions on Modeling and Computer Simulation, 22(1), 3.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rTrunS(10, 2, .6)
}
