\name{rTrunGamma}
\alias{rTrunGamma}
\title{Simulation from the truncated gamma distribution}
\description{
Simulates from the truncated gamma distribution.}
\usage{
rTrunGamma(n, t, mu, b = 1, step = 1)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{t}{
Parameter > 0.}
  \item{mu}{
Parameter > 0.}
  \item{b}{
Parameter > 0.}
  \item{step}{
Tuning parameter. The larger the step, the slower the rejection sampling, but the fewer the number of terms. See Hoefert (2011) or Section 4 in Grabchak (2019).
}
}
\details{
Simulates from the truncated gamma distribution. This distribution has Laplace transform

L(z) = exp( t int_0^b (e^(-xz)-1) x^(-1)e^(-mu*x) dx),   z>0

and Levy measure

M(dx) = t x^(-1) e^(-mu*x) 1(0<x<b) dx. 

The simulation is performed by applying rejection sampling (Algorithm 4.4 in Dassios, Qu, Lim (2020)) to the generalized Dickman distribution. We simulate from the latter using Algorithm 3.1 in Dassios, Qu, Lim (2019).
}
\value{
Returns a vector of n random numbers.
}
\references{
A. Dassios, Y. Qu, J.W. Lim (2019). Exact simulation of generalised Vervaat perpetuities. Journal of Applied Probability, 56(1):57-75.

A. Dassios, Y. Qu, J.W. Lim (2020). Exact simulation of a truncated Levy subordinator. ACM Transactions on Modeling and Computer Simulation, 30(10), 17.

M. Grabchak (2019). Rejection sampling for tempered Levy processes. Statistics and Computing, 29(3):549-558

M. Hofert (2011). Sampling exponentially tilted stable distributions. ACM Transactions on Modeling and Computer Simulation, 22(1), 3.}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rTrunGamma(10, 2, 1)
}
