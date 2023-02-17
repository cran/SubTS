\name{SubTS-package}
\alias{SubTS-package}
\alias{SubTS}
\docType{package}
\title{
\packageTitle{SubTS}
}
\description{
\packageDescription{SubTS}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{SubTS}
\packageIndices{SubTS}
}
\author{
\packageAuthor{SubTS}

Maintainer: \packageMaintainer{SubTS}
}
\references{
A. Dassios, Y. Qu, J.W. Lim (2019). Exact simulation of generalised Vervaat perpetuities. Journal of Applied Probability, 56(1):57-75.

A. Dassios, Y. Qu, J.W. Lim (2020). Exact simulation of a truncated Levy subordinator. ACM Transactions on Modeling and Computer Simulation, 30(10), 17.

M. Grabchak (2016). Tempered Stable Distributions: Stochastic Models for Multiscale Processes. Springer, Cham.

M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\keyword{Tempered stable distributions}
\keyword{Subordinators}
\keyword{Dickman distribution}
\keyword{Truncated distributions}
\examples{
rPRDTS(20, 2, 1, .7, 2)
rPRDTS(20, 2, 1, 0, 2)
rPRDTS(20, 2, 1, -.7, 2)
rDickman(10, 1)
rTrunGamma(10, 2, 1)
rPGamma(20, 2, 2, 2)
rTrunS(10, 2, .6)
rTrunTS(10, 2, 2, .6)
}
