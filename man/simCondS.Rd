\name{simCondS}
\alias{simCondS}
\title{
Simulation from a conditioned stable distribution.
}
\description{
Implements Algorithm 4.2 in Dassios, Qu, and Lim (2020) to simulate from a stable distribution conditioned on an appropriate event.}
\usage{
simCondS(t, alpha)
}
\arguments{
  \item{t}{
Parameter > 0.
}
  \item{alpha}{
Parameter in the open interval (0,1).
}
}
\details{
Implements Algorithm 4.2 in Dassios, Qu, and Lim (2020) to simulate from a stable distribution conditioned on an appropriate event. There are some typos in this algorithm, which are corrected in Grabchak (2021). These random variables are needed to simulate truncated stable, truncated tempered stable, and p-RDTS random variables. 
}
\value{
Returns one random number.
}
\references{
A. Dassios, Y. Qu, J.W. Lim (2020). Exact simulation of a truncated Levy subordinator. ACM Transactions on Modeling and Computer Simulation, 30(10), 17.

M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
simCondS(2, .7)
}
