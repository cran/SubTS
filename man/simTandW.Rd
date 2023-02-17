\name{simTandW}
\alias{simTandW}
\title{
Simulation of hitting time and overshoot.
}
\description{
Simulates the hitting time T and the overshoot W of a stable process by implimenting Algorithm 4.1 in Dassios, Qu, and Lim (2020). This is important for simulating other distribution.
}
\usage{
simTandW(alpha)
}
\arguments{
  \item{alpha}{
Parameter in the open interval (0,1).
}
}
\value{
Returns one pair of random numbers. The first is T and the second is W.
}
\references{
A. Dassios, Y. Qu, J.W. Lim (2020). Exact simulation of a truncated Levy subordinator. ACM Transactions on Modeling and Computer Simulation, 30(10), 17.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
simTandW(.6)
}
