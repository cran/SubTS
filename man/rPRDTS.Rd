\name{rPRDTS}
\alias{rPRDTS}
\title{
Simulation from p-RDTS distributions.
}
\description{
Simulates from p-rapidly decreasing tempered stable (p-RDTS) distributions.
}
\usage{
rPRDTS(n, t, mu, alpha, p, step = 1)
}
\arguments{
  \item{n}{
Number of observations.
}
  \item{t}{
Parameter >0.
}
  \item{mu}{
Parameter >0.
}
  \item{alpha}{
Parameter in (-infty,1)
}
  \item{p}{
Parameter >1 if 0<=alpha<1, >0 if alpha<0.
}
\item{step}{
Tuning parameter. The larger the step, the slower the rejection sampling, but the fewer the number of terms. See Hoefert (2011) or Section 4 in Grabchak (2019).}
}
\details{
Simulates from a p-RDTS distribution. When alpha >=0, uses Theorem 1 in Grabchak (2021) and when alpha<0 uses the method in Section 4 of Grabchak (2021). This distribution has Laplace transform

L(z) = exp( t int_0^infty (e^(-xz)-1)e^(-(mu*x)^p) x^(-1-alpha) dx ),   z>0

and Levy measure

M(dx) = t e^(-(mu*x)^p) x^(-1-alpha) 1(x>0)dx.
}
\value{
Returns a vector of n random numbers.
}
\references{
M. Grabchak (2019). Rejection sampling for tempered Levy processes. Statistics and Computing, 29(3):549-558

M. Grabchak (2021). An exact method for simulating rapidly decreasing tempered stable distributions. Statistics and Probability Letters, 170: Article 109015.

M. Hofert (2011). Sampling exponentially tilted stable distributions. ACM Transactions on Modeling and Computer Simulation, 22(1), 3.
}
\author{Michael Grabchak and Lijuan Cao}
\examples{
rPRDTS(20, 2, 1, .7, 2)
rPRDTS(20, 2, 1, 0, 2)
rPRDTS(20, 2, 1, -.7, 2)
}