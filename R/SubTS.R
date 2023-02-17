library(gsl)
library(copula)
library(tweedie)

dGGa<- function(x, a, p, b) {
  if(a <= 0)
    stop("a must be positive")
  if(p <= 0)
    stop("p must be positive")
  if(b <= 0)
    stop("b must be positive")
  ret = rep(0,length(x))
  ret[x>0] = exp(-b*x[x>0]^p)*x[x>0]^(a-1)*p*b^(a/p)/gamma(a/p)
  return(ret)
}

rGGa <- function(n, a, p, b) {
  if(n < 1 || is.integer(n))
    stop("n must be an integer greater than or equal to 1")
  if(a <= 0)
    stop("a must be positive")
  if(p <= 0)
    stop("p must be positive")
  if(b <= 0)
    stop("b must be positive")
  return((rgamma(n, shape=a/p, rate=b))^(1/p))
}

dF2forInteg <- function(x, p) {
  ( exp(-x^p)- exp(-x) )/x
}

getk1 <- function(alpha,p){
  if(alpha < 0)
    stop("alpha must be nonnegative")
  if(p <= 1)
    stop("p must be greater than 1")
  if(alpha>0)
    return((exp(-1)-gamma_inc(1-alpha/p,1))/alpha)
  else return(gamma_inc(0,1)/p)
}

getk2 <- function(alpha,p){
  if(alpha < 0 | alpha >= 1)
    stop("a must be in the interval [0,1)")
  if(p <= 1)
    stop("p must be greater than 1")
  if(alpha>0)
    return( (gamma(1-alpha) - gamma_inc(1-alpha,1)-gamma(1-alpha/p)+gamma_inc(1-alpha/p,1))/alpha )
  else return( integrate(dF2forInteg,lower=0,upper=1,p=p)$value)
}

dF1 <- function(x, a, p) {
  if(a < 0)
    stop("a must be >= 0")
  if(p <= 1)
    stop("p must be greater than 1")
  ret = rep(0,length(x))
  ret[x>1] = exp(-x[x>1]^p)*x[x>1]^(-a-1)/getk1(a,p)

  return(ret)
}

rF1 <- function(n, a, p) {
  if(n < 1 || is.integer(n))
    stop("n must be an integer greater than or equal to 1")
  if(a < 0)
    stop("a must be >= 0")
  if(p <= 1)
    stop("p must be greater than 1")

  .C("rF1", as.integer(n), as.double(a), as.double(p), as.double(vector("double", n)))[[4]]
}

dF2 <- function(x,a, p) {
  if(a < 0 | a>=1)
    stop("a must be in the interval [0,1)")
  if(p <= 1)
    stop("p must be greater than 1")
  ret = rep(0,length(x))
  loc = (x>0)&(x<1)
  ret[loc] = ( exp(-x[loc]^p)- exp(-x[loc]) )*x[loc]^(-a-1)/getk2(a,p)

  return(ret)
}

rF2 <- function(n, a, p) {
  if(n < 1 || is.integer(n))
    stop("n must be an integer greater than or equal to 1")
  if(a < 0 | a>=1)
    stop("a must be in the interval [0,1)")
  if(p <= 1)
    stop("p must be greater than 1")

  .C("rF2", as.integer(n), as.double(a), as.double(p), as.double(vector("double", n)))[[4]]
}

rSubCTS <- function(n, alpha, c, ell, method=NULL) {
  if(n < 1 || is.integer(n))
    stop("n must be an integer greater than or equal to 1")
  if(alpha < 0 | alpha >= 1)
    stop("alpha must be in the interval [0,1)")
  if(c <= 0)
    stop("c must be positive")
  if(ell <= 0)
    stop("ell must be positive")

  if(alpha==0){return(rgamma(n,c,scale=ell))}
  else return(retstable(alpha, rep(c*gamma(1-alpha)/alpha, n), h=1/ell, method=method))
}

dSubCTS <- function(x, alpha, c, ell) {
  if(alpha < 0 | alpha >= 1)
    stop("alpha must be in the interval [0,1)")
  if(c <= 0)
    stop("c must be positive")
  if(ell <= 0)
    stop("ell must be positive")

  if(alpha==0){return(dgamma(x,c,scale=ell))}
  else{
    xi = (2-alpha)/(1-alpha)
    mu = c*(gamma(1-alpha)/alpha)*alpha*ell^(1-alpha)
    phi = (c*gamma(1-alpha)/alpha*alpha)^(1-xi)*(1-alpha)

    ret = rep(0,length(x))

    ret[x>=0] = dtweedie(x[x>=0],xi,mu,phi)
    return(ret)
  }
}

#optimized in simTandW
CC <- function(x,alpha){
  zeta=1/gamma(1-alpha)
  A0 = (1-alpha)*alpha^(alpha/(1-alpha))
  A0*exp( zeta^(1/alpha)*x^(1-1/alpha)*alpha*(1-alpha)^(1/alpha-1) )*(A0-x)^(alpha-2)
}

logCC <- function(x,alpha){ #returns log(CC)-log(A0)
  zeta=1/gamma(1-alpha)
  A0 = (1-alpha)*alpha^(alpha/(1-alpha))
  zeta^(1/alpha)*x^(1-1/alpha)*alpha*(1-alpha)^(1/alpha-1) + (alpha-2)*log(A0-x)
}

simTandW <-function(alpha){  #Implements Alg 4.1 in Dassios et al. (2020)

  if(alpha <= 0 | alpha >= 1)
    stop("alpha must be in the open interval (0,1)")

  a0 = (1-alpha)*alpha^(alpha/(1-alpha))
  temp = optimize(f=logCC,alpha, interval=c(0,a0))
  lambda = temp$minimum
  cMax = a0*exp(temp$objective)

 .C("simTandW", as.double(alpha), as.double(cMax), as.double(lambda), as.double(vector("double", 2)))[[4]]

}

simCondS <-function(t, alpha){ #Implements Alg 4.2 in Dassios et al. (2020)

  if(alpha <= 0 | alpha >= 1)
    stop("alpha must be in the open interval (0,1)")
  if(t <= 0)
    stop("t must be positive")

  .C("simCondS", as.double(t), as.double(alpha), as.double(0))[[3]]
}

rTrunS <- function(n, t, alpha, b=1, step=1){ #Implements Alg 4.3 in Dassios et al. (2020)

  if(alpha <= 0 | alpha >= 1)
    stop("alpha must be in the open interval (0,1)")
  if(t <= 0)
    stop("t must be positive")
  if(b <= 0)
    stop("b must be positive")
  if(step <= 0)
    stop("step must be positive")

  t = t*b^(-alpha)

  a0 = (1-alpha)*alpha^(alpha/(1-alpha))
  temp = optimize(f=logCC,alpha, interval=c(0,a0))
  lambda = temp$minimum
  cMax = a0*exp(temp$objective)

  .C("rTrunSOptim", as.integer(n), as.double(t), as.double(alpha), as.double(cMax), as.double(lambda), as.double(step), as.double(vector("double", n)))[[7]] * b

}

rTrunTS <- function(n, t, mu, alpha, b=1, step=1){ #Implements Alg. 4.4 in Dassios et al. (2020), mu>0
  
  if(alpha <= 0 | alpha >= 1)
    stop("alpha must be in the open interval (0,1)")
  if(t <= 0)
    stop("t must be positive")
  if(b <= 0)
    stop("b must be positive")
  if(mu <= 0)
    stop("mu must be positive")
  
  t = t*b^(-alpha)
  mu = mu*b
  
  a0 = (1-alpha)*alpha^(alpha/(1-alpha))
  temp = optimize(f=logCC,alpha, interval=c(0,a0))
  lambda = temp$minimum
  cMax = a0*exp(temp$objective)
  
  .C("rTrunTS", as.integer(n), as.double(t), as.double(alpha), as.double(cMax), as.double(lambda), as.double(mu), as.double(step), as.double(vector("double", n)))[[8]] * b
  
}


rPRDTS <- function(n, t, mu, alpha, p, step=1) {

  if(alpha >= 1)
    stop("alpha must be <1")
  if(t <= 0)
    stop("t must be positive")
  if(mu <= 0)
    stop("mu must be positive")
  if(step <=0)
    stop("step must be positive")

  if(alpha>0){

    if(p <= 1)
      stop("p must be > 1 when alpha>0")

    a0 = (1-alpha)*alpha^(alpha/(1-alpha))
    temp = optimize(f=logCC,alpha, interval=c(0,a0))
    lambda = temp$minimum
    cMax = a0*exp(temp$objective)
    t = t*mu^alpha*gamma(1-alpha)/alpha

    k1 = getk1(alpha,p)
    k2 = getk2(alpha,p)

    return(.C("rPRDTS", as.integer(n), as.double(t), as.double(alpha), as.double(cMax), as.double(lambda), as.double(p), as.double(k1), as.double(k2), as.double(step), as.double(vector("double", n)))[[10]]/mu)
  }
  else{
    if(alpha==0){
      if(p <= 1)
        stop("p must be > 1 when alpha=0")
      return(rPGamma(n, t, mu, p, step))
    }
    else{
      if(p <= 0)
        stop("p must be > 0 when alpha<0")
        t = t*mu^alpha
        return(.C("rPRDTSneg", as.integer(n), as.double(t), as.double(alpha), as.double(p), as.double(vector("double", n)))[[5]]/mu)
      }
  }
}

rDickman <- function(n, t, b=1) {
  if(t <= 0)
    stop("t must be positive")
  if(b <= 0)
    stop("b must be positive")

  .C("rDickman", as.integer(n), as.double(t), as.double(vector("double", n)))[[3]]*b

}

rTrunGamma <- function(n, t, mu, b=1, step=1) {
  if(t <= 0)
    stop("t must be positive")
  if(mu <= 0)
    stop("mu must be positive")
  if(b <= 0)
    stop("b must be positive")
  if(step <= 0)
    stop("step must be positive")
  
  mu = mu * b
  .C("rTrunGamma", as.integer(n), as.double(t), as.double(mu), as.double(step), as.double(vector("double", n)))[[5]]*b
}

rPGamma <- function(n, t, mu, p, step=1) {
  if(t <= 0)
    stop("t must be positive")
  if(p < 1)
    stop("p must be >= 1")
  if(mu <= 0)
    stop("mu must be positive")
  if(step <= 0)
    stop("step must be positive")
  
  if(p==1){
    return(rgamma(n,shape=t,rate=mu))
  }
  else{
    k1 = getk1(0,p)
    k2 = getk2(0,p)
    
    return(.C("rPGamma", as.integer(n), as.double(t), as.double(p), as.double(k1), as.double(k2), as.double(step), as.double(vector("double", n)))[[7]]/mu)
  }
}