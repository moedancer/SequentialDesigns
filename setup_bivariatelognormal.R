#Setting up conditional hazard functions for bivariate log-normally distributed data 

require(MASS)

mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
rho <- 0

lambda.1a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.1a(x1,x1))
  else {
    x1.rescale <- (mu1-log(x1))/sigma1
    x2.rescale <- (mu2-log(x2))/sigma2
    nominator <- dnorm(x1.rescale)*pnorm((x2.rescale - rho*x1.rescale)/sqrt(1-rho^2))
    denominator <- sigma1 * x1 * pmvnorm(mean = c(0,0), corr = matrix(c(1,rho,rho,1),nrow = 2), lower = c(-Inf,-Inf), upper = c(x1.rescale,x2.rescale))
    return(nominator/denominator)
  }
}

lambda.1b <- function(x1,x2){
  nominator <- dnorm( (log(x1) - mu1 - (sigma1/sigma2) * rho * (log(x2)-mu2)) / (sigma1 * sqrt(1-rho^2)) )
  denominator <- sqrt(1-rho^2) * sigma1 * x1 * (1 - pnorm( (log(x1) - mu1 - (sigma1/sigma2) * rho * (log(x2)-mu2)) / (sigma1 * sqrt(1-rho^2)) ) )
  return(nominator/denominator)
}

#ATTENTION: If only one parameter is given, it shall logically be interpreted as the time parameter concerning the second event
lambda.2a <- function(x1,x2=NULL){
  if(is.null(x1)) return(lambda.2a(x1,x1))
  else{
    x1.rescale <- (mu1-log(x1))/sigma1
    x2.rescale <- (mu2-log(x2))/sigma2
    nominator <- dnorm(x2.rescale)*pnorm((x1.rescale - rho*x2.rescale)/sqrt(1-rho^1))
    denominator <- sigma2 * x2 * pmvnorm(mean = c(0,0), corr = matrix(c(1,rho,rho,1),nrow = 2), lower = c(-Inf,-Inf), upper = c(x1.rescale,x2.rescale))
    return(nominator/denominator)
  }
}

lambda.2b <- function(x1,x2){
  nominator <- dnorm( (log(x2) - mu2 - (sigma2/sigma1) * rho * (log(x1)-mu1)) / (sigma2 * sqrt(1-rho^2)) )
  denominator <- sqrt(1-rho^2) * sigma2 * x2 * (1 - pnorm( (log(x2) - mu1 - (sigma2/sigma1) * rho * (log(x1)-mu1)) / (sigma2 * sqrt(1-rho^2)) ) )
  return(nominator/denominator)
}

generate.event.times <- function(){
  cov.mat <- matrix(c(sigma1^2, sigma1*sigma2*rho, sigma1*sigma2*rho, sigma2^2), ncol = 2)
  return(raw.data <- exp(mvrnorm(n = n, mu = c(mu1, mu2), Sigma = cov.mat)))
}
