#Setting up CDFs and densities for both marginal distributions
#All parameters (a1, a2, b1, b2) have to be defined before running

require(stats)

#set up densities and cdf of marginal distributions according to the user's choice

if(marginal1 == "weibull"){
  density1 <- function(x) dweibull(x, shape = a1, scale = b1)
  cdf1 <- function(x) pweibull(x, shape = a1, scale = b1)
} else if(marginal1 == "gamma"){
  density1 <- function(x) dgamma(x, shape = a1, scale = b1)
  cdf1 <- function(x) pgamma(x, shape = a1, scale = b1)
} else if(marginal1 == "lnorm"){
  density1 <- function(x) dlnorm(x, shape = a1, scale = b1)
  cdf1 <- function(x) plnorm(x, shape = a1, scale = b1)
}

if(marginal2 == "weibull"){
  density2 <- function(x) dweibull(x, shape = a2, scale = b2)
  cdf2 <- function(x) pweibull(x, shape = a2, scale = b2)
} else if(marginal2 == "gamma"){
  density2 <- function(x) dgamma(x, shape = a2, scale = b2)
  cdf2 <- function(x) pgamma(x, shape = a2, scale = b2)
} else if(marginal2 == "lnorm"){
  density2<- function(x) dlnorm(x, shape = a2, scale = b2)
  cdf2 <- function(x) plnorm(x, shape = a2, scale = b2)
}

