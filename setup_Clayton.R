#Setting up conditional hazard functions for bivariate data whose dependency is given by a Clayton copula
#Marginals (CDFs and densities) must have been defined previously (via setup_marginals.R)
#Copula parameter a has to be defined before running

require(BivarP)

#Attention: x2 (2nd argument in the following functions) always refers to the event NOT named in the name of the function!

lambda.1a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.1a(x1,x1))
  else {
    
    x1.cdf <- cdf1(x1)
    x2.cdf <- cdf2(x2)
    x1.density <- density1(x1)
    
    nominator <- x1.density * (1 - x1.cdf^(-1-a) * (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-1-(1/a)))
    denominator <- 1 - x1.cdf - x2.cdf + (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-1/a)
    return(nominator/denominator)
  }
}

lambda.1b <- function(x1,x2){
  
  x1.cdf <- cdf1(x1)
  x2.cdf <- cdf2(x2)
  x1.density <- density1(x1)
  
  nominator <- (a+1) * x1.density * (x1.cdf^(-1-a) * x2.cdf^(-1-a) * (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-2-(1/a)))
  denominator <- 1 - x2.cdf^(-1-a) * (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-1-(1/a))
  return(nominator/denominator)
}

lambda.2a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.2a(x1,x1))
  else{
    
    x1.cdf <- cdf2(x1)
    x2.cdf <- cdf1(x2)
    x1.density <- density2(x1)
    
    nominator <- x1.density * (1 - x1.cdf^(-1-a) * (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-1-(1/a)))
    denominator <- 1 - x1.cdf - x2.cdf + (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-1/a)
    return(nominator/denominator)
  }
}

lambda.2b <- function(x1,x2){
  
  x1.cdf <- cdf2(x1)
  x2.cdf <- cdf1(x2)
  x1.density <- density2(x1)
  
  nominator <- (a+1) * x1.density * (x1.cdf^(-1-a) * x2.cdf^(-1-a) * (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-2-(1/a)))
  denominator <- 1 - x2.cdf^(-1-a) * (x1.cdf^(-a) + x2.cdf^(-a) - 1)^(-1-(1/a))
  return(nominator/denominator)
}

generate.event.times <- function(){
  raw.data <- BiCopGen(c(a1,b1,a2,b2,a), rodiny = c(marginal1,marginal2), rodina = "clayton", n, cens = FALSE, bicens = FALSE, digi = 5)
  #nicht so ideal, lieber Events zum Zeitpunkt 0 anders abfangen!
  return(cbind(pmax(raw.data$X, 0.0001),pmax(raw.data$Y, 0.0001)))
}
