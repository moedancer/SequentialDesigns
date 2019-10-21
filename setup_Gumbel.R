#Setting up conditional hazard functions for bivariate data whose dependency is given by a Gumbel copula
#Marginals (CDFs and densities) must have been defined previously (via setup_marginals.R)
#Copula parameter a has to be defined before running

require(BivarP)

#Attention: x2 (2nd argument in the following functions) always refers to the event NOT named in the name of the function!

jointcdf <- function(x1,x2){
  x1.cdf <- cdf1(x1)
  x2.cdf <- cdf2(x2)
  
  return( exp ( - ( (-log(x1.cdf))^a + (-log(x2.cdf))^a )^(1/a) ) )
}

lambda.1a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.1a(x1,x1))
  else {
    
    x1.cdf <- cdf1(x1)
    x2.cdf <- cdf2(x2)
    x1.density <- density1(x1)
    
    nominator <- x1.density * (1 - jointcdf(x1,x2) * ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^((1/a)-1) * (-log(x1.cdf))^(a-1) * (1/x1.cdf))
    denominator <- 1 - x1.cdf - x2.cdf + jointcdf(x1,x2)
    return(nominator/denominator)
  }
}

lambda.1b <- function(x1,x2){
  
  x1.cdf <- cdf1(x1)
  x2.cdf <- cdf2(x2)
  x1.density <- density1(x1)
  x2.density <- density2(x2)
  
  nominator <- ((x1.density * x2.density)/(x1.cdf * x2.cdf)) * jointcdf(x1,x2) * (-log(x1.cdf))^(a-1) * (-1 + a + ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^(1/a)) * ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^(-2+(1/a)) * (-log(x2.cdf))^(a-1)
  denominator <- x2.density * (1 - jointcdf(x1,x2) * ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^((1/a)-1) * (-log(x2.cdf))^(a-1) * (1/x2.cdf))
  return(nominator/denominator)
}

lambda.2a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.2a(x1,x1))
  else {
    
    x1.cdf <- cdf2(x1)
    x2.cdf <- cdf1(x2)
    x1.density <- density2(x1)
    
    nominator <- x1.density * (1 - jointcdf(x2,x1) * ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^((1/a)-1) * (-log(x1.cdf))^(a-1) * (1/x1.cdf))
    denominator <- 1 - x1.cdf - x2.cdf + jointcdf(x2,x1)
    return(nominator/denominator)
  }
}

lambda.2b <- function(x1,x2){
  
  x1.cdf <- cdf2(x1)
  x2.cdf <- cdf1(x2)
  x1.density <- density2(x1)
  x2.density <- density1(x2)
  
  nominator <- ((x1.density * x2.density)/(x1.cdf * x2.cdf)) * jointcdf(x2,x1) * (-log(x1.cdf))^(a-1) * (-1 + a + ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^(1/a)) * ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^(-2+(1/a)) * (-log(x2.cdf))^(a-1)
  denominator <- x2.density * (1 - jointcdf(x2,x1) * ((-log(x1.cdf))^a + (-log(x2.cdf))^a)^((1/a)-1) * (-log(x2.cdf))^(a-1) * (1/x2.cdf))
  return(nominator/denominator)
}

generate.event.times <- function(){
  raw.data <- BiCopGen(c(a1,b1,a2,b2,a), rodiny = c(marginal1,marginal2), rodina = "gumbel", n, cens = FALSE, bicens = FALSE, digi = 5)
  #nicht so ideal, lieber Events zum Zeitpunkt 0 anders abfangen!
  return(cbind(pmax(raw.data$X, 0.0001),pmax(raw.data$Y, 0.0001)))
}


