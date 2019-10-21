#Setting up conditional hazard functions for bivariate data whose dependency is given by a Frank copula
#Marginals (CDFs and densities) must have been defined previously (via setup_marginals.R)
#Copula parameter a has to be defined before running

require(BivarP)

#Attention: x2 (2nd argument in the following functions) always refers to the event NOT named in the name of the function!
#The terms nominator and denominator refer to the nominator and denominator of the final hazard formulas in Chapter 5 of the paper 

jointcdf <- function(x1,x2){
  x1.cdf <- cdf1(x1)
  x2.cdf <- cdf2(x2)
  
  return( (-1/a) * log (1 + ( ( (exp(-a*x1.cdf) - 1) * (exp(-a*x2.cdf) - 1) ) / (exp(-a)-1) ) ) )
}

lambda.1a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.1a(x1,x1))
  else {
    
    x1.cdf <- cdf1(x1)
    x2.cdf <- cdf2(x2)
    x1.density <- density1(x1)
    
    nominator <- ( exp(a*x1.cdf) * (exp(a) - exp(a*x2.cdf)) * x1.density ) / ( -exp(a) + (exp(a+a*x1.cdf)) - exp(a*(x1.cdf+x2.cdf)) + exp(a+a*x2.cdf) )
    denominator <- 1 - x1.cdf - x2.cdf + jointcdf(x1,x2)
    return(nominator/denominator)
  }
}

lambda.1b <- function(x1,x2){
  
  x1.cdf <- cdf1(x1)
  x2.cdf <- cdf2(x2)
  x1.density <- density1(x1)
  x2.density <- density2(x2)
  
  nominator <- ( a * exp(a*(1+x1.cdf+x2.cdf)) * (exp(a)-1) * x1.density * x2.density ) / (( exp(a) - (exp(a+a*x2.cdf)) + exp(a*(x2.cdf+x1.cdf)) - exp(a+a*x1.cdf) )^2)
  denominator <- ( exp(a*x2.cdf) * (exp(a) - exp(a*x1.cdf)) * x2.density ) / ( -exp(a) + (exp(a+a*x2.cdf)) - exp(a*(x2.cdf+x1.cdf)) + exp(a+a*x1.cdf) )
  return(nominator/denominator)
}

lambda.2a <- function(x1,x2=NULL){
  if(is.null(x2)) return(lambda.2a(x1,x1))
  else {
    
    x1.cdf <- cdf2(x1)
    x2.cdf <- cdf1(x2)
    x1.density <- density2(x1)
    
    nominator <- ( exp(a*x1.cdf) * (exp(a) - exp(a*x2.cdf)) * x1.density ) / ( -exp(a) + (exp(a+a*x1.cdf)) - exp(a*(x1.cdf+x2.cdf)) + exp(a+a*x2.cdf) )
    denominator <- 1 - x1.cdf - x2.cdf + jointcdf(x1,x2)
    return(nominator/denominator)
  }
}

lambda.2b <- function(x1,x2){
  
  x1.cdf <- cdf2(x1)
  x2.cdf <- cdf1(x2)
  x1.density <- density2(x1)
  x2.density <- density1(x2)
  
  nominator <- ( a * exp(a*(1+x1.cdf+x2.cdf)) * (exp(a)-1) * x1.density * x2.density ) / (( exp(a) - (exp(a+a*x2.cdf)) + exp(a*(x2.cdf+x1.cdf)) - exp(a+a*x1.cdf) )^2)
  denominator <- ( exp(a*x2.cdf) * (exp(a) - exp(a*x1.cdf)) * x2.density ) / ( -exp(a) + (exp(a+a*x2.cdf)) - exp(a*(x2.cdf+x1.cdf)) + exp(a+a*x1.cdf) )
  return(nominator/denominator)
}

generate.event.times <- function(){
  raw.data <- BiCopGen(c(a1,b1,a2,b2,a), rodiny = c(marginal1,marginal2), rodina = "frank", n, cens = FALSE, bicens = FALSE, digi = 5)
  #nicht so ideal, lieber Events zum Zeitpunkt 0 anders abfangen!
  return(cbind(pmax(raw.data$X, 0.0001),pmax(raw.data$Y, 0.0001)))
}


