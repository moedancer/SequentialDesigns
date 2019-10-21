#Simulation to compare results for different choices of time point for interim analyses
#"setup_[...].R" scripts will be sourced according to the chosen distribution 

require(Hmisc) #function: rcorr
#setwd() according to folder containing scripts

#####################################################################################################################

#set up distribution for simulation

#Simulate bivariate log-normal distribution? TRUE or FALSE
BivariateLogNormal <- FALSE
#if TRUE is chosen, set the parameters for the bivariate log-normal distribution
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
rho <- 0

#In case of no bivariate log-normal distribution: Which copula shall be used? "Clayton", "Gumbel" or "Frank"?
Sim.Copula <- "Frank"
#Choose copula parameter for chosen copula 
a <- 20
#Choose marginal distributions! Choose among "weibull", "gamma" and "lnorm"
marginal1 <- "weibull"
marginal2 <- "weibull"
#Choose parameters for marginal distributions: a1 and b1 refer to 1st marginal, a2 and b2 refer to 2nd marginal
  #In case of Weibull or Gamma-distribution: a is shape and b is scale parameter
  #In case of log-normal distribution: a is meanlog and b is sdlog
m1 <- 2
a1 <- 1
b1 <- m1 * log(2)^(-1/a1)

m2 <- 1
a2 <- 1
b2 <- m2 * log(2)^(-1/a2)

#####################################################################################################################

#do the setup according to the chosen distribution

if(BivariateLogNormal == TRUE){
  source("setup_bivariatelognormal.R")
} else {
  source("setup_marginals.R")
  source(paste("setup_", Sim.Copula, ".r", sep = ""))
}

#####################################################################################################################

t1 <- 2
t1_alt <- 2.5
t2 <- 5
accrual.period <- 3
followup <- 2

simulations <- 100000
resultscol <- 12
results <- data.frame(matrix(rep(0,simulations * resultscol),ncol = resultscol))
colnames(results) <- c("Event1.t1", "Var.Event1.t1", "Event2.t1", "Var.Event2.t1", "Event1.t1_alt", "Var.Event1.t1_alt", "Event2.t1_alt", "Var.Event2.t1_alt", "Event1.t2", "Var.Event1.t2", "Event2.t2", "Var.Event2.t2")

n <- 100

for(i in 1:simulations){
  
  tryCatch({
    
    arrival.times <- runif(n, min = 0, max = accrual.period)
    censoring.times <- rep(followup, n)
    event.times <- generate.event.times()
    trial.data <- data.frame(arrival.times, event.times, event.times + arrival.times, censoring.times + arrival.times)
    colnames(trial.data) <- c("Arrival", "Event1", "Event2", "Event1Calendar", "Event2Calendar", "Censoring")
    
    #Compute bounds for integrals in calculation of compensators
    
    trial.data$changestate.t1 <- pmin(trial.data$Event1, trial.data$Event2, followup, pmax(0,(t1 - trial.data$Arrival)))
    trial.data$upperbound.comp1.t1 <- pmin(trial.data$Event1, followup, pmax(0,(t1 - trial.data$Arrival)))
    trial.data$upperbound.comp2.t1 <- pmin(trial.data$Event2, followup, pmax(0,(t1 - trial.data$Arrival)))
    
    trial.data$changestate.t1_alt <- pmin(trial.data$Event1, trial.data$Event2, followup, pmax(0,(t1_alt - trial.data$Arrival)))
    trial.data$upperbound.comp1.t1_alt <- pmin(trial.data$Event1, followup, pmax(0,(t1_alt - trial.data$Arrival)))
    trial.data$upperbound.comp2.t1_alt <- pmin(trial.data$Event2, followup, pmax(0,(t1_alt - trial.data$Arrival)))
    
    trial.data$changestate.t2 <- pmin(trial.data$Event1, trial.data$Event2, followup, pmax(0,(t2 - trial.data$Arrival)))
    trial.data$upperbound.comp1.t2 <- pmin(trial.data$Event1, followup, pmax(0,(t2 - trial.data$Arrival)))
    trial.data$upperbound.comp2.t2 <- pmin(trial.data$Event2, followup, pmax(0,(t2 - trial.data$Arrival)))
    
    #calculate compensators
    
    trial.data$compensator.comp1.initialstate.t1 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.1a, u, v)), 0, trial.data$changestate.t1))
    trial.data$compensator.comp1.changedstate.t1 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.1b, u, v, x2 = w)), trial.data$changestate.t1, trial.data$upperbound.comp1.t1, trial.data$Event2))
    trial.data$compensator.comp2.initialstate.t1 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.2a, u, v)), 0, trial.data$changestate.t1))
    trial.data$compensator.comp2.changedstate.t1 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.2b, u, v, x2 = w)), trial.data$changestate.t1, trial.data$upperbound.comp2.t1, trial.data$Event1))
    
    trial.data$compensator.comp1.initialstate.t1_alt <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.1a, u, v)), 0, trial.data$changestate.t1_alt))
    trial.data$compensator.comp1.changedstate.t1_alt <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.1b, u, v, x2 = w)), trial.data$changestate.t1_alt, trial.data$upperbound.comp1.t1_alt, trial.data$Event2))
    trial.data$compensator.comp2.initialstate.t1_alt <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.2a, u, v)), 0, trial.data$changestate.t1_alt))
    trial.data$compensator.comp2.changedstate.t1_alt <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.2b, u, v, x2 = w)), trial.data$changestate.t1_alt, trial.data$upperbound.comp2.t1_alt, trial.data$Event1))
    
    trial.data$compensator.comp1.initialstate.t2 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.1a, u, v)), 0, trial.data$changestate.t2))
    trial.data$compensator.comp1.changedstate.t2 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.1b, u, v, x2 = w)), trial.data$changestate.t2, trial.data$upperbound.comp1.t2, trial.data$Event2))
    trial.data$compensator.comp2.initialstate.t2 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.2a, u, v)), 0, trial.data$changestate.t2))
    trial.data$compensator.comp2.changedstate.t2 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.2b, u, v, x2 = w)), trial.data$changestate.t2, trial.data$upperbound.comp2.t2, trial.data$Event1))
    
    #calculate martingales on patient level
    trial.data$martingale.comp1.t1 <- ifelse(trial.data$Event1 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp1.initialstate.t1 - trial.data$compensator.comp1.changedstate.t1
    trial.data$martingale.comp2.t1 <- ifelse(trial.data$Event2 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp2.initialstate.t1 - trial.data$compensator.comp2.changedstate.t1 
    
    trial.data$martingale.comp1.t1_alt <- ifelse(trial.data$Event1 <= pmin(followup, pmax(0,t1_alt - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp1.initialstate.t1_alt - trial.data$compensator.comp1.changedstate.t1_alt
    trial.data$martingale.comp2.t1_alt <- ifelse(trial.data$Event2 <= pmin(followup, pmax(0,t1_alt - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp2.initialstate.t1_alt - trial.data$compensator.comp2.changedstate.t1_alt 
    
    trial.data$martingale.comp1.t2 <- ifelse(trial.data$Event1 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp1.initialstate.t2 - trial.data$compensator.comp1.changedstate.t2
    trial.data$martingale.comp2.t2 <- ifelse(trial.data$Event2 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp2.initialstate.t2 - trial.data$compensator.comp2.changedstate.t2
    
    #calculate overall martingales
    results$Event1.t1[i] <- n^(-1/2) * sum(trial.data$martingale.comp1.t1)
    results$Event2.t1[i] <- n^(-1/2) * sum(trial.data$martingale.comp2.t1)
    results$Event1.t1_alt[i] <- n^(-1/2) * sum(trial.data$martingale.comp1.t1_alt)
    results$Event2.t1_alt[i] <- n^(-1/2) * sum(trial.data$martingale.comp2.t1_alt)
    results$Event1.t2[i] <- n^(-1/2) * sum(trial.data$martingale.comp1.t2)
    results$Event2.t2[i] <- n^(-1/2) * sum(trial.data$martingale.comp2.t2)
    
    results$Var.Event1.t1[i] <- sum(trial.data$Event1 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)))
    results$Var.Event2.t1[i] <- sum(trial.data$Event2 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)))
    results$Var.Event1.t1_alt[i] <- sum(trial.data$Event1 <= pmin(followup, pmax(0,t1_alt - trial.data$Arrival)))
    results$Var.Event2.t1_alt[i] <- sum(trial.data$Event2 <= pmin(followup, pmax(0,t1_alt - trial.data$Arrival)))
    results$Var.Event1.t2[i] <- sum(trial.data$Event1 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)))
    results$Var.Event2.t2[i] <- sum(trial.data$Event2 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)))
    
  }, error = function(e){cat("ERROR in simulation ", i, ": ", conditionMessage(e), sep = "", "\n")})
  
}

#Calculate final statistics

results$Z1.t1 <- sqrt(n) * results$Event1.t1/sqrt(results$Var.Event1.t1)
results$Z2.t1 <- sqrt(n) * results$Event2.t1/sqrt(results$Var.Event2.t1)
results$Z1.t1_alt <- sqrt(n) * results$Event1.t1_alt/sqrt(results$Var.Event1.t1_alt)
results$Z2.t1_alt <- sqrt(n) * results$Event2.t1_alt/sqrt(results$Var.Event2.t1_alt)
results$Z1.t2 <- sqrt(n) * (results$Event1.t2 - results$Event1.t1)/sqrt(results$Var.Event1.t2 - results$Var.Event1.t1)
results$Z2.t2 <- sqrt(n) * (results$Event2.t2 - results$Event2.t1)/sqrt(results$Var.Event2.t2 - results$Var.Event2.t1)
results$Z1.t2_alt <- sqrt(n) * (results$Event1.t2 - results$Event1.t1_alt)/sqrt(results$Var.Event1.t2 - results$Var.Event1.t1_alt)
results$Z2.t2_alt <- sqrt(n) * (results$Event2.t2 - results$Event2.t1_alt)/sqrt(results$Var.Event2.t2 - results$Var.Event2.t1_alt)

#####################################################################################################################
#Evaluation of empirical distribution of test-statistics

#compute empirical correlations
corr_matrix <- diag(4)
test.statistics <- c("Z1.t1_alt", "Z2.t1_alt", "Z1.t2_alt", "Z2.t2_alt")
for(i in 1:length(test.statistics)){
  for(j in 1:length(test.statistics)){
    corr_matrix[i,j] = rcorr(results[,test.statistics[i]], results[,test.statistics[j]])$r[2]
  }
}

corr_matrix <- round(corr_matrix, digits = 4)
colnames(corr_matrix) <- c(expression(Z[1]^{(1)}), expression(Z[1]^{(2)}), expression(Z[2]^{(1)}), expression(Z[2]^{(2)}))
rownames(corr_matrix) <- c(expression(Z[1]^{(1)}), expression(Z[1]^{(2)}), expression(Z[2]^{(1)}), expression(Z[2]^{(2)}))


#compute empirical error level for L1 and L2 rejection regions
alpha <- 0.05

u <- qnorm( 1 - ( 1 - (1-alpha)^(1/4) ) / 2)
results$rejection.L1<- 1-(abs(results$Z1.t1)<u)*(abs(results$Z2.t1)<u)*(abs(results$Z1.t2)<u)*(abs(results$Z2.t2)<u)
results$rejection.L1_alt<- 1-(abs(results$Z1.t1_alt)<u)*(abs(results$Z2.t1_alt)<u)*(abs(results$Z1.t2_alt)<u)*(abs(results$Z2.t2_alt)<u)
print("Maximum norm, Pocock bounds, t1 = 2:")
sum(results$rejection.L1)/simulations
print("Maximum norm, Pocock bounds, t1 = 2.5:")
sum(results$rejection.L1_alt)/simulations

u <- qchisq( 1 - ( 1 - (1-alpha)^(1/2) ), df = 2)
results$rejection.L2 <- 1-((results$Z1.t1)^2 + (results$Z2.t1)^2 < u)*((results$Z1.t2)^2 + (results$Z2.t2)^2 < u)
results$rejection.L2_alt <- 1-((results$Z1.t1_alt)^2 + (results$Z2.t1_alt)^2 < u)*((results$Z1.t2_alt)^2 + (results$Z2.t2_alt)^2 < u)
print("L2 norm, Pocock bounds, t1 = 2:")
sum(results$rejection.L2)/simulations
print("L2 norm, Pocock bounds, t1 = 2.5:")
sum(results$rejection.L2_alt)/simulations

weight.factor <- 0.5

alpha_star <- ((1+weight.factor)/(2*weight.factor)) - sqrt(((1+weight.factor)/(2*weight.factor))^2 - (1/weight.factor) * (1 - sqrt(1-alpha)))
u1 <- qnorm(1 - weight.factor*alpha_star/2)
u2 <- qnorm(1 - alpha_star/2)
results$rejection.weightedbounds.L1<- 1-(abs(results$Z1.t1)<u1)*(abs(results$Z2.t1)<u1)*(abs(results$Z1.t2)<u2)*(abs(results$Z2.t2)<u2)
results$rejection.weightedbounds.L1_alt<- 1-(abs(results$Z1.t1_alt)<u1)*(abs(results$Z2.t1_alt)<u1)*(abs(results$Z1.t2_alt)<u2)*(abs(results$Z2.t2_alt)<u2)
print("Maximum norm, weighted bounds, t1 = 2:")
sum(results$rejection.weightedbounds.L1)/simulations
print("Maximum norm, weighted bounds, t1 = 2.5:")
sum(results$rejection.weightedbounds.L1_alt)/simulations

alpha_star <- ((1+weight.factor)/(2*weight.factor)) - sqrt(((1+weight.factor)/(2*weight.factor))^2 - (1/weight.factor) * alpha)
u1 <- qchisq(1 - weight.factor*alpha_star, df = 2)
u2 <- qchisq(1 - alpha_star, df = 2)
results$rejection.weightedbounds.L2 <- 1-((results$Z1.t1)^2 + (results$Z2.t1)^2 < u1)*((results$Z1.t2)^2 + (results$Z2.t2)^2 < u2)
results$rejection.weightedbounds.L2_alt <- 1-((results$Z1.t1_alt)^2 + (results$Z2.t1_alt)^2 < u1)*((results$Z1.t2_alt)^2 + (results$Z2.t2_alt)^2 < u2)
print("L2 norm, weighted bounds, t1 = 2:")
sum(results$rejection.weightedbounds.L2)/simulations
print("L2 norm, weighted bounds, t1 = 2.5:")
sum(results$rejection.weightedbounds.L2_alt)/simulations


#show histograms with standard normal density
par(mfrow = c(2,2), pty = "m")

hist(results$Z1.t1, freq = FALSE,
     main = expression( paste("Marginal distribution of ", Z[1]^{(1)} )),
     xlim = c(-4,4),
     xlab = expression( Z[1]^{(1)} ),
     breaks = seq(floor(min(results$Z1.t1)), ceiling(max(results$Z1.t1)), by = 0.1))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

hist(results$Z2.t1, freq = FALSE,
     main = expression( paste("Marginal distribution of ", Z[1]^{(2)} )),
     xlim = c(-4,4),
     xlab = expression( Z[1]^{(2)} ),
     breaks = seq(floor(min(results$Z2.t1)), ceiling(max(results$Z2.t1)), by = 0.1))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

hist(results$Z1.t2, freq = FALSE,
     main = expression( paste("Marginal distribution of ", Z[2]^{(1)} )),
     xlim = c(-4,4),
     xlab = expression( Z[2]^{(1)} ),
     breaks = seq(floor(min(results$Z1.t2)), ceiling(max(results$Z1.t2)), by = 0.1))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

hist(results$Z2.t2, freq = FALSE,
     main = expression( paste("Marginal distribution of ", Z[2]^{(2)} )),
     xlim = c(-4,4),
     xlab = expression( Z[2]^{(2)} ),
     breaks = seq(floor(min(results$Z2.t2)), ceiling(max(results$Z2.t2)), by = 0.1))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")


#show histograms for alternative choice of interim analysis with standard normal density
par(mfrow = c(2,2), pty = "m")

hist(results$Z1.t1_alt, freq = FALSE,
     main = expression( paste("Empirical distribution of ", Z[1]^{(1)})),
     xlim = c(-4,4),
     xlab = expression( Z[1]^{(1)} ),
     breaks = seq(floor(min(results$Z1.t1_alt)), ceiling(max(results$Z1.t1_alt)), by = 0.2))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

hist(results$Z2.t1_alt, freq = FALSE,
     main = expression( paste("Empirical distribution of ", Z[1]^{(2)})),
     xlim = c(-4,4),
     xlab = expression( Z[1]^{(2)} ),
     breaks = seq(floor(min(results$Z2.t1_alt)), ceiling(max(results$Z2.t1_alt)), by = 0.2))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

hist(results$Z1.t2_alt, freq = FALSE,
     main = expression( paste("Empirical distribution of ", Z[2]^{(1)} )),
     xlim = c(-4,4),
     xlab = expression( Z[2]^{(1)} ),
     breaks = seq(floor(min(results$Z1.t2_alt)), ceiling(max(results$Z1.t2_alt)), by = 0.2))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

hist(results$Z2.t2_alt, freq = FALSE,
     main = expression( paste("Empirical distribution of ", Z[2]^{(2)} )),
     xlim = c(-4,4),
     xlab = expression( Z[2]^{(2)} ),
     breaks = seq(floor(min(results$Z2.t2_alt)), ceiling(max(results$Z2.t2_alt)), by = 0.2))
curve(dnorm(x, mean=0, sd=1), col="black", lwd=2, add=TRUE, yaxt="n")

#show scatterplots for illustrating covariance and overview table for covariances
par(mfrow = c(2,2), pty = "s")

plot(results$Z2.t1_alt~results$Z1.t1_alt,
     xlim=c(-5,5), ylim=c(-5,5), 
     xlab=expression( Z[1]^{(1)}), ylab=expression( Z[1]^{(2)}), 
     main = expression( paste( "Joint empirical distribution of ", Z[1]^{(1)}, " and ", Z[1]^{(2)})),
     pch = 3, cex = 0.2)

plot(results$Z1.t2_alt~results$Z1.t1_alt,
     xlim=c(-5,5), ylim=c(-5,5), 
     xlab=expression( Z[1]^{(1)}), ylab=expression( Z[2]^{(1)}), 
     main = expression( paste( "Joint empirical distribution of ", Z[1]^{(1)}, " and ", Z[2]^{(1)})),
     pch = 3, cex = 0.2)

plot(results$Z2.t2_alt~results$Z1.t1_alt,
     xlim=c(-5,5), ylim=c(-5,5), 
     xlab=expression( Z[1]^{(1)}), ylab=expression( Z[2]^{(2)}), 
     main = expression( paste( "Joint empirical distribution of ", Z[1]^{(1)}, " and ", Z[2]^{(2)})),
     pch = 3, cex = 0.2)

frame()
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)), rowhead = list(fg_params = list(parse = TRUE)))
grid.table(corr_matrix, theme = tt)
