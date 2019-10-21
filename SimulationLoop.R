#Execute simulations for all scnearios to be considered in the mauscript
#"setup_[...].R" scripts will be sourced in the loop according to the currently investigated distribution 

require(Hmisc) #function: rcorr
#setwd() according to folder containing scripts
set.seed(1)
#####################################################################################################################

#set up loop for simulation with different parameters
copulas.vector <- c("Clayton", "Gumbel", "Frank")
intensities <- c("low", "medium", "high")
num.copulas <- length(copulas.vector)
copula.parameters <- matrix(c(1, 5, 15, 1.5, 4, 8, 5, 20, 50), ncol = num.copulas)

colnames(copula.parameters) <- copulas.vector
rownames(copula.parameters) <- intensities

fixed.median <- 2
fixed.shape <- 1

alternative.median <- 1
alternative.shapes <- c(1/2, 1, 2)

sample.sizes <- c(50, 100, 250, 500)

emp.alpha.level <- matrix(rep(0, num.copulas * length(intensities) * length(alternative.median) * length(alternative.shapes) * length(sample.sizes) * 13), ncol=13)
colnames(emp.alpha.level) <- c("Sample size", "m1", "k1", "m2", "k2", "Copula", "Copula parameter", "alpha.maxnorm", "alpha.L2", "alpha.maxnorm.weight1", "alpha.L2.weight1", "alpha.maxnorm.weight2", "alpha.L2.weight2")

#Calculate rejection levels and rejection bounds for standard and weighted cases
alpha <- 0.05
weight1 <- 1/2
weight2 <- 1/4

u.maxnorm <- qnorm( 1 - ( 1 - (1-alpha)^(1/4) ) / 2)
u.L2 <- qchisq( 1 - ( 1 - (1-alpha)^(1/2) ), df = 2)

alpha_star.maxnorm.weight1 <- ((1+weight1)/(2*weight1)) - sqrt(((1+weight1)/(2*weight1))^2 - (1/weight1) * (1 - sqrt(1-alpha)))
u1.maxnorm.weight1 <- qnorm(1 - weight1*alpha_star.maxnorm.weight1/2)
u2.maxnorm.weight1 <- qnorm(1 - alpha_star.maxnorm.weight1/2)

alpha_star.L2.weight1 <- ((1+weight1)/(2*weight1)) - sqrt(((1+weight1)/(2*weight1))^2 - (1/weight1) * alpha)
u1.L2.weight1 <- qchisq(1 - weight1*alpha_star.L2.weight1, df = 2)
u2.L2.weight1 <- qchisq(1 - alpha_star.L2.weight1, df = 2)

alpha_star.maxnorm.weight2 <- ((1+weight2)/(2*weight2)) - sqrt(((1+weight2)/(2*weight2))^2 - (1/weight2) * (1 - sqrt(1-alpha)))
u1.maxnorm.weight2 <- qnorm(1 - weight2*alpha_star.maxnorm.weight2/2)
u2.maxnorm.weight2 <- qnorm(1 - alpha_star.maxnorm.weight2/2)

alpha_star.L2.weight2 <- ((1+weight2)/(2*weight2)) - sqrt(((1+weight2)/(2*weight2))^2 - (1/weight2) * alpha)
u1.L2.weight2 <- qchisq(1 - weight2*alpha_star.L2.weight2, df = 2)
u2.L2.weight2 <- qchisq(1 - alpha_star.L2.weight2, df = 2)

current.row.to.fill <- 0

for(sample.size.temp in sample.sizes){
  for(med2.temp in alternative.median){
    for(kappa2.temp in alternative.shapes){
      for(copula.temp in copulas.vector){
        for(intensity.temp in intensities){
      
          current.row.to.fill <- current.row.to.fill + 1
          print(current.row.to.fill)
          
          #set up distribution for simulation
          
          Sim.Copula <- copula.temp
          #Choose copula parameter for chosen copula 
          a <- copula.parameters[intensity.temp, copula.temp]
          #Choose marginal distributions! Choose among "weibull", "gamma" and "lnorm"
          marginal1 <- "weibull"
          marginal2 <- "weibull"
          #Choose parameters for marginal distributions: a1 and b1 refer to 1st marginal, a2 and b2 refer to 2nd marginal
            #In case of Weibull or Gamma-distribution: a is shape and b is scale parameter
            #In case of log-normal distribution: a is meanlog and b is sdlog
          m1 <- fixed.median
          a1 <- fixed.shape
          b1 <- m1 * log(2)^(-1/a1)
          
          m2 <- med2.temp
          a2 <- kappa2.temp
          b2 <- m2 * log(2)^(-1/a2)
          
          #####################################################################################################################
          
          #do the setup according to the chosen distribution
          
          source("setup_marginals.R")
          source(paste("setup_", Sim.Copula, ".r", sep = ""))
          
          #####################################################################################################################
          
          t1 <- 2.5
          t2 <- 5
          accrual.period <- 3
          followup <- 2
          
          simulations <- 10000
          resultscol <- 12 #number of results needed for each simulation step (2 events; 2 timesteps; martingale and variation; 4 columns for final statistics)
          results <- data.frame(matrix(rep(0,simulations * resultscol),ncol = resultscol))
          colnames(results) <- c("Event1.t1", "Var.Event1.t1", "Event2.t1", "Var.Event2.t1", "Event1.t2", "Var.Event1.t2", "Event2.t2", "Var.Event2.t2", "Z1.t1", "Z2.t1", "Z1.t2", "Z2.t2")
          
          n <- sample.size.temp
          
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
                
              trial.data$changestate.t2 <- pmin(trial.data$Event1, trial.data$Event2, followup, pmax(0,(t2 - trial.data$Arrival)))
              trial.data$upperbound.comp1.t2 <- pmin(trial.data$Event1, followup, pmax(0,(t2 - trial.data$Arrival)))
              trial.data$upperbound.comp2.t2 <- pmin(trial.data$Event2, followup, pmax(0,(t2 - trial.data$Arrival)))
              
              #calculate compensators
              
              trial.data$compensator.comp1.initialstate.t1 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.1a, u, v, rel.tol = 1e-10)), 0, trial.data$changestate.t1))
              trial.data$compensator.comp1.changedstate.t1 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.1b, u, v, x2 = w)), trial.data$changestate.t1, trial.data$upperbound.comp1.t1, trial.data$Event2))
              trial.data$compensator.comp2.initialstate.t1 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.2a, u, v)), 0, trial.data$changestate.t1))
              trial.data$compensator.comp2.changedstate.t1 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.2b, u, v, x2 = w)), trial.data$changestate.t1, trial.data$upperbound.comp2.t1, trial.data$Event1))
              
              trial.data$compensator.comp1.initialstate.t2 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.1a, u, v, rel.tol = 1e-10)), 0, trial.data$changestate.t2))
              trial.data$compensator.comp1.changedstate.t2 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.1b, u, v, x2 = w)), trial.data$changestate.t2, trial.data$upperbound.comp1.t2, trial.data$Event2))
              trial.data$compensator.comp2.initialstate.t2 <- unlist(mapply(function(u,v) ifelse(u==v, 0, integrate(lambda.2a, u, v)), 0, trial.data$changestate.t2))
              trial.data$compensator.comp2.changedstate.t2 <- unlist(mapply(function(u,v,w) ifelse(u==v, 0, integrate(lambda.2b, u, v, x2 = w)), trial.data$changestate.t2, trial.data$upperbound.comp2.t2, trial.data$Event1))
              
              #calculate martingales on patient level
              trial.data$martingale.comp1.t1 <- ifelse(trial.data$Event1 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp1.initialstate.t1 - trial.data$compensator.comp1.changedstate.t1
              trial.data$martingale.comp2.t1 <- ifelse(trial.data$Event2 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp2.initialstate.t1 - trial.data$compensator.comp2.changedstate.t1 
              
              trial.data$martingale.comp1.t2 <- ifelse(trial.data$Event1 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp1.initialstate.t2 - trial.data$compensator.comp1.changedstate.t2
              trial.data$martingale.comp2.t2 <- ifelse(trial.data$Event2 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)), 1, 0) - trial.data$compensator.comp2.initialstate.t2 - trial.data$compensator.comp2.changedstate.t2
              
              #calculate overall martingales
              results$Event1.t1[i] <- n^(-1/2) * sum(trial.data$martingale.comp1.t1)
              results$Event2.t1[i] <- n^(-1/2) * sum(trial.data$martingale.comp2.t1)
              results$Event1.t2[i] <- n^(-1/2) * sum(trial.data$martingale.comp1.t2)
              results$Event2.t2[i] <- n^(-1/2) * sum(trial.data$martingale.comp2.t2)
              
              results$Var.Event1.t1[i] <- sum(trial.data$Event1 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)))
              results$Var.Event2.t1[i] <- sum(trial.data$Event2 <= pmin(followup, pmax(0,t1 - trial.data$Arrival)))
              results$Var.Event1.t2[i] <- sum(trial.data$Event1 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)))
              results$Var.Event2.t2[i] <- sum(trial.data$Event2 <= pmin(followup, pmax(0,t2 - trial.data$Arrival)))
              
            }, error = function(e){cat("ERROR in simulation ", i, ": ", conditionMessage(e), sep = "", "\n")})
            
          }
          
          
          #Calculate final statistics
          
          results$Z1.t1 <- sqrt(n) * results$Event1.t1/sqrt(results$Var.Event1.t1)
          results$Z2.t1 <- sqrt(n) * results$Event2.t1/sqrt(results$Var.Event2.t1)
          results$Z1.t2 <- sqrt(n) * (results$Event1.t2 - results$Event1.t1)/sqrt(results$Var.Event1.t2 - results$Var.Event1.t1)
          results$Z2.t2 <- sqrt(n) * (results$Event2.t2 - results$Event2.t1)/sqrt(results$Var.Event2.t2 - results$Var.Event2.t1)
          
          
          emp.alpha.level[current.row.to.fill,"Sample size"] <- sample.size.temp
          emp.alpha.level[current.row.to.fill,"m1"] <- fixed.median
          emp.alpha.level[current.row.to.fill,"k1"] <- fixed.shape
          emp.alpha.level[current.row.to.fill,"m2"] <- med2.temp
          emp.alpha.level[current.row.to.fill,"k2"] <- kappa2.temp
          emp.alpha.level[current.row.to.fill,"Copula"] <- copula.temp
          emp.alpha.level[current.row.to.fill,"Copula parameter"] <- copula.parameters[intensity.temp, copula.temp]
          
          #compute empirical error level for standard and weighted maxnorm and L2 rejection regions and fill the corresponding cells
          
          results$rejection.maxnorm <- 1-(abs(results$Z1.t1)<u.maxnorm)*(abs(results$Z2.t1)<u.maxnorm)*(abs(results$Z1.t2)<u.maxnorm)*(abs(results$Z2.t2)<u.maxnorm)
          emp.alpha.level[current.row.to.fill,"alpha.maxnorm"] <- sum(results$rejection.maxnorm, na.rm = TRUE)/simulations
          
          results$rejection.L2 <- 1-((results$Z1.t1)^2 + (results$Z2.t1)^2 < u.L2)*((results$Z1.t2)^2 + (results$Z2.t2)^2 < u.L2)
          emp.alpha.level[current.row.to.fill,"alpha.L2"] <- sum(results$rejection.L2, na.rm = TRUE)/simulations
          
          results$rejection.maxnorm.weight1 <- 1-(abs(results$Z1.t1)<u1.maxnorm.weight1)*(abs(results$Z2.t1)<u1.maxnorm.weight1)*(abs(results$Z1.t2)<u2.maxnorm.weight1)*(abs(results$Z2.t2)<u2.maxnorm.weight1)
          emp.alpha.level[current.row.to.fill,"alpha.maxnorm.weight1"] <- sum(results$rejection.maxnorm.weight1, na.rm = TRUE)/simulations
          
          results$rejection.L2.weight1 <- 1-((results$Z1.t1)^2 + (results$Z2.t1)^2 < u1.L2.weight1)*((results$Z1.t2)^2 + (results$Z2.t2)^2 < u2.L2.weight1)
          emp.alpha.level[current.row.to.fill,"alpha.L2.weight1"] <- sum(results$rejection.L2.weight1, na.rm = TRUE)/simulations
          
          results$rejection.maxnorm.weight2 <- 1-(abs(results$Z1.t1)<u1.maxnorm.weight2)*(abs(results$Z2.t1)<u1.maxnorm.weight2)*(abs(results$Z1.t2)<u2.maxnorm.weight2)*(abs(results$Z2.t2)<u2.maxnorm.weight2)
          emp.alpha.level[current.row.to.fill,"alpha.maxnorm.weight2"] <- sum(results$rejection.maxnorm.weight2, na.rm = TRUE)/simulations
          
          results$rejection.L2.weight2 <- 1-((results$Z1.t1)^2 + (results$Z2.t1)^2 < u1.L2.weight2)*((results$Z1.t2)^2 + (results$Z2.t2)^2 < u2.L2.weight2)
          emp.alpha.level[current.row.to.fill,"alpha.L2.weight2"] <- sum(results$rejection.L2.weight2, na.rm = TRUE)/simulations
          
        }
      }
    }
  }
}