#============================================================#
# Simulation: Consistent Mediation with "Regular" Prop. Mediated
# R version 4.3.2 (2023-10-31 ucrt)
# RStudio 2023.09.1+494 "Desert Sunflower" Release (cd7011dce393115d3a7c3db799dda4b1c7e88711, 2023-10-16) for windows
# Last updated = 2023-11-14
#============================================================#
### Load packages
library(lavaan) #v0.6-16
library(semTools) #v0.5-6
library(tidyverse) #v2.0.0

# For full code runtime calculation
fullrun.start <- Sys.time()

#------------------------------------------------------------#
### Specify Population and Estimation Models
#------------------------------------------------------------#
# Population Models for data generation
pop.m1 <- ' m ~ 0.138564*x
            y ~ 0.138564*m + 0.1408*x  '
pop.m2 <- ' m ~ 0.195959*x
            y ~ 0.195959*m + 0.1216*x  '
pop.m3 <- ' m ~ 0.246577*x
            y ~ 0.246577*m + 0.0992*x  '
pop.m4 <- ' m ~ 0.3020*x
            y ~ 0.3020*m + 0.0688*x  '
pop.m5 <- ' m ~ 0.348712*x
            y ~ 0.348712*m + 0.0384*x  '

# Estimation Model -- syntax guide @ https://lavaan.ugent.be/tutorial/
estim <- '  m ~ a*x
            y ~ b*m + c*x 
            PM := (a*b)/(c+a*b)'     # direct effect labeled 'c' for convenience. For absPM use: 'AbsPM := abs(a*b)/(abs(c) + abs(a*b))'

# Assign population PM Values
pop1.PM <- 0.12
pop2.PM <- 0.24
pop3.PM <- 0.38
pop4.PM <- 0.57
pop5.PM <- 0.76

# Specify number of replications for the simulation function and the Monte Carlo method.
# For each Model x Sample Size combination, we run the simulation 'n' x 'rep' times.
Nrep = 1000 # number of simulation replications. Typically 1000
MCrep = 10000 # number of Monte Carlo replications used to find the 95% CI

# List all models and sample sizes combinations
model.list = c("m1n","m2n","m3n","m4n","m5n")
samplesize.list = c(100,300,500,1000)
sim.list = c(outer(model.list, samplesize.list, paste, sep = "")) # list of model x sample sizes combinations
stats.list = c(outer(sim.list, ".stats", paste, sep = "")) # list of output stats to combine.

# Set seed globally
set.seed(555)

#--------------------------------------------------------------------------------#
### Set up simulation functions and run for 'Nrep' times.
# It'd be more efficient to use loop/map/apply() to iterate through various 
# sample sizes and population models, but I didn't have time to figure this out.
#--------------------------------------------------------------------------------#
#=== Population Model 1 ===#
# N = 100
m1n100.PM.func = function(  ){
  m1n100.dta <- simulateData(pop.m1, sample.nobs=100) # generate data
  fit.m1n100 <- sem(estim, data=m1n100.dta) # analyze
  m1n100.PM <- monteCarloCI(fit.m1n100, rep = MCrep) }
# N = 300
m1n300.PM.func = function(  ){
  m1n300.dta <- simulateData(pop.m1, sample.nobs=300)
  fit.m1n300 <- sem(estim, data=m1n300.dta)
  m1n300.PM <- monteCarloCI(fit.m1n300, rep = MCrep) }
# N = 500
m1n500.PM.func = function(  ){
  m1n500.dta <- simulateData(pop.m1, sample.nobs=500)
  fit.m1n500 <- sem(estim, data=m1n500.dta)
  m1n500.PM <- monteCarloCI(fit.m1n500, rep = MCrep) }
# N = 1000
m1n1000.PM.func = function(  ){
  m1n1000.dta <- simulateData(pop.m1, sample.nobs=1000)
  fit.m1n1000 <- sem(estim, data=m1n1000.dta)
  m1n1000.PM <- monteCarloCI(fit.m1n1000, rep = MCrep) }
# Perform the simulations
# If simplify = TRUE, then you get a matrix of values instead of a list. Each column in the matrix is the output from one run of the function, so it's harder to extract individual values.
m1n100reps.PM = replicate(n = Nrep, m1n100.PM.func(), simplify = F)
m1n300reps.PM = replicate(n = Nrep, m1n300.PM.func(), simplify = F)
m1n500reps.PM = replicate(n = Nrep, m1n500.PM.func(), simplify = F)
m1n1000reps.PM = replicate(n = Nrep, m1n1000.PM.func(), simplify = F)

#=== Population Model 2 ===#
# N = 100
m2n100.PM.func = function(  ){
  m2n100.dta <- simulateData(pop.m2, sample.nobs=100)
  fit.m2n100 <- sem(estim, data=m2n100.dta)
  m2n100.PM <- monteCarloCI(fit.m2n100, rep = MCrep) }
# N = 300
m2n300.PM.func = function(  ){
  m2n300.dta <- simulateData(pop.m2, sample.nobs=300)
  fit.m2n300 <- sem(estim, data=m2n300.dta)
  m2n300.PM <- monteCarloCI(fit.m2n300, rep = MCrep) }
# N = 500
m2n500.PM.func = function(  ){
  m2n500.dta <- simulateData(pop.m2, sample.nobs=500)
  fit.m2n500 <- sem(estim, data=m2n500.dta)
  m2n500.PM <- monteCarloCI(fit.m2n500, rep = MCrep) }
# N = 1000
m2n1000.PM.func = function(  ){
  m2n1000.dta <- simulateData(pop.m2, sample.nobs=1000)
  fit.m2n1000 <- sem(estim, data=m2n1000.dta)
  m2n1000.PM <- monteCarloCI(fit.m2n1000, rep = MCrep) }
# Perform the simulations
m2n100reps.PM = replicate(n = Nrep, m2n100.PM.func(), simplify = F)
m2n300reps.PM = replicate(n = Nrep, m2n300.PM.func(), simplify = F)
m2n500reps.PM = replicate(n = Nrep, m2n500.PM.func(), simplify = F)
m2n1000reps.PM = replicate(n = Nrep, m2n1000.PM.func(), simplify = F)
    
#=== Population Model 3 ===#
# N = 100
m3n100.PM.func = function(  ){
  m3n100.dta <- simulateData(pop.m3, sample.nobs=100)
  fit.m3n100 <- sem(estim, data=m3n100.dta)
  m3n100.PM <- monteCarloCI(fit.m3n100, rep = MCrep) }
# N = 300
m3n300.PM.func = function(  ){
  m3n300.dta <- simulateData(pop.m3, sample.nobs=300)
  fit.m3n300 <- sem(estim, data=m3n300.dta)
  m3n300.PM <- monteCarloCI(fit.m3n300, rep = MCrep) }
# N = 500
m3n500.PM.func = function(  ){
  m3n500.dta <- simulateData(pop.m3, sample.nobs=500)
  fit.m3n500 <- sem(estim, data=m3n500.dta)
  m3n500.PM <- monteCarloCI(fit.m3n500, rep = MCrep) }
# N = 1000
m3n1000.PM.func = function(  ){
  m3n1000.dta <- simulateData(pop.m3, sample.nobs=1000)
  fit.m3n1000 <- sem(estim, data=m3n1000.dta)
  m3n1000.PM <- monteCarloCI(fit.m3n1000, rep = MCrep) }
# Perform the simulations
m3n100reps.PM = replicate(n = Nrep, m3n100.PM.func(), simplify = F)
m3n300reps.PM = replicate(n = Nrep, m3n300.PM.func(), simplify = F)
m3n500reps.PM = replicate(n = Nrep, m3n500.PM.func(), simplify = F)
m3n1000reps.PM = replicate(n = Nrep, m3n1000.PM.func(), simplify = F)

#=== Population Model 4 ===#
# N = 100
m4n100.PM.func = function(  ){
  m4n100.dta <- simulateData(pop.m4, sample.nobs=100)
  fit.m4n100 <- sem(estim, data=m4n100.dta)
  m4n100.PM <- monteCarloCI(fit.m4n100, rep = MCrep) }
# N = 300
m4n300.PM.func = function(  ){
  m4n300.dta <- simulateData(pop.m4, sample.nobs=300)
  fit.m4n300 <- sem(estim, data=m4n300.dta)
  m4n300.PM <- monteCarloCI(fit.m4n300, rep = MCrep) }
# N = 500
m4n500.PM.func = function(  ){
  m4n500.dta <- simulateData(pop.m4, sample.nobs=500)
  fit.m4n500 <- sem(estim, data=m4n500.dta)
  m4n500.PM <- monteCarloCI(fit.m4n500, rep = MCrep) }
# N = 1000
m4n1000.PM.func = function(  ){
  m4n1000.dta <- simulateData(pop.m4, sample.nobs=1000)
  fit.m4n1000 <- sem(estim, data=m4n1000.dta)
  m4n1000.PM <- monteCarloCI(fit.m4n1000, rep = MCrep) }
# Perform the simulations
m4n100reps.PM = replicate(n = Nrep, m4n100.PM.func(), simplify = F)
m4n300reps.PM = replicate(n = Nrep, m4n300.PM.func(), simplify = F)
m4n500reps.PM = replicate(n = Nrep, m4n500.PM.func(), simplify = F)
m4n1000reps.PM = replicate(n = Nrep, m4n1000.PM.func(), simplify = F)

#=== Population Model 5 ===#
# N = 100
m5n100.PM.func = function(  ){
  m5n100.dta <- simulateData(pop.m5, sample.nobs=100)
  fit.m5n100 <- sem(estim, data=m5n100.dta)
  m5n100.PM <- monteCarloCI(fit.m5n100, rep = MCrep) }
# N = 300
m5n300.PM.func = function(  ){
  m5n300.dta <- simulateData(pop.m5, sample.nobs=300)
  fit.m5n300 <- sem(estim, data=m5n300.dta)
  m5n300.PM <- monteCarloCI(fit.m5n300, rep = MCrep) }
# N = 500
m5n500.PM.func = function(  ){
  m5n500.dta <- simulateData(pop.m5, sample.nobs=500)
  fit.m5n500 <- sem(estim, data=m5n500.dta)
  m5n500.PM <- monteCarloCI(fit.m5n500, rep = MCrep) }
# N = 1000
m5n1000.PM.func = function(  ){
  m5n1000.dta <- simulateData(pop.m5, sample.nobs=1000)
  fit.m5n1000 <- sem(estim, data=m5n1000.dta)
  m5n1000.PM <- monteCarloCI(fit.m5n1000, rep = MCrep) }
# Perform the simulations
m5n100reps.PM = replicate(n = Nrep, m5n100.PM.func(), simplify = F)
m5n300reps.PM = replicate(n = Nrep, m5n300.PM.func(), simplify = F)
m5n500reps.PM = replicate(n = Nrep, m5n500.PM.func(), simplify = F)
m5n1000reps.PM = replicate(n = Nrep, m5n1000.PM.func(), simplify = F)

#-------------------------------------------------------------#
### Model 1 Extract Sim Info - Basic PM
#-------------------------------------------------------------#
#=== N = 100 ===#
m1n100.PMest = as.data.frame(t(sapply(m1n100reps.PM, unlist))) # convert list to dataframe. t() transposes the df from wide to long.
m1n100.PMest = m1n100.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM via dplyr filter()
m1n100.PM.bias = mean(m1n100.PMest$est) - pop1.PM # Avg. bias
m1n100.PM.relbias = (mean(m1n100.PMest$est) - pop1.PM)/(pop1.PM) # Relative bias
m1n100.PM.sd = sd(m1n100.PMest$est) # Efficiency (SD)
m1n100.PMest$ci.width = m1n100.PMest$ci.upper - m1n100.PMest$ci.lower # Calculate CI width
m1n100.power = sum(m1n100.PMest$ci.lower>0) / nrow(m1n100.PMest) # Power = proportion of reps in which the CI do not include zero when the true value is nonzero
# We divide by nrow since we don't know beforehand how many rows will have PMs out of bound.
m1n100.cover = sum(cbind(m1n100.PMest$ci.lower <= pop1.PM & m1n100.PMest$ci.upper >= pop1.PM), na.rm = TRUE) / nrow(m1n100.PMest) # Coverage = proportion of reps in which CIs include the true value.
### Store all info in a new data.frame
m1n100.stats = data.frame(m1n100.PM.bias, m1n100.PM.relbias, m1n100.PM.sd, m1n100.power, m1n100.cover, mean(m1n100.PMest$ci.width), median(m1n100.PMest$ci.width), mean(m1n100.PMest$ci.width)/pop1.PM)
rownames(m1n100.stats) <- c("m1n100")
#=== N = 300 ===#
m1n300.PMest = as.data.frame(t(sapply(m1n300reps.PM, unlist))) # Convert list to dataframe
m1n300.PMest = m1n300.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m1n300.PM.bias = mean(m1n300.PMest$est) - pop1.PM # Avg. bias
m1n300.PM.relbias = (mean(m1n300.PMest$est) - pop1.PM)/(pop1.PM) # Relative bias
m1n300.PM.sd = sd(m1n300.PMest$est) # Efficiency (SD)
m1n300.PMest$ci.width = m1n300.PMest$ci.upper - m1n300.PMest$ci.lower # CI width
m1n300.power = sum(m1n300.PMest$ci.lower>0) / nrow(m1n300.PMest) # Power
m1n300.cover = sum(cbind(m1n300.PMest$ci.lower <= pop1.PM & m1n300.PMest$ci.upper >= pop1.PM), na.rm = TRUE) / nrow(m1n300.PMest) # Coverage
### Store all info in a new data.frame
m1n300.stats = data.frame(m1n300.PM.bias, m1n300.PM.relbias, m1n300.PM.sd, m1n300.power, m1n300.cover, mean(m1n300.PMest$ci.width), median(m1n300.PMest$ci.width), mean(m1n300.PMest$ci.width)/pop1.PM)
rownames(m1n300.stats) <- c("m1n300")
#=== N = 500 ===#
m1n500.PMest = as.data.frame(t(sapply(m1n500reps.PM, unlist))) # Convert list to dataframe
m1n500.PMest = m1n500.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m1n500.PM.bias = mean(m1n500.PMest$est) - pop1.PM # Avg. bias
m1n500.PM.relbias = (mean(m1n500.PMest$est) - pop1.PM)/(pop1.PM) # Relative bias
m1n500.PM.sd = sd(m1n500.PMest$est) # Efficiency (SD)
m1n500.PMest$ci.width = m1n500.PMest$ci.upper - m1n500.PMest$ci.lower # CI width
m1n500.power = sum(m1n500.PMest$ci.lower>0) / nrow(m1n500.PMest) # Power
m1n500.cover = sum(cbind(m1n500.PMest$ci.lower <= pop1.PM & m1n500.PMest$ci.upper >= pop1.PM), na.rm = TRUE) / nrow(m1n500.PMest) # Coverage
### Store all info in a new data.frame
m1n500.stats = data.frame(m1n500.PM.bias, m1n500.PM.relbias, m1n500.PM.sd, m1n500.power, m1n500.cover, mean(m1n500.PMest$ci.width), median(m1n500.PMest$ci.width), mean(m1n500.PMest$ci.width)/pop1.PM)
rownames(m1n500.stats) <- c("m1n500")
#=== N = 1000 ===#
m1n1000.PMest = as.data.frame(t(sapply(m1n1000reps.PM, unlist))) # Convert list to dataframe
m1n1000.PMest = m1n1000.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m1n1000.PM.bias = mean(m1n1000.PMest$est) - pop1.PM # Avg. bias
m1n1000.PM.relbias = (mean(m1n1000.PMest$est) - pop1.PM)/(pop1.PM) # Relative bias
m1n1000.PM.sd = sd(m1n1000.PMest$est) # Efficiency (SD)
m1n1000.PMest$ci.width = m1n1000.PMest$ci.upper - m1n1000.PMest$ci.lower # CI width
m1n1000.power = sum(m1n1000.PMest$ci.lower>0) / nrow(m1n1000.PMest) # Power
m1n1000.cover = sum(cbind(m1n1000.PMest$ci.lower <= pop1.PM & m1n1000.PMest$ci.upper >= pop1.PM), na.rm = TRUE) / nrow(m1n1000.PMest) # Coverage
### Store all info in a new data.frame
m1n1000.stats = data.frame(m1n1000.PM.bias, m1n1000.PM.relbias, m1n1000.PM.sd, m1n1000.power, m1n1000.cover, mean(m1n1000.PMest$ci.width), median(m1n1000.PMest$ci.width), mean(m1n1000.PMest$ci.width)/pop1.PM)
rownames(m1n1000.stats) <- c("m1n1000")

#-------------------------------------------------------------#
### Model 2 Extract Sim Info - Basic PM
#-------------------------------------------------------------#
#=== N = 100 ===#
m2n100.PMest = as.data.frame(t(sapply(m2n100reps.PM, unlist))) # Convert list to dataframe
m2n100.PMest = m2n100.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m2n100.PM.bias = mean(m2n100.PMest$est) - pop2.PM # Avg. bias
m2n100.PM.relbias = (mean(m2n100.PMest$est) - pop2.PM)/(pop2.PM) # Relative bias
m2n100.PM.sd = sd(m2n100.PMest$est) # Efficiency (SD)
m2n100.PMest$ci.width = m2n100.PMest$ci.upper - m2n100.PMest$ci.lower # Calculate CI width
m2n100.power = sum(m2n100.PMest$ci.lower>0) / nrow(m2n100.PMest) # Power = proportion of reps in which the CI do not include zero when the true value is nonzero
# We divide by nrow since we don't know beforehand how many rows will have PMs out of bound.
m2n100.cover = sum(cbind(m2n100.PMest$ci.lower <= pop2.PM & m2n100.PMest$ci.upper >= pop2.PM), na.rm = TRUE) / nrow(m2n100.PMest) # Coverage = proportion of reps in which CIs include the true value.
### Store all info in a new data.frame
m2n100.stats = data.frame(m2n100.PM.bias, m2n100.PM.relbias, m2n100.PM.sd, m2n100.power, m2n100.cover, mean(m2n100.PMest$ci.width), median(m2n100.PMest$ci.width), mean(m2n100.PMest$ci.width)/pop2.PM)
rownames(m2n100.stats) <- c("m2n100")
#=== N = 300 ===#
m2n300.PMest = as.data.frame(t(sapply(m2n300reps.PM, unlist))) # Convert list to dataframe
m2n300.PMest = m2n300.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m2n300.PM.bias = mean(m2n300.PMest$est) - pop2.PM # Avg. bias
m2n300.PM.relbias = (mean(m2n300.PMest$est) - pop2.PM)/(pop2.PM) # Relative bias
m2n300.PM.sd = sd(m2n300.PMest$est) # Efficiency (SD)
m2n300.PMest$ci.width = m2n300.PMest$ci.upper - m2n300.PMest$ci.lower # CI width
m2n300.power = sum(m2n300.PMest$ci.lower>0) / nrow(m2n300.PMest) # Power
m2n300.cover = sum(cbind(m2n300.PMest$ci.lower <= pop2.PM & m2n300.PMest$ci.upper >= pop2.PM), na.rm = TRUE) / nrow(m2n300.PMest) # Coverage
### Store all info in a new data.frame
m2n300.stats = data.frame(m2n300.PM.bias, m2n300.PM.relbias, m2n300.PM.sd, m2n300.power, m2n300.cover, mean(m2n300.PMest$ci.width), median(m2n300.PMest$ci.width), mean(m2n300.PMest$ci.width)/pop2.PM)
rownames(m2n300.stats) <- c("m2n300")
#=== N = 500 ===#
m2n500.PMest = as.data.frame(t(sapply(m2n500reps.PM, unlist))) # Convert list to dataframe
m2n500.PMest = m2n500.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m2n500.PM.bias = mean(m2n500.PMest$est) - pop2.PM # Avg. bias
m2n500.PM.relbias = (mean(m2n500.PMest$est) - pop2.PM)/(pop2.PM) # Relative bias
m2n500.PM.sd = sd(m2n500.PMest$est) # Efficiency (SD)
m2n500.PMest$ci.width = m2n500.PMest$ci.upper - m2n500.PMest$ci.lower # CI width
m2n500.power = sum(m2n500.PMest$ci.lower>0) / nrow(m2n500.PMest) # Power
m2n500.cover = sum(cbind(m2n500.PMest$ci.lower <= pop2.PM & m2n500.PMest$ci.upper >= pop2.PM), na.rm = TRUE) / nrow(m2n500.PMest) # Coverage
### Store all info in a new data.frame
m2n500.stats = data.frame(m2n500.PM.bias, m2n500.PM.relbias, m2n500.PM.sd, m2n500.power, m2n500.cover, mean(m2n500.PMest$ci.width), median(m2n500.PMest$ci.width), mean(m2n500.PMest$ci.width)/pop2.PM)
rownames(m2n500.stats) <- c("m2n500")
#=== N = 1000 ===#
m2n1000.PMest = as.data.frame(t(sapply(m2n1000reps.PM, unlist))) # Convert list to dataframe
m2n1000.PMest = m2n1000.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m2n1000.PM.bias = mean(m2n1000.PMest$est) - pop2.PM # Avg. bias
m2n1000.PM.relbias = (mean(m2n1000.PMest$est) - pop2.PM)/(pop2.PM) # Relative bias
m2n1000.PM.sd = sd(m2n1000.PMest$est) # Efficiency (SD)
m2n1000.PMest$ci.width = m2n1000.PMest$ci.upper - m2n1000.PMest$ci.lower # CI width
m2n1000.power = sum(m2n1000.PMest$ci.lower>0) / nrow(m2n1000.PMest) # Power
m2n1000.cover = sum(cbind(m2n1000.PMest$ci.lower <= pop2.PM & m2n1000.PMest$ci.upper >= pop2.PM), na.rm = TRUE) / nrow(m2n1000.PMest) # Coverage
### Store all info in a new data.frame
m2n1000.stats = data.frame(m2n1000.PM.bias, m2n1000.PM.relbias, m2n1000.PM.sd, m2n1000.power, m2n1000.cover, mean(m2n1000.PMest$ci.width), median(m2n1000.PMest$ci.width), mean(m2n1000.PMest$ci.width)/pop2.PM)
rownames(m2n1000.stats) <- c("m2n1000")

#-------------------------------------------------------------#
### Model 3 Extract Sim Info - Basic PM
#-------------------------------------------------------------#
#=== N = 100 ===#
m3n100.PMest = as.data.frame(t(sapply(m3n100reps.PM, unlist))) # Convert list to dataframe
m3n100.PMest = m3n100.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m3n100.PM.bias = mean(m3n100.PMest$est) - pop3.PM # Avg. bias
m3n100.PM.relbias = (mean(m3n100.PMest$est) - pop3.PM)/(pop3.PM) # Relative bias
m3n100.PM.sd = sd(m3n100.PMest$est) # Efficiency (SD)
m3n100.PMest$ci.width = m3n100.PMest$ci.upper - m3n100.PMest$ci.lower # Calculate CI width
m3n100.power = sum(m3n100.PMest$ci.lower>0) / nrow(m3n100.PMest) # Power
m3n100.cover = sum(cbind(m3n100.PMest$ci.lower <= pop3.PM & m3n100.PMest$ci.upper >= pop3.PM), na.rm = TRUE) / nrow(m3n100.PMest) # Coverage
### Store all info in a new data.frame
m3n100.stats = data.frame(m3n100.PM.bias, m3n100.PM.relbias, m3n100.PM.sd, m3n100.power, m3n100.cover, mean(m3n100.PMest$ci.width), median(m3n100.PMest$ci.width), mean(m3n100.PMest$ci.width)/pop3.PM)
rownames(m3n100.stats) <- c("m3n100")
#=== N = 300 ===#
m3n300.PMest = as.data.frame(t(sapply(m3n300reps.PM, unlist))) # Convert list to dataframe
m3n300.PMest = m3n300.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m3n300.PM.bias = mean(m3n300.PMest$est) - pop3.PM # Avg. bias
m3n300.PM.relbias = (mean(m3n300.PMest$est) - pop3.PM)/(pop3.PM) # Relative bias
m3n300.PM.sd = sd(m3n300.PMest$est) # Efficiency (SD)
m3n300.PMest$ci.width = m3n300.PMest$ci.upper - m3n300.PMest$ci.lower # CI width
m3n300.power = sum(m3n300.PMest$ci.lower>0) / nrow(m3n300.PMest) # Power
m3n300.cover = sum(cbind(m3n300.PMest$ci.lower <= pop3.PM & m3n300.PMest$ci.upper >= pop3.PM), na.rm = TRUE) / nrow(m3n300.PMest) # Coverage
### Store all info in a new data.frame
m3n300.stats = data.frame(m3n300.PM.bias, m3n300.PM.relbias, m3n300.PM.sd, m3n300.power, m3n300.cover, mean(m3n300.PMest$ci.width), median(m3n300.PMest$ci.width), mean(m3n300.PMest$ci.width)/pop3.PM)
rownames(m3n300.stats) <- c("m3n300")
#=== N = 500 ===#
m3n500.PMest = as.data.frame(t(sapply(m3n500reps.PM, unlist))) # # Convert list to dataframe
m3n500.PMest = m3n500.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m3n500.PM.bias = mean(m3n500.PMest$est) - pop3.PM # Avg. bias
m3n500.PM.relbias = (mean(m3n500.PMest$est) - pop3.PM)/(pop3.PM) # Relative bias
m3n500.PM.sd = sd(m3n500.PMest$est) # Efficiency (SD)
m3n500.PMest$ci.width = m3n500.PMest$ci.upper - m3n500.PMest$ci.lower # CI width
m3n500.power = sum(m3n500.PMest$ci.lower>0) / nrow(m3n500.PMest) # Power
m3n500.cover = sum(cbind(m3n500.PMest$ci.lower <= pop3.PM & m3n500.PMest$ci.upper >= pop3.PM), na.rm = TRUE) / nrow(m3n500.PMest) # Coverage
### Store all info in a new data.frame
m3n500.stats = data.frame(m3n500.PM.bias, m3n500.PM.relbias, m3n500.PM.sd, m3n500.power, m3n500.cover, mean(m3n500.PMest$ci.width), median(m3n500.PMest$ci.width), mean(m3n500.PMest$ci.width)/pop3.PM)
rownames(m3n500.stats) <- c("m3n500")
#=== N = 1000 ===#
m3n1000.PMest = as.data.frame(t(sapply(m3n1000reps.PM, unlist))) # Convert list to dataframe
m3n1000.PMest = m3n1000.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m3n1000.PM.bias = mean(m3n1000.PMest$est) - pop3.PM # Avg. bias
m3n1000.PM.relbias = (mean(m3n1000.PMest$est) - pop3.PM)/(pop3.PM) # Relative bias
m3n1000.PM.sd = sd(m3n1000.PMest$est) # Efficiency (SD)
m3n1000.PMest$ci.width = m3n1000.PMest$ci.upper - m3n1000.PMest$ci.lower # CI width
m3n1000.power = sum(m3n1000.PMest$ci.lower>0) / nrow(m3n1000.PMest) # Power
m3n1000.cover = sum(cbind(m3n1000.PMest$ci.lower <= pop3.PM & m3n1000.PMest$ci.upper >= pop3.PM), na.rm = TRUE) / nrow(m3n1000.PMest) # Coverage
### Store all info in a new data.frame
m3n1000.stats = data.frame(m3n1000.PM.bias, m3n1000.PM.relbias, m3n1000.PM.sd, m3n1000.power, m3n1000.cover, mean(m3n1000.PMest$ci.width), median(m3n1000.PMest$ci.width), mean(m3n1000.PMest$ci.width)/pop3.PM)
rownames(m3n1000.stats) <- c("m3n1000")

#-------------------------------------------------------------#
### Model 4 Extract Sim Info - Basic PM
#-------------------------------------------------------------#
#=== N = 100 ===#
m4n100.PMest = as.data.frame(t(sapply(m4n100reps.PM, unlist))) # Convert list to dataframe
m4n100.PMest = m4n100.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m4n100.PM.bias = mean(m4n100.PMest$est) - pop4.PM # Avg. bias
m4n100.PM.relbias = (mean(m4n100.PMest$est) - pop4.PM)/(pop4.PM) # Relative bias
m4n100.PM.sd = sd(m4n100.PMest$est) # Efficiency (SD)
m4n100.PMest$ci.width = m4n100.PMest$ci.upper - m4n100.PMest$ci.lower # Calculate CI width
m4n100.power = sum(m4n100.PMest$ci.lower>0) / nrow(m4n100.PMest) # Power
m4n100.cover = sum(cbind(m4n100.PMest$ci.lower <= pop4.PM & m4n100.PMest$ci.upper >= pop4.PM), na.rm = TRUE) / nrow(m4n100.PMest) # Coverage = proportion of reps in which CIs include the true value.
### Store all info in a new data.frame
m4n100.stats = data.frame(m4n100.PM.bias, m4n100.PM.relbias, m4n100.PM.sd, m4n100.power, m4n100.cover, mean(m4n100.PMest$ci.width), median(m4n100.PMest$ci.width), mean(m4n100.PMest$ci.width)/pop4.PM)
rownames(m4n100.stats) <- c("m4n100")
#=== N = 300 ===#
m4n300.PMest = as.data.frame(t(sapply(m4n300reps.PM, unlist))) # Convert list to dataframe
m4n300.PMest = m4n300.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m4n300.PM.bias = mean(m4n300.PMest$est) - pop4.PM # Avg. bias
m4n300.PM.relbias = (mean(m4n300.PMest$est) - pop4.PM)/(pop4.PM) # Relative bias
m4n300.PM.sd = sd(m4n300.PMest$est) # Efficiency (SD)
m4n300.PMest$ci.width = m4n300.PMest$ci.upper - m4n300.PMest$ci.lower # CI width
m4n300.power = sum(m4n300.PMest$ci.lower>0) / nrow(m4n300.PMest) # Power
m4n300.cover = sum(cbind(m4n300.PMest$ci.lower <= pop4.PM & m4n300.PMest$ci.upper >= pop4.PM), na.rm = TRUE) / nrow(m4n300.PMest) # Coverage
### Store all info in a new data.frame
m4n300.stats = data.frame(m4n300.PM.bias, m4n300.PM.relbias, m4n300.PM.sd, m4n300.power, m4n300.cover, mean(m4n300.PMest$ci.width), median(m4n300.PMest$ci.width), mean(m4n300.PMest$ci.width)/pop4.PM)
rownames(m4n300.stats) <- c("m4n300")
#=== N = 500 ===#
m4n500.PMest = as.data.frame(t(sapply(m4n500reps.PM, unlist))) # # Convert list to dataframe
m4n500.PMest = m4n500.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m4n500.PM.bias = mean(m4n500.PMest$est) - pop4.PM # Avg. bias
m4n500.PM.relbias = (mean(m4n500.PMest$est) - pop4.PM)/(pop4.PM) # Relative bias
m4n500.PM.sd = sd(m4n500.PMest$est) # Efficiency (SD)
m4n500.PMest$ci.width = m4n500.PMest$ci.upper - m4n500.PMest$ci.lower # CI width
m4n500.power = sum(m4n500.PMest$ci.lower>0) / nrow(m4n500.PMest) # Power
m4n500.cover = sum(cbind(m4n500.PMest$ci.lower <= pop4.PM & m4n500.PMest$ci.upper >= pop4.PM), na.rm = TRUE) / nrow(m4n500.PMest) # Coverage
### Store all info in a new data.frame
m4n500.stats = data.frame(m4n500.PM.bias, m4n500.PM.relbias, m4n500.PM.sd, m4n500.power, m4n500.cover, mean(m4n500.PMest$ci.width), median(m4n500.PMest$ci.width), mean(m4n500.PMest$ci.width)/pop4.PM)
rownames(m4n500.stats) <- c("m4n500")
#=== N = 1000 ===#
m4n1000.PMest = as.data.frame(t(sapply(m4n1000reps.PM, unlist))) # Convert list to dataframe
m4n1000.PMest = m4n1000.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m4n1000.PM.bias = mean(m4n1000.PMest$est) - pop4.PM # Avg. bias
m4n1000.PM.relbias = (mean(m4n1000.PMest$est) - pop4.PM)/(pop4.PM) # Relative bias
m4n1000.PM.sd = sd(m4n1000.PMest$est) # Efficiency (SD)
m4n1000.PMest$ci.width = m4n1000.PMest$ci.upper - m4n1000.PMest$ci.lower # CI width
m4n1000.power = sum(m4n1000.PMest$ci.lower>0) / nrow(m4n1000.PMest) # Power
m4n1000.cover = sum(cbind(m4n1000.PMest$ci.lower <= pop4.PM & m4n1000.PMest$ci.upper >= pop4.PM), na.rm = TRUE) / nrow(m4n1000.PMest) # Coverage
### Store all info in a new data.frame
m4n1000.stats = data.frame(m4n1000.PM.bias, m4n1000.PM.relbias, m4n1000.PM.sd, m4n1000.power, m4n1000.cover, mean(m4n1000.PMest$ci.width), median(m4n1000.PMest$ci.width), mean(m4n1000.PMest$ci.width)/pop4.PM)
rownames(m4n1000.stats) <- c("m4n1000")

#-------------------------------------------------------------#
### Model 5 Extract Sim Info - Basic PM
#-------------------------------------------------------------#
#=== N = 100 ===#
m5n100.PMest = as.data.frame(t(sapply(m5n100reps.PM, unlist))) # Convert list to dataframe
m5n100.PMest = m5n100.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m5n100.PM.bias = mean(m5n100.PMest$est) - pop5.PM # Avg. bias
m5n100.PM.relbias = (mean(m5n100.PMest$est) - pop5.PM)/(pop5.PM) # Relative bias
m5n100.PM.sd = sd(m5n100.PMest$est) # Efficiency (SD)
m5n100.PMest$ci.width = m5n100.PMest$ci.upper - m5n100.PMest$ci.lower # Calculate CI width
m5n100.power = sum(m5n100.PMest$ci.lower>0) / nrow(m5n100.PMest) # Power
m5n100.cover = sum(cbind(m5n100.PMest$ci.lower <= pop5.PM & m5n100.PMest$ci.upper >= pop5.PM), na.rm = TRUE) / nrow(m5n100.PMest) # Coverage = proportion of reps in which CIs include the true value.
### Store all info in a new data.frame
m5n100.stats = data.frame(m5n100.PM.bias, m5n100.PM.relbias, m5n100.PM.sd, m5n100.power, m5n100.cover, mean(m5n100.PMest$ci.width), median(m5n100.PMest$ci.width), mean(m5n100.PMest$ci.width)/pop5.PM)
rownames(m5n100.stats) <- c("m5n100")
#=== N = 300 ===#
m5n300.PMest = as.data.frame(t(sapply(m5n300reps.PM, unlist))) # Convert list to dataframe
m5n300.PMest = m5n300.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m5n300.PM.bias = mean(m5n300.PMest$est) - pop5.PM # Avg. bias
m5n300.PM.relbias = (mean(m5n300.PMest$est) - pop5.PM)/(pop5.PM) # Relative bias
m5n300.PM.sd = sd(m5n300.PMest$est) # Efficiency (SD)
m5n300.PMest$ci.width = m5n300.PMest$ci.upper - m5n300.PMest$ci.lower # CI width
m5n300.power = sum(m5n300.PMest$ci.lower>0) / nrow(m5n300.PMest) # Power
m5n300.cover = sum(cbind(m5n300.PMest$ci.lower <= pop5.PM & m5n300.PMest$ci.upper >= pop5.PM), na.rm = TRUE) / nrow(m5n300.PMest) # Coverage
### Store all info in a new data.frame
m5n300.stats = data.frame(m5n300.PM.bias, m5n300.PM.relbias, m5n300.PM.sd, m5n300.power, m5n300.cover, mean(m5n300.PMest$ci.width), median(m5n300.PMest$ci.width), mean(m5n300.PMest$ci.width)/pop5.PM)
rownames(m5n300.stats) <- c("m5n300")
#=== N = 500 ===#
m5n500.PMest = as.data.frame(t(sapply(m5n500reps.PM, unlist))) # # Convert list to dataframe
m5n500.PMest = m5n500.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m5n500.PM.bias = mean(m5n500.PMest$est) - pop5.PM # Avg. bias
m5n500.PM.relbias = (mean(m5n500.PMest$est) - pop5.PM)/(pop5.PM) # Relative bias
m5n500.PM.sd = sd(m5n500.PMest$est) # Efficiency (SD)
m5n500.PMest$ci.width = m5n500.PMest$ci.upper - m5n500.PMest$ci.lower # CI width
m5n500.power = sum(m5n500.PMest$ci.lower>0) / nrow(m5n500.PMest) # Power
m5n500.cover = sum(cbind(m5n500.PMest$ci.lower <= pop5.PM & m5n500.PMest$ci.upper >= pop5.PM), na.rm = TRUE) / nrow(m5n500.PMest) # Coverage
### Store all info in a new data.frame
m5n500.stats = data.frame(m5n500.PM.bias, m5n500.PM.relbias, m5n500.PM.sd, m5n500.power, m5n500.cover, mean(m5n500.PMest$ci.width), median(m5n500.PMest$ci.width), mean(m5n500.PMest$ci.width)/pop5.PM)
rownames(m5n500.stats) <- c("m5n500")
#=== N = 1000 ===#
m5n1000.PMest = as.data.frame(t(sapply(m5n1000reps.PM, unlist))) # Convert list to dataframe
m5n1000.PMest = m5n1000.PMest %>% filter(est>=0 & est <=1) # Remove out of bound PM
m5n1000.PM.bias = mean(m5n1000.PMest$est) - pop5.PM # Avg. bias
m5n1000.PM.relbias = (mean(m5n1000.PMest$est) - pop5.PM)/(pop5.PM) # Relative bias
m5n1000.PM.sd = sd(m5n1000.PMest$est) # Efficiency (SD)
m5n1000.PMest$ci.width = m5n1000.PMest$ci.upper - m5n1000.PMest$ci.lower # CI width
m5n1000.power = sum(m5n1000.PMest$ci.lower>0) / nrow(m5n1000.PMest) # Power
m5n1000.cover = sum(cbind(m5n1000.PMest$ci.lower <= pop5.PM & m5n1000.PMest$ci.upper >= pop5.PM), na.rm = TRUE) / nrow(m5n1000.PMest) # Coverage
### Store all info in a new data.frame
m5n1000.stats = data.frame(m5n1000.PM.bias, m5n1000.PM.relbias, m5n1000.PM.sd, m5n1000.power, m5n1000.cover, mean(m5n1000.PMest$ci.width), median(m5n1000.PMest$ci.width), mean(m5n1000.PMest$ci.width)/pop5.PM)
rownames(m5n1000.stats) <- c("m5n1000")

#----------------------------------------#
### Combine all the simulation results ###
#----------------------------------------#
for (i in stats.list){
  x = get(i)
  colnames(x) = c("bias", "rel.bias","SD","power","coverage","MeanCIWidth","MedianCIWidth","Rel.MeanCIWidth")
  assign(i,x)
} # Mass rename columns for all the combined PM stats data frames. Need to create [stats.list] first (see code line:55)

PM.SimStats.Consistent = bind_rows(m1n100.stats, m1n300.stats, m1n500.stats, m1n1000.stats,
                                   m2n100.stats, m2n300.stats, m2n500.stats, m2n1000.stats,
                                   m3n100.stats, m3n300.stats, m3n500.stats, m3n1000.stats,
                                   m4n100.stats, m4n300.stats, m4n500.stats, m4n1000.stats,
                                   m5n100.stats, m5n300.stats, m5n500.stats, m5n1000.stats) # Combine everything into 1 table.

# Calculate sim run time
fullrun.end = Sys.time() - fullrun.start
fullrun.end # 13.18553 mins

### Export results as CSV
setwd("directory here")
write.csv(PM.SimStats.Consistent, file="PM.SimStats.Consistent.csv", row.names = TRUE)
### Save R Data for replication / backtracking
save.image("PM.SimStats.Consistent_Results.RData")

