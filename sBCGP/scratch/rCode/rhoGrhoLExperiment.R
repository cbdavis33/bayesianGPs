rm(list = ls())
cat("\014")

library(plyr)
library(tidyverse)
library(stringr)
library(rstan)
library(reshape2)
library(bayesplot)

rstan_options(auto_write = TRUE)
numCores <- 4                     # Find the number of cores on your machine
options(mc.cores = numCores)      # Then use one less than that for MCMC sampling

D <- 3
alphaW <- 6
betaW <- 4
alphaRhoG <- c(1, 1, 1)
betaRhoG <- c(1, 1, 1)
alphaRhoL <- c(1, 8, 3)
betaRhoL <- c(1, 2, 6)

stanData <- list(d = D, alphaW = alphaW, betaW = betaW, alphaRhoG = alphaRhoG, betaRhoG = betaRhoG,
                 alphaRhoL = alphaRhoL, betaRhoL = betaRhoL)

# stanModel <- stan_model(file = 'sBCGP/scratch/stanCode/rhoGRhoL1.stan')

samples <- sampling(stanModel, 
                    data = stanData, 
                    cores = 1, 
                    chains = 4, 
                    iter = 1000, 
                    control = list(adapt_delta = 0.90))

parsToMonitor <- c("w", "wRaw", "rhoG", "rhoL", "rhoLRaw")
print(samples, pars = parsToMonitor, digits = 4)

samps <- extract(samples)

wRaw <- seq(0, 1, .005)
fwRaw <- dbeta(wRaw, alphaW, betaW)

par(mfrow = c(3,1))
hist(samps$w)
hist(samps$wRaw)
plot(w, fw, type = 'l')

print(round(alphaW/(alphaW + betaW),4))
print(round(sqrt(alphaW*betaW/((alphaW + betaW)^2 * (alphaW + betaW + 1))), 4))
      