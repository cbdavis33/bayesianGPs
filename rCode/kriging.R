rm(list = ls())
cat("\014")

library(plyr)
library(tidyverse)
library(stringr)
library(rstan)
library(reshape2)
library(bayesplot)

# minor edit
# setwd("~/Documents/bayesianGPs/")

rstan_options(auto_write = TRUE)
numCores <- 4                     # Find the number of cores on your machine
options(mc.cores = numCores)      # Then use one less than that for MCMC sampling



sim_data_model <- stan_model(file = 'stanCode/generateKrigingData.stan')

dat_list <- list(N = 1000, alpha = 1, length_scale = 0.15, sigma = sqrt(0.1))
set <- sample(1:dat_list$N, size = 30, replace = F)
draw <- sampling(sim_data_model, iter = 1, algorithm = 'Fixed_param', chains = 1, data = dat_list,
                 seed = 363360090)

samps <- rstan::extract(draw)
plt_df = with(samps,data.frame(x = x[1,], y = y[1,], f = f[1,]))

ggplot(data = plt_df[set,], aes(x=x, y=y)) +
  geom_point(aes(colour = 'Realized data')) +
  geom_line(data = plt_df, aes(x = x, y = f, colour = 'Latent mean function')) +
  theme_bw() + theme(legend.position="bottom") +
  scale_color_manual(name = '', values = c('Realized data'='black','Latent mean function'='red')) +
  xlab('X') +
  ylab('y') +
  ggtitle(paste0('N=',length(set),' from length-scale = 0.15, alpha = 1, sigma = 0.32'))

stan_data <- list(N = length(set), N_pred = dat_list$N - length(set),
                  zeros = rep(0,length(set)), x = samps$x[1,set], y = samps$y[1,set],
                  x_pred = samps$x[1,-set], f_pred = samps$f[1,-set])

comp_gp_mod_lat <- stan_model(file = 'stanCode/kriging.stan')
gp_mod_lat <- sampling(comp_gp_mod_lat, 
                       data = stan_data, 
                       cores = 4, 
                       chains = 4, 
                       iter = 500, 
                       control = list(adapt_delta = 0.999))

samps_gp_mod_lat <- extract(gp_mod_lat)
post_pred <- data.frame(x = stan_data$x_pred,
                        pred_mu = colMeans(samps_gp_mod_lat$f_pred))
plt_df_rt = data.frame(x = stan_data$x_pred, f = t(samps_gp_mod_lat$f_pred))

plt_df_rt_melt = melt(plt_df_rt,id.vars = 'x')

p <- ggplot(data = plt_df[set,], aes(x=x, y=y)) +
  geom_line(data = plt_df_rt_melt, aes(x = x, y = value, group = variable, color = "Posterior draws")) + 
  geom_point(aes(colour = 'Realized data')) +
  geom_line(data = plt_df, aes(x = x, y = f, colour = 'Latent mean function')) +
  geom_line(data = post_pred, aes(x = x, y = pred_mu, colour = 'Posterior mean function')) +
  theme_bw() + 
  theme(legend.position="bottom") +
  scale_color_manual(name = '', values = c('Realized data'='black',
                                           'Latent mean function'='red',
                                           'Posterior draws' = 'blue',
                                           'Posterior mean function' = 'green')) +  
  xlab('X') +
  ylab('y') +
  ggtitle(paste0('N = ',length(set),' from length-scale = 0.15, alpha = 1, sigma = 0.32'))
p


ppc_interval_df <- function(yrep, y) {
  q_95 <- apply(yrep,2,quantile,0.95)
  q_75 <- apply(yrep,2,quantile,0.75)
  q_50 <- apply(yrep,2,median)
  q_25 <- apply(yrep,2,quantile,0.25)
  q_05 <- apply(yrep,2,quantile,0.05)
  mu <- colMeans(yrep)
  df_post_pred <- data.frame(y_obs = y,
                             q_95 = q_95,
                             q_75 = q_75,
                             q_50 = q_50,
                             q_25 = q_25,
                             q_05 = q_05,
                             mu = mu)
  return(df_post_pred)
}
ppc_interval_norm_df <- function(means, sds, y) {
  q_95 <- qnorm(0.95,mean = means, sd = sds)
  q_75 <- qnorm(0.75,mean = means, sd = sds)
  q_50 <- qnorm(0.5,mean = means, sd = sds)
  q_25 <- qnorm(0.25,mean = means, sd = sds)
  q_05 <- qnorm(0.05,mean = means, sd = sds)
  df_post_pred <- data.frame(y_obs = y,
                             q_95 = q_95,
                             q_75 = q_75,
                             q_50 = q_50,
                             q_25 = q_25,
                             q_05 = q_05,
                             mu = means)
  return(df_post_pred)
}
interval_cover <- function(upper, lower, elements) {
  return(mean(upper >= elements & lower <= elements))
}
ppc_full_bayes <- ppc_interval_df(samps_gp_mod_lat$y_pred, samps$y[1,-set])

coverage90 <- interval_cover(upper = ppc_full_bayes$q_95,
                             lower = ppc_full_bayes$q_05,
                             elements = ppc_full_bayes$y_obs)
coverage50 <- interval_cover(upper = ppc_full_bayes$q_75,
                             lower = ppc_full_bayes$q_25,
                             elements = ppc_full_bayes$y_obs)

print(c(coverage90, coverage50))

post <- as.data.frame(gp_mod_lat)
tmp <- select(post, alpha, length_scale, sigma)
bayesplot::mcmc_recover_hist(tmp, true = c(1, 0.15, 0.32),
                             facet_args = list(ncol = 1))

ppc_full_bayes$x <- samps$x[1,-set]
ggplot(data = ppc_full_bayes, aes(x = x, y = y_obs)) +
  geom_ribbon(aes(ymax = q_95, ymin = q_05,alpha=0.5, colour = '90% predictive interval')) + 
  geom_point(data = plt_df[set,], aes(x = x, y = y, colour='Observed data')) + 
  geom_point(data = plt_df[-set,], aes(x = x, y = y, colour='Out-of-sample data')) + 
  theme(legend.position="bottom") + 
  geom_line(data = ppc_full_bayes, aes(x = x, y = mu, colour = 'Posterior predictive mean')) +
  scale_color_manual(name = '', values = c('Observed data' = 'red',
                                           '90% predictive interval' = 'blue',
                                           'Out-of-sample data' = 'black',
                                           'Posterior predictive mean' = 'green')) +
  xlab('X') +
  ylab('y') +
  ggtitle(str_c('Full Bayes PP intervals for N = ',length(set),' from length-scale = 0.15, alpha = 1, sigma = 0.32'))
