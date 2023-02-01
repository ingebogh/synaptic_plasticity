
# load all functions for inference
source("all_functions.R")

# generating data ####

# setting values
# note: we have chosen to use high spike rates and a short dataset to increase speed of the inference in this example
A <- 0.005
tau <- 0.02
bs <- 0.002
w0 <- 0.6
spike_rate1 <- 60
spike_rate2 <- 60
b1 <- logit(spike_rate1*bs)
b2 <- logit(spike_rate2*bs)
sec <- 5 
beta_1 <- 0.5
X <- matrix(cos((1:(sec/bs))/100)) # some covariate

# simulating data
data0 <- generate_data(Ap = A, Am = A*1.05, taup = tau, taum = tau, std = 0.001, w0 = w0, 
                       sec = sec, lag = 1,
                       b1 = b1, b2 = b2, beta = beta_1, X = X,
                       binsize = bs, rule = "add", seed = 1)



# creating a list with the data needed for inference
# to use your own data, create a list with the same elements
data1 <- list(
  s1 = data0$s1,
  s2 = data0$s2,
  t = data0$t,
  lag = data0$lag,
  X = data0$X, # if covariates, if not do not include this
  binsize = data0$binsize,
  timesteps = data0$timesteps
)

# list of settings for the inference ####
inf_dat <- list(
  shapes_prior = c(2, 4), # shape parameters for prior on A [1] and tau [2]
  rates_prior = c(200, 100), # rate parameters for prior on A [1] and tau [2]
  it = 2000, # number of iterations in total
  bi = 500, # number of samples to use as burn-in (defaults to half the samples)
  P = 200, # number of particles in the particle filter
  Usim = 100, # how often do we adjust the variance of A, tau, b2 and beta
  Usim2 = 100, # how often do we adjust the variance of w0
  Hsim = 100, # how many samples do we go back when adjusting the variance of A, tau, b2 and beta
  Hsim2 = 100, # how many samples do we go back when adjusting the variance of w0
  w0_prior = list(mean = sign(w0), sd = 5, min_w0_sd = 0.1), # the prior for w0
  w0_sign = sign(w0), # the sign of the weight trajectory (and thus w0)
  b2_prior = list(mean = 0, sd = 10), # the prior of b2
  beta_prior = list(mean = c(0), sd = c(5)), # the prior of the covariate parameters
  w_min = if (w0 > 0) 0 else -1, # parameter for the multiplicative learning rule
  w_max = if (w0 > 0) 1 else 0 # parameter for the multiplicative learning rule
)
inf_dat$infstd <- inf_dat$shapes_prior[1]/inf_dat$rates_prior[1]/3 # fixed value of the noise of the weight trajectory
inf_dat$rule <- "add" # the learning rule we choose

# fit model ####
res <- mh_alg(data1, inf_dat, store_less = T,
              b2_init = inf_dat$b2_prior$mean, beta_init = NULL, 
              theta_init = inf_dat$shapes_prior/inf_dat$rates_prior)


# some results #### 

A_post <- res$theta[(inf_dat$bi+1):inf_dat$it, 1]
tau_post <- res$theta[(inf_dat$bi+1):inf_dat$it, 2]
# computing mean and standard deviation of weight trajectory
w_post <- res$w_summary$w/res$w_summary$no
w_sd_post <- sqrt((res$w_summary$w2/res$w_summary$no - (res$w_summary$w/res$w_summary$no)^2)*(res$w_summary$no/(res$w_summary$no-1)))

# graphs

library(ggplot2)
ggplot(data.frame(A = A_post), aes(x = A, y = ..density..)) +
  geom_histogram(aes(fill = "Posterior"), alpha = 0.4, col = NA, bins = 30) + 
  geom_density(aes(col = "Posterior"), fill = NA) +
  geom_function(fun = dgamma, args = list(shape = inf_dat$shapes_prior[1], rate = inf_dat$rates_prior[1]),
                inherit.aes = F, mapping = aes(col = "Prior")) +
  scale_color_manual(values = c("Prior" = "blue", "Posterior" = "black"), aesthetics = c("color", "fill"),
                     name = element_blank()) +
  ggtitle("Prior and posterior of A")
ggplot(data.frame(tau = tau_post), aes(x = tau, y = ..density..)) +
  geom_histogram(aes(fill = "Posterior"), alpha = 0.4, col = NA, bins = 30) + 
  geom_density(aes(col = "Posterior"), fill = NA) +
  geom_function(fun = dgamma, args = list(shape = inf_dat$shapes_prior[2], rate = inf_dat$rates_prior[2]),
                inherit.aes = F, mapping = aes(col = "Prior")) +
  scale_color_manual(values = c("Prior" = "blue", "Posterior" = "black"), aesthetics = c("color", "fill"),
                     name = element_blank()) +
  ggtitle("Prior and posterior of tau")

# computing the additive learning rule for a value of A and tau
comp_add_lr <- function(A, tau, dt){
  Ap <- A
  Am <- A*1.05
  lr <- rep(0, length(dt))
  lr[dt <= 0] <- lr[dt <= 0] - Am*exp(dt[dt <= 0]/tau)
  lr[dt >= 0] <- lr[dt >= 0] + Ap*exp(-dt[dt >= 0]/tau)
  return(data.frame(dt = dt, lr = lr))
}
ggplot(comp_add_lr(median(A_post), median(tau_post), seq(-0.2, 0.2, length.out = 1001)), aes(x = dt, y = lr)) +
  geom_line() +
  ggtitle("Learning rule from posterior medians")

ggplot(data.frame(t = data1$t, w_truth = data0$W, w = w_post, w_sd = w_sd_post), aes(x = t)) +
  geom_ribbon(aes(ymin = w-w_sd, ymax = w+w_sd, fill = "Posterior"), alpha = 0.4) +
  geom_line(aes(y = w, col = "Posterior")) + 
  geom_line(aes(y = w_truth, col = "Truth")) +
  scale_color_manual(values = c("Truth" = "red", "Posterior" = "black"), aesthetics = c("color", "fill"),
                     name = element_blank()) +
  ggtitle("Posterior and true weight trajectory")





