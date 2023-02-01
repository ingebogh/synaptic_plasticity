


logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1 + exp(-x)) # inverse logit

# ind: timestep number (integer)
# t: time vector, so that t[ind] = the current time
learning_rule <- function(s1, s2, Ap, Am, taup, taum, t, ind, binsize){
  
  # which timesteps to actually use, do not use all of them from the start
  l <- ind - ceiling(10*taup/binsize)
  inds <- max(l, 1):(ind)
  
  if (s1[ind] == 0 && s2[ind] == 0){
    return(0)
  }
  
  lrm <- 0
  if (s1[ind] != 0){
    lrm <- s1[ind] * sum(s2[inds] * exp((t[inds] - t[ind])/taum))
  }
  
  lrp <- 0
  if (s2[ind] != 0){
    lrp <- s2[ind] * sum(s1[inds] * exp((t[inds] - t[ind])/taup))
  }
  
  return(Ap*lrp - Am*lrm)
  
}

# ind: timestep number (integer)
# t: time vector, so that t[ind] = the current time
learning_rule_mult <- function(s1, s2, w_curr, Ap, Am, taup, taum, t, ind, binsize, w_min, w_max){
  
  # which timesteps to actually use, do not use all of them from the start
  l <- ind - ceiling(10*taup/binsize)
  inds <- max(l, 1):(ind)
  
  if (s1[ind] == 0 && s2[ind] == 0){
    return(0)
  }
  
  lrm <- 0
  if (s1[ind] != 0){
    lrm <- s1[ind] * sum(s2[inds] * exp((t[inds] - t[ind])/taum))
  }
  
  lrp <- 0
  if (s2[ind] != 0){
    lrp <- s2[ind] * sum(s1[inds] * exp((t[inds] - t[ind])/taup))
  }
  
  return(Ap*lrp*(w_max - w_curr) - Am*lrm*(w_curr - w_min))
  
}

### generate data ###
# A and tau: learning rule parameters
# std: standard deviation of the noise of the weight trajectory
# w0: initial value of weight trajectory
# sec: number of seconds we want to generate
# lag: time lag 
# b1 and b2: value of background noise for neuron 1 and 2
# beta and X: value of the p covariate parameters and covariates, as a vector of length p and column matrix of dim. (sec/bs, p)
# binsize: size of time bins
# rule: which learning rule to use when generating data
# w_min and w_max: for multiplicative learning rule only
# seed (= 1 by default): ensure reproducible results, some integer, see set.seed()
generate_data <- function(Ap = 0.005, Am = 0.005*1.05, 
                          taup = 0.02, taum,
                           std = 0.001, 
                           w0 = 1, sec = 1, 
                           lag = 1,
                           b1 = -2, b2 = -2, beta = NULL, X = NULL,
                           binsize = 1/500,
                           rule = "add",
                           w_min = 0, w_max = 1,
                           seed = 1){
  set.seed(seed)
  
  t <- seq(0, sec-binsize, by = binsize)
  if (!is.null(X) && length(t) != nrow(X)) stop("Covariate(s) must have same length as requested time series.")
  
  w0_sign <- if (w0 < 0) -1 else if (w0 > 0) 1 else stop("w0 cannot be 0")
  s1 <- s2 <- rep(0, length(t))
  W <- rep(w0, lag)
  rate_1 <- expit(b1)
  s1[1] <- rbinom(1, 1, rate_1)
  b2beta <- b2 + if (!is.null(beta)) X %*% beta else rep(0, length(t))
  rate_2 <- expit(b2beta[1])
  
  lr <- 0
  
  for (ind in (lag+1):(sec/binsize)){
    
    s1[ind] <- s2[ind] <- 0 
    if (rule == "add") { # additive
      lr <- learning_rule(s1, s2, Ap, Am, taup, taum, t, ind-1, binsize)
    } else if (rule == "mult") { # multiplicative
      lr <- learning_rule_mult(s1, s2, W[ind-1], Ap, Am, taup, taum, t, ind-1, binsize, w_min, w_max)
    } else stop("This rule is not implemented")
    
    step <- W[ind-1] + lr + rnorm(1, 0, std)
    if (rule == "add"){
      W[ind] <- if (w0_sign == 1) max(step, 0) else min(step, 0)
    } else {
      W[ind] <- max(min(step, w_max), w_min)
    }
    
    rate_1[ind] <- expit(b1)
    rate_2[ind] <- expit(W[ind]*s1[ind-lag] + b2beta[ind])
    s2[ind] <- rbinom(1, 1, rate_2[ind])
    s1[ind] <- rbinom(1, 1, rate_1[ind])
    
  }
  
  res <- list(
    s1 = s1,
    s2 = s2,
    t = t,
    lag = lag,
    X = X,
    sec = sec,
    binsize = binsize,
    timesteps = length(s1),
    W = W # true values, only known due to simulated data
  )
  return(res)
  
}

### main function for fitting model ###
# data is a list with the data (must contain: s1, s2, t, lag, X (if covariates), binsize, timesteps) (see run_model.R-file)
# inf_data is a list with algorithm/inference settings/data (see run_model.R-file)
# seed (= 1 by default): ensure reproducible results, some integer, see set.seed()
# store_less = TRUE: only the mean and variance of the weight trajectory is stored to save space and memory for long time series
# initial values: if NULL (default), a value is drawn from the prior
# theta_init: initial value for A and tau (vector of length 2)
# b2_init: initial value of intercept (background noise)
# beta_init: initial value of covariate parameters (vector of same length as number of covariates)
# initial value of w0 is drawn from the prior (one value for each particle) and cannot be set
mh_alg <- function(data, inf_data, seed = 1, store_less = F, theta_init = NULL, b2_init = NULL, beta_init = NULL){
  
  cat("Start.\n")
  
  stopifnot(inf_data$w0_sign %in% c(-1, 1))
  resample <- TRUE
  if (is.null(data$lag)){
    data$lag <- 1
    warning("Setting lag to 1!")
  }
  
  likelihood <- 0
  
  # by default, half is taken as burn-in
  bi <- if (is.null(inf_data$bi)) inf_data$it/2 else inf_data$bi
  
  particle_filter <- if (inf_data$rule == "mult") 
    particle_filter_mult else if (inf_data$rule == "add") 
      particle_filter_add else stop("Not valid rule")
  
  time_start <- Sys.time()
  
  set.seed(seed)
  
  data$start_ind <- 1 + data$lag
  
  w0_m <- inf_data$w0_prior$mean
  w0_s <- w0_s_old <- inf_data$w0_prior$sd
  
  b2_param <- if (is.null(b2_init)) b2_prior(inf_data$b2_prior$mean, inf_data$b2_prior$sd) else b2_init
  b2_sd <- b2_sd_old <- inf_data$b2_prior$sd
  
  n_betas <- ncol(data$X)
  if (is.null(n_betas)) n_betas <- 0
  beta_param <- matrix(0, nrow = inf_data$it, ncol = n_betas)
  beta_sd <- beta_sd_old <- 0
  for (beta_ind in seq_len(n_betas)){
    beta_param[1,beta_ind] <- if (is.null(beta_init)) beta_prior(inf_data$beta_prior$mean[beta_ind], inf_data$beta_prior$sd[beta_ind]) else beta_init[beta_ind]
    beta_sd[beta_ind] <- beta_sd_old[beta_ind] <- inf_data$beta_prior$sd[beta_ind]
  }
  shapes <- shapes_old <- inf_data$shapes_prior
  theta_prior <- if (is.null(theta_init)) rgamma(2, shape = shapes, rate = inf_data$rates_prior) else theta_init
  theta <- theta_prior # store all thetas
  
  tmp <- try(particle_filter(theta_prior, data, inf_data, ind2 = 1, w0_m, w0_s, resample, b2_param[1], beta_param[1,]))
  wp <- tmp$w # note: this is just w
  w0_vec <- w0_vec_new <- tmp$w0_vec
  vp <- vp_new <- tmp$vp
  w0_means <- sum(w0_vec*normalize(vp))
  # summing w and squares of w so we can compute the mean and variance using less memory
  wp_mean2 <- list(w = rep(0, data$timesteps), w2 = rep(0, data$timesteps), no = 0)
  if (!store_less) {
    wp_mean <- matrix(NA, nrow = inf_data$it, ncol = length(wp))
    wp_mean[1,] <- wp
  }
  
  old_log_post <- tmp$log_post
  b_loglik_old <- likelihood_full(data$lag, wp, data$s1, data$s2, b2_param[1], beta_param[1,], data$X)
  
  log_post_vec <- old_log_post
  
  w0_s_vec <- w0_s
  
  for (ind in 2:inf_data$it){

    # first update fixed effect parameters one by one
    # b2_sd is the value we use now, while b2_sd_old is the st.dev. used for the previous accepted sample
    tmp <- propose_b2_value(data, wp, b2_param[ind-1], b_loglik_old, b2_sd, b2_sd_old, inf_data, beta_param[ind-1,])
    b2_param[ind] <- tmp$b2
    b_loglik_old <- tmp$loglik
    b2_sd_old <- tmp$b2_sd 
    
    for (beta_ind in seq_len(n_betas)){
      
      tmp <- propose_beta_value(data, wp, beta_ind, b_loglik_old, beta_param[ind-1,], beta_sd, beta_sd_old, inf_data, b2_param[ind])
      beta_param[ind, beta_ind] <- tmp$betas[beta_ind]
      b_loglik_old <- tmp$loglik
      beta_sd_old <- tmp$beta_sd
      
    }
    
    if (ind %% inf_data$Usim == 0){
      b2_sd <- adjust_b2_sd(b2_param, b2_sd, inf_data)$b2_sd
      for (beta_ind in seq_len(n_betas)){
        beta_sd[beta_ind] <- adjust_b2_sd(beta_param[,beta_ind], b2_sd[beta_ind], inf_data)$b2_sd
      } 
    }
    
    if (ind %% inf_data$Usim == 0){
      tmp <- adjust_variance(theta, shapes, inf_data, ind)
      shapes <- tmp$shapes
      theta_next <- tmp$theta
    } else {
      theta_next <- rgamma(2, shape = shapes, scale = theta_prior/shapes)
    }
    
    if (any(theta_next < 0)) {
      stop("theta < 0")
    }
    
    w0_means[ind] <- sum(w0_vec*normalize(vp))
    w0_m <- w0_means[ind]
    
    if (ind %% inf_data$Usim2 == 0){
      w0_s <- adjust_w0_sd(w0_means, w0_s, inf_data, ind)$w0_s
      w0_s_vec <- c(w0_s_vec, w0_s)
    }
    
    tmp <- try(particle_filter(theta_next, data, inf_data, ind2 = ind, w0_m, w0_s, resample, b2_param[ind], beta_param[ind,]))
    if (class(tmp) == "try-error"){ # crashing
      return(list(theta = theta, log_post = log_post_vec, 
                  w = if (!store_less) wp_mean else NULL,
                  w22 = wp_mean2,
                  w0_mean = w0_means,
                  b2 = b2_param,
                  beta_param = beta_param,
                  w0_s_vec = w0_s_vec,
                  inf_data = inf_data,
                  time = Sys.time() - time_start
      ))
    }
    wp_new <- tmp$w
    w0_vec_new <- tmp$w0_vec
    vp_new <- tmp$vp
    new_log_post <- tmp$log_post
    w0_prop <- sum(w0_vec_new*normalize(vp_new))
    
    tmp <- scaled2_spike_prob_log(old_log_post, new_log_post)
    prob_old <- tmp[1]
    prob_next <- tmp[2]
    
    r <- ratio_log(prob_old, prob_next, shapes, shapes_old, theta_next, theta_prior, w0_prop, w0_m, w0_s, w0_s_old, inf_data)
    
    if (runif(1) < r){
      theta_choice <- theta_next
      old_log_post <- new_log_post
      wp <- wp_new
      w0_vec <- w0_vec_new
      vp <- vp_new
      shapes_old <- shapes
      w0_s_old <- w0_s
    } else {
      theta_choice <- theta_prior
    }
    theta <- rbind(theta, theta_choice)
    theta_prior <- theta_choice
    log_post_vec <- c(log_post_vec, old_log_post)
    
    if (ind > bi) {
      wp_mean2$w <- wp_mean2$w + wp
      wp_mean2$w2 <- wp_mean2$w2 + wp^2
      wp_mean2$no <- wp_mean2$no + 1
    }
    if (!store_less) {
      wp_mean[ind,] <- wp
    }
    
    # compute likelihood
    likelihood <- c(likelihood, likelihood_full(data$lag, wp, data$s1, data$s2, b2_param[ind], beta_param[ind,], data$X))
    
    if (ind %% round(inf_data$it/10) == 0) {
      cat(paste0(ind, ", ", format(Sys.time() - time_start), "\n"))
      flush.console()
    }
    
    if (ind == round(inf_data$it/10)){
      oo <- Sys.time() - time_start
      cat(
        paste0("First ", round(inf_data$it/10), " iterations took ", format(oo),
               ". Estimating that ", inf_data$it, " iterations takes ", 
               round(as.numeric(strsplit(format(oo), " ")[[1]][1])*inf_data$it/(round(inf_data$it/10))), " ",
               strsplit(format(oo), " ")[[1]][2],
               ".\nThis depends on among others mixing and memory.\n")
      )
    }
    
  }
  
  cat(paste0("Total time used: ", format(Sys.time() - time_start), ".\n"))
  
  return(list(theta = theta, log_post = log_post_vec, 
              w = if (!store_less) wp_mean else NULL,
              w_summary = wp_mean2,
              w0_mean = w0_means,
              b2 = b2_param,
              beta_param = beta_param,
              likelihood = likelihood,
              inf_data = inf_data,
              time = Sys.time() - time_start
  ))
  
}

# the rest of the functions are not described, and are performing operations for the main algorithm

particle_filter_add <- function(theta, data, inf_data, ind2 = -1, w0_m, w0_s, resample, b2_param, beta_param){
  
  t <- data$t
  start_ind <- data$start_ind
  
  t <- data$t
  
  wp_curr <- w0_prior(inf_data$P, w0_m, w0_s, inf_data$w0_sign)
  wp <- matrix(0, ncol = length(t), nrow = inf_data$P)
  wp[,start_ind-1] <- wp_curr
  index_list <- list()
  
  vp <- rep(1, inf_data$P)
  # we do not have a likelihood value for the first time-step(s), because it needs s_{1, t = 0}, and thus
  # we start at t = lag
  log_posterior <- 0 
  
  b2beta <- b2_param + if (!is.null(data$X)) data$X %*% beta_param else rep(0, length(t))
  
  for (ind in start_ind:length(t)){
    v_normalized <- normalize(vp)
    if (resample) {
      perplexity <- perplexity_func(v_normalized)
      if (is.na(perplexity)) {
        stop(paste0("Iteration ", ind2, ", timestep ", ind))
      }
      if (perplexity < 0.66){
        tmp <- resampling_full2(v_normalized, wp_curr)
        wp_curr <- tmp$wp
        index_list <- c(index_list, list(list(indexes = tmp$indexes, ind = ind - 1))) # -1 because we resample previous w_curr
        vp <- rep(1/inf_data$P, inf_data$P)
        v_normalized <- normalize(vp)
      }
    }
    # computing w_t = w_{t-1} + lr(t-1) + noise
    lr <- learning_rule(data$s1, data$s2, 
                         theta[1], theta[1]*1.05, theta[2], theta[2], 
                         t, ind-1, data$binsize)

    step <- wp_curr + lr + rnorm(inf_data$P, 0, inf_data$infstd)
    if (inf_data$w0_sign == 1) step[step < 0] <- 0 else step[step > 0] <- 0
    wp_curr <- step
    ls <- likelihood_step(data$s1[ind-data$lag], data$s2[ind], wp_curr, b2beta[ind]) # using w_t to compute likelihood
    wp[,ind] <- wp_curr
    vp <- ls*v_normalized
    
    log_posterior <- log_posterior + log(sum(vp)/inf_data$P)
    
  }
  
  # do all the resampling now
  
  for (ind in rev(seq_along(index_list))){
    if (ind == length(index_list)){
      index_list[[ind]]$indexes2 <- index_list[[ind]]$indexes
      next
    }
    index_list[[ind]]$indexes2 <- index_list[[ind]]$indexes[index_list[[ind+1]]$indexes2]
  }
  
  ind1 <- 1
  for (ind in seq_along(index_list)){
    ind3 <- index_list[[ind]]$ind
    wp[,ind1:ind3] <- wp[index_list[[ind]]$indexes2,ind1:ind3]
    ind1 <- ind3 + 1
  }
  
  w_mean <- colSums(wp*normalize(vp))

  return(list(w = w_mean, log_post = log_posterior, w0_vec = wp[,start_ind-1], vp = vp))
  
}

particle_filter_mult <- function(theta, data, inf_data, ind2 = -1, w0_m, w0_s, resample, b2_param, beta_param){
  
  t <- data$t
  start_ind <- data$start_ind
  
  t <- data$t
  
  wp_curr <- w0_prior(inf_data$P, w0_m, w0_s, inf_data$w0_sign)
  wp <- matrix(0, ncol = length(t), nrow = inf_data$P)
  wp[,start_ind-1] <- wp_curr
  index_list <- list()
  
  vp <- rep(1, inf_data$P)
  # we do not have a likelihood value for the first time-step(s), because it needs s_{1, t = 0}, and thus
  # we start at t = lag
  log_posterior <- 0 
  
  b2beta <- b2_param + if (!is.null(data$X)) data$X %*% beta_param else rep(0, length(t))
  
  for (ind in start_ind:length(t)){
    v_normalized <- normalize(vp)
    if (resample) {
      perplexity <- perplexity_func(v_normalized)
      if (is.na(perplexity)) {
        stop(paste0("Iteration ", ind2, ", timestep ", ind))
      }
      if (perplexity < 0.66){
        tmp <- resampling_full2(v_normalized, wp_curr)
        wp_curr <- tmp$wp
        index_list <- c(index_list, list(list(indexes = tmp$indexes, ind = ind - 1))) # -1 because we resample previous w_curr
        vp <- rep(1/inf_data$P, inf_data$P)
        v_normalized <- normalize(vp)
      }
    }
    # computing w_t = w_{t-1} + lr(t-1) + noise
    lr <- learning_rule_mult(data1$s1, data$s2, wp_curr,
                             theta[1], theta[1]*1.05, theta[2], theta[2], 
                             t, ind-1, data$binsize, w_min = inf_data$w_min, w_max = inf_data$w_max)

    step <- wp_curr + lr + rnorm(inf_data$P, 0, inf_data$infstd)
    if (inf_data$w0_sign == 1) step[step < 0] <- 0 else step[step > 0] <- 0
    wp_curr <- step
    ls <- likelihood_step(data$s1[ind-data$lag], data$s2[ind], wp_curr, b2beta[ind]) # using w_t to compute likelihood
    wp[,ind] <- wp_curr
    vp <- ls*v_normalized
    
    log_posterior <- log_posterior + log(sum(vp)/inf_data$P)
  }
  
  # do all the resampling now
  
  for (ind in rev(seq_along(index_list))){
    if (ind == length(index_list)){
      index_list[[ind]]$indexes2 <- index_list[[ind]]$indexes
      next
    }
    index_list[[ind]]$indexes2 <- index_list[[ind]]$indexes[index_list[[ind+1]]$indexes2]
  }
  
  ind1 <- 1
  for (ind in seq_along(index_list)){
    ind3 <- index_list[[ind]]$ind
    wp[,ind1:ind3] <- wp[index_list[[ind]]$indexes2,ind1:ind3]
    ind1 <- ind3 + 1
  }
  
  w_mean <- colSums(wp*normalize(vp))

  return(list(w = w_mean, log_post = log_posterior, w0_vec = wp[,start_ind-1], vp = vp))
  
}

ratio_log <- function(prob_old, prob_next, shapes, shapes_old, theta_next, theta_prior, w0_prop, w0_old, w0_sd, w0_sd_old, inf_data){
  
  spike_prob_ratio <- prob_next - prob_old
  prior_ratio <- 0
  proposal_ratio <- 0

  # A and tau (= theta)
  for (ind in 1:2){
    prior_ratio <- prior_ratio +
      dgamma(theta_next[ind], shape = inf_data$shapes_prior[ind], rate = inf_data$rates_prior[ind], log = T) -
      dgamma(theta_prior[ind], shape = inf_data$shapes_prior[ind], rate = inf_data$rates_prior[ind], log = T)
    proposal_ratio <- proposal_ratio +
      dgamma(theta_prior[ind], shape = shapes[ind], scale = theta_next[ind]/shapes[ind], log = T) -
      dgamma(theta_next[ind], shape = shapes_old[ind], scale = theta_prior[ind]/shapes_old[ind], log = T)
  }

  if (inf_data$w0_sign > 0){
    # w0
    prior_ratio <- prior_ratio +
      sum( # if we have several initial values (for when time series is split)
        log(truncnorm::dtruncnorm(w0_prop, 0, Inf, inf_data$w0_prior$mean, inf_data$w0_prior$sd)) -
          log(truncnorm::dtruncnorm(w0_old, 0, Inf, inf_data$w0_prior$mean, inf_data$w0_prior$sd))
      )
    
    proposal_ratio <- proposal_ratio +
      sum( # if we have several initial values (for when time series is split)
        log(truncnorm::dtruncnorm(w0_old, 0, Inf, w0_prop, w0_sd)) -
          log(truncnorm::dtruncnorm(w0_prop, 0, Inf, w0_old, w0_sd_old))
      )
  } else {
    # w0
    prior_ratio <- prior_ratio +
      sum( # if we have several initial values (for when time series is split)
        log(truncnorm::dtruncnorm(w0_prop, -Inf, 0, inf_data$w0_prior$mean, inf_data$w0_prior$sd)) -
          log(truncnorm::dtruncnorm(w0_old, -Inf, 0, inf_data$w0_prior$mean, inf_data$w0_prior$sd))
      )
    
    proposal_ratio <- proposal_ratio +
      sum( # if we have several initial values (for when time series is split)
        log(truncnorm::dtruncnorm(w0_old, -Inf, 0, w0_prop, w0_sd)) -
          log(truncnorm::dtruncnorm(w0_prop, -Inf, 0, w0_old, w0_sd_old))
      )
  }
  
  
  return(exp(spike_prob_ratio + prior_ratio + proposal_ratio))
  
}

scaled2_spike_prob_log <- function(old, new) c(old - min(old, new), new - min(old, new))

adjust_variance <- function(theta, shapes, inf_data, ind_it){
  
  var_new <- c(0, 0)
  h_temp <- inf_data$Hsim
  while (any(var_new == 0)){
    inds <- if (h_temp >= nrow(theta)) 1:nrow(theta) else -c(1:(nrow(theta)-h_temp))
    thetas <- theta[inds,]
    means <- colMeans(thetas)
    var_new <- apply(thetas, 2, var) * (2.4/sqrt(2))^2
    h_temp <- h_temp + 50
    if (h_temp > inf_data$it){
      return(list(shapes = shapes, theta = rgamma(2, shape = shapes, scale = theta[nrow(theta),]/shapes)))
    }
  }
  new_shapes <- means^2/var_new
  proposal <- rgamma(2, shape = new_shapes, scale = theta[nrow(theta),]/new_shapes)
  
  return(list(shapes = new_shapes, theta = proposal))
}

adjust_variance_std <- function(par, shape, inf_data){
  
  var_new <- 0
  h_temp <- inf_data$Hsim
  while (var_new == 0){
    inds <- if (h_temp >= length(par)) 1:length(par) else -c(1:(length(par)-h_temp))
    mean_new <- mean(par[inds])
    var_new <- var(par[inds]) * 2.4^2
    h_temp <- h_temp + 50
    if (h_temp > inf_data$it){
      return(list(shape = shape))
    }
  }
  new_shape <- mean_new^2/var_new
  
  return(list(shape = new_shape))
  
}

adjust_w0_sd <- function(w0_means, w0_s, inf_data, ind_it){
  
  var_new <- 0
  h_temp <- inf_data$Hsim2
  while (var_new == 0){
    inds <- if (h_temp >= length(w0_means)) 1:length(w0_means) else -c(1:(length(w0_means)-h_temp))
    w0s <- w0_means[inds]
    var_new <- var(w0s) * 2.4^2
    h_temp <- h_temp + 50
    if (h_temp > inf_data$it){
      return(list(w0_s = w0_s))
    }
  }
  w0_sd_new <- sqrt(var_new)
  
  # min_w0_sd is the smallest standard deviation we allow
  return(list(w0_s = max(w0_sd_new, inf_data$w0_prior$min_w0_sd)))
  
}

w0_prior <- function(P, w0_m = 1, w0_s = 1, w0_sign = 1){
  if (w0_sign > 0){
    truncnorm::rtruncnorm(P, 0, Inf, w0_m, w0_s)
  } else {
    truncnorm::rtruncnorm(P, -Inf, 0, w0_m, w0_s)
  }
}

b2_prior <- function(b2, b2_sd) rnorm(1, b2, b2_sd)

propose_b2_value <- function(data, w_it, b2_old, loglik_old, b2_sd, b2_sd_old, inf_data, beta_param){
  
  b2_prop <- b2_prior(b2_old, b2_sd)
  
  loglik_new <- likelihood_full(data$lag, w_it, data$s1, data$s2, b2_prop, beta_param, data$X)
  
  # likelihood
  acc_rate <- loglik_new - loglik_old
  
  # prior
  acc_rate <- acc_rate +
    dnorm(b2_prop, inf_data$b2_prior$mean, inf_data$b2_prior$sd, log = T) -
    dnorm(b2_old, inf_data$b2_prior$mean, inf_data$b2_prior$sd, log = T)
  
  # proposal
  acc_rate <- acc_rate +
    dnorm(b2_old, b2_prop, b2_sd, log = T) - 
    dnorm(b2_prop, b2_old, b2_sd_old, log = T)
  
  if (runif(1) < exp(acc_rate)){ 
    return(list(b2 = b2_prop, loglik = loglik_new, b2_sd = b2_sd))
  } else {
    return(list(b2 = b2_old, loglik = loglik_old, b2_sd = b2_sd_old))
  }
  
}

adjust_b2_sd <- function(b2, b2_sd, inf_data){
  
  var_new <- 0
  h_temp <- inf_data$Hsim
  while (var_new == 0){
    inds <- if (h_temp >= length(b2)) 1:length(b2) else -c(1:(length(b2)-h_temp))
    var_new <- var(b2[inds]) * 2.4^2
    h_temp <- h_temp + 50
    if (h_temp > inf_data$it){
      return(list(b2_sd = b2_sd))
    }
  }
  b2_sd_new <- sqrt(var_new)
  
  return(list(b2_sd = b2_sd_new))
  
}

likelihood_full <- function(lag, w, s1, s2, b2, betas = NULL, X_mat = NULL){
  
  # Note: Do not have likelihood for t = lag and less, because that requires s1-data for t <= 0
  eta <- w[(lag+1):length(s1)]*s1[1:(length(s1)-lag)] + b2
  if (length(betas) > 0){
    eta <- eta + c(X_mat %*% betas)[(lag+1):length(s1)]
  }
  mu <- expit(eta)
  loglik <- sum(s2[-c(1:lag)]*log(mu) + (1-s2[-c(1:lag)]) * log1p(-mu))
  if (any(mu == 1)) loglik <- sum(log(mu^s2[-c(1:lag)]) + log((1-mu)^(1-s2[-c(1:lag)])))
  if (is.na(loglik)) {
    loglik <- -Inf
    warning(paste0("b2 value ", b2, " gives NaN in loglik, setting loglik to -Inf"))
  }
  
  return(loglik)
  
}

normalize <- function(x) x/sum(x)

perplexity_func <- function(x){
  # for x == 0, 0*log(0) is (in the limit and for Shannon entropy) 
  # taken to be 0, which we also get when x = 1, so setting 0 to 1
  x[x == 0] <- 1
  h <- -sum(x*log(x))
  return(exp(h)/length(x))
}

resampling_full2 <- function(v_norm, wp){
  indexes <- sample(1:length(v_norm), size = length(v_norm), replace = TRUE, prob = v_norm)
  return(list(wp = wp[indexes], indexes = indexes))
}

likelihood_step <- function(s1prev, s2next, wcurr, b2beta){
  
  if (s2next == 0){
    return(exp(-log1p(exp(wcurr*s1prev + b2beta))))
  } else if (s2next == 1){
    return(exp(-log1p(exp(-(wcurr*s1prev + b2beta)))))
  } else stop("Observations must be 0 or 1")
  
}

beta_prior <- function(beta, beta_sd) rnorm(1, beta, beta_sd)

propose_beta_value <- function(data, w_it, beta_num, loglik_old, all_betas_old, beta_sd, beta_sd_old, inf_data, b2_param){
  
  beta_old <- all_betas_old[beta_num]
  beta_prop <- beta_prior(beta_old, beta_sd)
  all_betas_prop <- all_betas_old
  all_betas_prop[beta_num] <- beta_prop
  
  loglik_new <- likelihood_full(data$lag, w_it, data$s1, data$s2, b2_param, all_betas_prop, data$X)
  
  # likelihood
  acc_rate <- loglik_new - loglik_old
  
  # prior
  acc_rate <- acc_rate +
    dnorm(beta_prop, inf_data$beta_prior$mean[beta_num], inf_data$beta_prior$sd[beta_num], log = T) -
    dnorm(beta_old, inf_data$beta_prior$mean[beta_num], inf_data$beta_prior$sd[beta_num], log = T)
  
  # proposal
  acc_rate <- acc_rate +
    dnorm(beta_old, beta_prop, beta_sd[beta_num], log = T) -
    dnorm(beta_prop, beta_old, beta_sd_old[beta_num], log = T)
  
  if (runif(1) < exp(acc_rate)){ # min(acc_rate, 1)
    return(list(betas = all_betas_prop, loglik = loglik_new, beta_sd = beta_sd))#, a = acc_rate, val = b2_prop, val2 = loglik_new))
  } else {
    return(list(betas = all_betas_old, loglik = loglik_old, beta_sd = beta_sd_old))#, a = acc_rate, val = b2_prop, val2 = loglik_new))
  }  
  
}











