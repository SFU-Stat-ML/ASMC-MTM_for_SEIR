# List of seeds
seeds <- c(...)

# Iterate over each seed
for (seed in seeds) {
  # Set the seed for reproducibility
  set.seed(seed)
  
  # bisection function
  
  # - recursive implementation of bisection algorithm
  bisection <- function(low, high, W, u, eta){
    mid <- (low + high)/2
    f.low <- rCESS(W, u, low, eta)
    f.mid <- rCESS(W, u, mid, eta)
    f.high <- rCESS(W, u, high, eta)
    
    if(f.low * f.high > 0)
      stop('Invalid endpoint for bisection.')
    
    try({if( low >= high )
      stop('bisection overlap')
      
      if((abs(f.mid) < 1e-10)||((high - low)/2<1e-10))
        return(mid)
      if((f.low * f.mid) < 0)
        return(bisection(low, mid, W, u, eta))
      if((f.high * f.mid) < 0)
        return(bisection(mid, high, W, u, eta))
    })
    
    stop('bisection flawed')
  }
  
  # log-sum-exponential evaluation of log(sum(w)) (Avoids numerical underflow and overflow)
  logsum <- function(logw){
    logmax = max(logw)
    log(sum(exp(logw-logmax))) + logmax
  }
  
  rCESS <- function(W, u, a, eta) {
    logw <- a * u          # weight update
    exp(2*logsum(log(W)+logw) - logsum(log(W)+2*logw)) - eta
  }
  
  # effective sample size (Conditional effective sample size)
  rESS <- function(logW){
    K <- length(logW) # number of particles
    logWmax <- max(logW)
    logRESS <- -(2*logWmax + log(sum(exp(2*logW-2*logWmax)))) - log(K)
    return(exp(logRESS))
  } 
  
  systematic_resample <- function(W){
    K <- length(W)
    U <- runif(1,0,1/K) + 0:(K-1)/K # these evenly spaced points determine which weights are selected
    W.sum <- cumsum(W)
    N1 <- rep(NA,K)
    j <- 1
    
    # re-sampling loop
    for( i in 1:K )
    {
      found = F
      while( !found )
      {
        if( U[i]>W.sum[j] )
          j <- j+1
        else
          found = T
      }
      N1[i] <- j 
    }
    return( N1 ) 
  }
  
  # ASMC function ----------------------------------------------------------------
  
  #' Annealed Sequential Monte Carlo
  #'
  #' Implementation of the annealed sequential Monte Carlo algorithm with adaptively determined annealing scheme
  #' @param K Number of particles
  #' @param theta_cov Either a proposal covariance matrix or a string "adaptive" to use an adaptive proposal scheme
  #' @param tuning_param List of tuning parameters for the ASMC algorithm: 1) eps - re-sampling threshold,
  #'  2) eta - target rCESS for adaptive re-sampling scheme, 3) alpha - annealing scheme. Only supply one of eta or alpha, not both.
  #' @param init Initialization function for the parameters
  #' @param data Matrix of data used in the likelihood
  #' @param likelihood Function to evaluate log-likelihood. Takes data and parameters as arguments.
  #' @param prior Function to evaluate the log-prior. Takes parameters as arguments.
  #' @param reference Function to evaluate the log-reference distribution. Takes parameters as arguments.
  #' @returns List of particles, weights, marginal log-likelihood, and effective sample size
  #' @export
  #' @examples
  #' 
  asmc <- function(K, theta_cov, tuning_param, init, data, likelihood, prior, reference){
    
    # initialize storage list
    alpha <- list() # a list for tempering parameters
    ESS <- list()   # a list for ESS
    logZ <- list()  # list for log-normalizing constant
    particles <- list()
    W <- list()
    
    # check if adaptive annealing scheme is being used
    adaptive_alpha <- "eta" %in% names(tuning_param) 
    
    r <- 1                
    logZ[[r]] <- 0 # log-normalizing constant for the 1st iteration is set to 0.
    
    # initialize the annealing parameter 
    if(adaptive_alpha){
      alpha[[r]] <- 0    # tempering parameters
      alphaDiff <- 0     # difference b/w two successive tempering parameters
    } else{
      alpha <- as.list(tuning_param$alpha)
    }
    
    particles[[r]] <- lapply(1:K, function(k) list(theta = init(), accept = F))
    
    weighted_cov <- diag(d)
    W[[r]] <- rep(1/K,K)
    logW <- log(W[[r]])
    
    # un-normalized weights
    w <- rep(1,K) # initially all the un-normalized weights are set to 1.
    logw <- rep(0,K)
    
    # begin iterative process
    while( alpha[[r]]<1 )   # repeat this step if the tempering parameter is less than 1
    {
      cat("iteration:",r,"\n") # prints the current iteration number
      r <- r+1  # increment iteration
      # evaluate the log-likelihood for updating alpha
      u <- rep(0, K)   # incremental log-importance weights
      
      # evaluate the log-likelihood for each particle using current data and previous theta.
      u <- sapply(1:K, function(k){
        logL <- likelihood(data, particles[[r-1]][[k]]$theta)
        return(logL)
      })
      
      if(adaptive_alpha){
        alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$eta)
        alpha[[r]] <- alpha[[r-1]] + alphaDiff
      } else{
        alphaDiff <- alpha[[r]] - alpha[[r-1]]
      }
      
      cat("annealing parameter:",alpha[[r]],"\n") 
      
      if( alpha[[r]]>1 ){
        alpha[[r]] <- 1
        alphaDiff <- 1-alpha[[r-1]] # alphaDiff is recalculated
      }
      
      # Updated MCMCmove Function
      MCMCmove <- function(particle, theta_cov){
        # Adaptive MTM
        MH <- function(theta_old, data, likelihood, prior, alpha){
          theta_proposals <- matrix(NA, nrow = J, ncol = length(theta_old))
          
          for (j in 1:J) {
            if(r <= 2 * d){
              theta_cov <- 0.1^2 * diag(d)/d
              theta_proposals[j, ] <- rtmvnorm(1, theta_old, theta_cov, lower = rep(0, d))
            }
            else{
              if(runif(1) < b){
                theta_cov <- 2.38^2 * weighted_cov / d
                theta_proposals[j, ] <- rtmvnorm(1, theta_old, theta_cov, lower = rep(0, d))
              }
              else{
                theta_cov <- 0.1^2 * diag(d) / d
                theta_proposals[j, ] <- rtmvnorm(1, theta_old, theta_cov, lower = rep(0, d))
              }
            }
          }
          
          # Compute likelihoods and priors for all proposals
          likelihoods <- alpha * apply(theta_proposals, 1, likelihood, data = data)
          priors <- apply(theta_proposals, 1, prior)
          
          # Compute weights (unnormalized probabilities) for all proposals
          logwyx <- likelihoods + priors
          
          # Normalize weights to get probabilities
          probabilities <- exp((logwyx) - log(sum(exp(logwyx))))
          
          # Step 2: Select Y from the proposals based on the computed probabilities
          index <- sample(1:J, size = 1, prob = probabilities)
          theta_new <- theta_proposals[index, ]
          
          # Generate k-1 new proposals from T(theta_new, .) and compute weights for the new set
          new_proposals <- matrix(NA, nrow = J - 1, ncol = length(theta_old))
          
          for (j in 1:(J - 1)) {
            if(r <= 2 * d){
              theta_cov <- 0.1^2 * diag(d)/d
              new_proposals[j, ] <- rtmvnorm(1, theta_new, theta_cov, lower = rep(0, d)) 
            }
            else{
              if(runif(1) < b){
                theta_cov <- 2.38^2 * weighted_cov / d
                new_proposals[j, ] <- rtmvnorm(1, theta_new, theta_cov, lower = rep(0, d))  
              }
              else{
                theta_cov <- 0.1^2 * diag(d)/d
                new_proposals[j, ] <- rtmvnorm(1, theta_new, theta_cov, lower = rep(0, d)) 
              }
            }
          }
          
          New_proposals <- rbind(new_proposals, theta_old)
          new_likelihoods <- alpha * apply(New_proposals, 1, likelihood, data = data)
          new_priors <- apply(New_proposals, 1, prior)
          
          logwxy <- new_likelihoods + new_priors 
          
          # Step 3: Accept or reject Y based on the acceptance ratio
          ratio <- log(sum(exp(logwyx)))- log(sum(exp(logwxy)))
          
          # accept/reject step
          if (log(runif(1, 0, 1)) < ratio){
            return(list(theta = theta_new, accept = T))
          } else{
            return(list(theta = theta_old, accept = F))
          }
        }
        # Calls the MH function
        MH(particle$theta, data, likelihood, prior, alpha[[r]])
      }   
      
      particles[[r]] <- lapply(particles[[r-1]],  MCMCmove, theta_cov)
      
      ## Converting the particles[[r]] list to a matrix with K rows and d columns
      theta_values_r_1 <- lapply(particles[[r - 1]], function(particle) particle$theta)
      theta_matrix_r_1 <- do.call(rbind, theta_values_r_1)
      
      # compute the ESS
      log_incremental_w <- alphaDiff * u
      logw <- log_incremental_w + logw  # log un-normalized weights are updated by adding incremental
      # weight to the previous log un-normalized weight
      logmax <- max(logw) # find the maximum of the updated log un-normalized weights
      logZ[[r]] <- logZ[[r-1]] + logsum(log_incremental_w + logW)  # update the log normalizing constant 
      W[[r]] <- exp(logw-logmax)/sum(exp(logw-logmax))   # normalized weights
      logW <- log(W[[r]]) # log normalized weights
      
      ESS[[r]] <- rESS(logW) # Calculate the ESS 
      
      # re-sample if ESS below threshold
      if( ESS[[r]]<tuning_param$eps )
      {
        cat("Resample: ESS=", ESS[[r]], '\n')
        ancestors <- systematic_resample( W[[r]] )
        particles[[r]] <-  particles[[r]][ancestors] #particles are then updated by selecting the particles 
        # corresponding to the ancestors indices. This step effectively replaces particles with low weights 
        # with duplicates of particles with higher weights
        W[[r]] <- rep(1/K,K) # resets the normalized weights
        logW <- log(W[[r]]) # re-calculate the logarithm of the normalized weights
        w <- rep(1,K) # un-normalized weights are reset to 1
        logw <- rep(0,K) # logarithm of un-normalized weights
      }
      
      # Calculate the weighted covariance matrix
      cov_weighted_l <- cov.wt(theta_matrix_r_1, W[[r - 1]])
      weighted_cov <- cov_weighted_l[[1]]
    }
    
    # these particles are the primary output of the algorithm. They represent the distribution of interest.
    # logZ is used to calculate the marginal likelihood
    output <- list(particles = particles,
                   alpha = alpha,
                   ESS = ESS,
                   logZ = logZ,
                   W = W)
    return(output)
  }
  
  # Solving I(t) DE trajectory using ODE solver
  # Solving I(t) DE trajectory using ODE solver
  library(deSolve)
  library(mvtnorm)
  library(MASS)
  library(coda)
  library(MCMCpack)
  library(truncnorm)
  library(tmvtnorm)
  
  # initial (state) values for SIR model
  N <- 22e6
  i0 <- 6
  x0 <- c(S = N - i0, E = 0.4 * i0, I = 0.6 * i0, R = 0)
  
  # vector of time steps
  times <- seq(0, 149, by = 1)
  
  # vector of parameters used in SIR model
  params <- c(beta = 0.1447, delta = 1/3.2, gamma = 1/8.5, logphi = log(0.3824282))
  
  SEIR <- function(t, x, params) {
    with(as.list(c(params, x)), {
      dS <- -beta * S * I / N 
      dE <- beta * S * I / N - delta * E
      dI <- delta * E - gamma * I
      dR <- gamma * I 
      return(list(c(dS, dE, dI, dR)))
    })
  }
  
  likelihood_NB <- list()
  likelihood_NB$dFun <- function(data, fit, ODEpars, likeliParam, log = TRUE) {
    with(as.list(likeliParam), {
      phi = exp(logphi)
      temp = dnbinom(data$y, size = phi, mu = fit[,"I"], log = log)
      return(sum(temp[!is.nan(temp)]))
    })
  }
  likelihood_NB$getMean <- function(fit, pars) {
    with(as.list(pars), {
      mu = fit[,"I"]
      return(data.frame(t = fit[,"time"], mu = mu))
    })
  }
  
  # Single simulation function
  sim_data <- function(n, likelihood, ODEmodel, initial_states, times, pars, ...) {
    fit <- rk4(y = initial_states, times = times, func = ODEmodel, parms = pars, ...)
    means <- likelihood$getMean(fit, pars)
    sampled_indices <- seq(1, nrow(means), length.out = n)
    phi <- exp(pars["logphi"])
    
    list(
      y = rnbinom(n, size = phi, mu = means$mu[sampled_indices]),
      t = means$t[sampled_indices],
      n_DE = length(initial_states),
      DE_val_names = names(initial_states),
      initial_state = x0,
      partial_likelihood = FALSE,
      fit = fit # Return the full ODE solution
    )
  }
  
  # Simulate data
  simulated_data <- sim_data(n = 150, likelihood = likelihood_NB, ODEmodel = SEIR, initial_states = x0, times = times, pars = params)
  
  # Extract simulated values
  simulated_y <- simulated_data$y
  simulated_t <- simulated_data$t
  
  # Define the Likelihood
likelihood <- function(data, theta) {
  beta <- theta[1] # transmission rate
  delta <- theta[2] # incubation rate
  gamma <- theta[3] # recovery rate
  phi <- theta[4] # dispersion parameter of the negative binomial model
  
  params2 <- c(beta, delta, gamma)
  
  SEIR_model <- function(t, x, params2) {
    with(as.list(c(params2, x)), {
      dS <- -beta * S * I / N 
      dE <- beta * S * I / N - delta * E
      dI <- delta * E - gamma * I
      dR <- gamma * I 
      list(c(dS, dE, dI, dR))
    })
  }
  S0 <- N - i0
  E0 <- 0.4 * i0
  I0 <- 0.6 * i0
  x_0 <- c(S = S0, E = E0, I = I0, R = N - (E0 + I0 + S0))
  
  fit <- ode(y = x_0, times = times, func = SEIR_model, parms = params2)
  fit <- as.data.frame(fit)
  mu <- fit$I
  
  # Calculate the log-likelihood
  Loglikelihood <- sum(dnbinom(simulated_y, size = phi, mu = mu, log = TRUE))
  return(Loglikelihood)
}

prior <- function(theta){
  beta <- theta[1] 
  delta <- theta[2]
  gamma <- theta[3] 
  phi <- theta[4] 
  
  log_prior_beta <- dlnorm(beta, log(0.1447), 0.2, log = TRUE)  
  log_prior_delta <- dlnorm(delta, log(1/3.2), 0.15, log = TRUE)  
  log_prior_gamma <- dlnorm(gamma, log(1/8.5), 0.2, log = TRUE)
  log_prior_phi <- dgamma(phi, shape = 382.4282, rate =1000, log = TRUE)
  
  return(log_prior_beta + log_prior_delta + log_prior_gamma + log_prior_phi)
}

reference <- prior

J <- 5
K <- 100
d <- 4
init <- function(){c(rlnorm(1,log(0.1447), 0.2),
                     rlnorm(1,log(1/3.2), 0.15),
                     rlnorm(1,log(1/8.5), 0.2), 
                     rgamma(1, shape = 382.4282, rate =1000)
)}

b <- 0.95
tuning_param <- list(eps = 0.5, eta = 0.9)
  
  # Run ASMC
  asmc_output <- asmc(K, theta_cov, tuning_param, init, data, likelihood, prior, reference)
  R1 <- length(asmc_output$particles)
  
  theta_R1 <- t(sapply(asmc_output$particles[[R1]], `[[`, "theta"))
  
  W_R1 <- asmc_output$W[[R1]]
  
  # Post-run resampling
  resampled_indices <- sample(x = 1:K, size = K, replace = TRUE, prob = W_R1)
  resampled_particles <- theta_R1[resampled_indices,]
  
  # Estimation of the parameters
  mean_theta <- colMeans(resampled_particles)
  
  # Initialize a matrix to store the credible intervals
  credible_intervals <- matrix(nrow = 2, ncol = ncol(resampled_particles))
  rownames(credible_intervals) <- c("Lower", "Upper")
  colnames(credible_intervals) <- colnames(resampled_particles)
  
  # Calculate 95% credible intervals for each parameter
  for(i in 1:ncol(resampled_particles)) {
    credible_intervals[, i] <- quantile(resampled_particles[, i], probs = c(0.025, 0.975))
  }
  
  Parameters <- c("Beta", "Delta", "Gamma", "Phi")
  result <- data.frame(Parameters, mean_theta, t(credible_intervals))
  
  ASMC_AMTM_rsp <- unlist(resampled_particles)
  alpha_values <- unlist(asmc_output$alpha)
  ess_values <- unlist(asmc_output$ESS)
  W_values <-unlist(asmc_output$W)
  logZ_values <- unlist(asmc_output$logZ)
  
  ASMC_AMTM_ESS <- data.frame(ESS = ess_values)
  ASMC_AMTM_alpha <- data.frame(alpha_val = alpha_values)
  ASMC_AMTM_W <- data.frame(W = W_values)
  ASMC_AMTM_logZ <- data.frame(logZ = logZ_values)
  
  # Save the outputs for this seed
 # write.table(ASMC_AMTM_ESS, file = paste0("/scratch/theniw8/ASMC_AMTM/ASMC_AMTM_SEIR_", seed, "_ESS.csv"), sep = "," ,row.names = FALSE)
 # write.table(ASMC_AMTM_alpha, file = paste0("/scratch/theniw8/ASMC_AMTM/ASMC_AMTM_SEIR_", seed, "_alpha.csv"), sep = "," ,row.names = FALSE)
  write.table(ASMC_AMTM_rsp, file = paste0("/scratch/theniw8/ASMC_AMTM/ASMC_AMTM_SEIR_", seed, "_rsp.csv"), sep = ",")
  
  # Optionally, uncomment these lines if you need these outputs
  # write.table(ASMC_AMTM_W, file = paste0("/scratch/theniw8/ASMC_AMTM/ASMC_AMTM_SEIR_", seed, "_W.csv"), sep = "," ,row.names = FALSE)
  # write.table(ASMC_AMTM_logZ, file = paste0("/scratch/theniw8/ASMC_AMTM/ASMC_AMTM_SEIR_", seed, "_logZ.csv"), sep = "," ,row.names = FALSE)
}
