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

# relative conditional effective sample size
# W - vector of weights
# u - vector of same length as W, represents some adjustment or transformation to the weights
# a - scalar value used to scale the transformation u
# phi - user specified threshold (eta)
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
} # Higher ESS are better

# systematic re-sampling algorithm # Commonly used in particle filters and other Monte Carlo methods
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
        j <- j + 1
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
  
  # initialize values
  # initialize the iterations
  r <- 1                  # SMC iteration # SMC process starts with the first iteration.
  #initialize the log-normalizing constant
  logZ[[r]] <- 0 # log-normalizing constant for the 1st iteration is set to 0.
  
  # initialize the annealing parameter 
  if(adaptive_alpha){
    alpha[[r]] <- 0    # tempering parameters
    alphaDiff <- 0     # difference b/w two successive tempering parameters
  } else{
    alpha <- as.list(tuning_param$alpha)
  }
  
  # initialize particles
  # Create a list with two elements for each particle (k). First element is initial values of theta and
  # second element is a logical value set to False). init() is a user defined function.
  particles[[r]] <- lapply(1:K, function(k) list(theta = init(), accept = F))
  
  # initialize the weighted covariance matrix
  weighted_cov <- diag(d)
  ## Initialize weights
  # normalized weights
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
      # if adaptive_alpha is True, alpha is updated using the bisection function based on the target "phi"
      # in tuning_param list
      alphaDiff <- bisection( 0, 1,  W[[r-1]], u, tuning_param$eta)
      alpha[[r]] <- alpha[[r-1]] + alphaDiff
    } else{
      alphaDiff <- alpha[[r]] - alpha[[r-1]]
    }
    
    cat("annealing parameter:",alpha[[r]],"\n") # prints the current alpha value for the ongoing iteration.
    # if alpha is set greater than 1, fix by setting to 1
    if( alpha[[r]]>1 ){
      alpha[[r]] <- 1
      alphaDiff <- 1-alpha[[r-1]] # alphaDiff is recalculated
    }
    
    # Updated MCMCmove Function
    MCMCmove <- function(particle, theta_cov){
      # Adaptive Metropolis-Hastings
      MH <- function(theta_old, data, likelihood, prior, alpha){
        
        theta_old <- as.numeric(theta_old)
        # Applying Adaptive Metropolis
        if(r <= 2 *d){
          theta_cov <- 0.1^2 * diag(d) / d
          theta_new <- rtmvnorm(1, theta_old, theta_cov, lower = rep(0, d))
        }
        else{
          if(runif(1) < b){
            theta_cov <- 2.38^2 * weighted_cov / d
            theta_new <- rtmvnorm(1, theta_old, theta_cov, lower = rep(0, d))
          }
          else{
            theta_cov <- 0.1^2 * diag(d) / d
            theta_new <- rtmvnorm(1, theta_old, theta_cov, lower = rep(0, d))
          }
        }
        
        # compute numerator and denominator of the MH ratio         
        num <- (likelihood(data, theta_new)) * alpha +  prior(theta_new)
        den <- (likelihood(data, theta_old)) * alpha + prior(theta_old)
        
        # compute acceptance ratio
        ratio <- num - den # min(1, exp(num - den))
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
      particles[[r]] <-  particles[[r]][ancestors] 
      W[[r]] <- rep(1/K,K) # resets the normalized weights
      logW <- log(W[[r]]) # re-calculate the logarithm of the normalized weights
      w <- rep(1,K) # un-normalized weights are reset to 1
      logw <- rep(0,K) # logarithm of un-normalized weights
    }
    
    # Calculate the weighted covariance matrix
    cov_weighted_l <- cov.wt(theta_matrix_r_1, W[[r - 1]])
    weighted_cov <- cov_weighted_l[[1]]
  }
  
  # logZ is used to calculate the marginal likelihood
  output <- list(particles = particles,
                 alpha = alpha,
                 ESS = ESS,
                 logZ = logZ,
                 W = W)
  return(output)
}


set.seed(444)

# Solving I(t) DE trajectory using ODE solver
library(deSolve)
library(mvtnorm)
library(MASS)
library(coda)
library(MCMCpack)
library(truncnorm)
library(tmvtnorm)


# initial (state) values for SIR model
N <- 1000
x0 <- c(S = N-14, E = 8,  I = 6, R = 0)

# vector of time steps
times <- 0:59

# vector of parameters used in SIR model
params <- c(beta = 0.2976, delta = 1/3.2, gamma = 1/8.5)

SEIR <- function(t, x, params) {
  with(as.list(c(params, x)), {
    dS <- -beta * S * I / N 
    dE <- beta * S * I / N - delta * E
    dI <- delta * E - gamma * I
    dR <- gamma * I 
    list(c(dS, dE, dI, dR))
  })
}

r1 <-  rk4(x0, times, SEIR, params)

I_values <- r1[, "I"]

mu_t <-  I_values

h <- 100
p1 <- h/(h + mu_t)

data <- rnbinom(60, size = h, mu = mu_t)

data_table <- data.frame(Time = times, I_Values = round(I_values),Y = data)

# Define the Likelihood
likelihood <- function(data,theta) {
  beta <- theta[1] # transmission rate
  delta <- theta[2] # incubation rate
  gamma <- theta[3] # recovery rate
  phi <- theta[4] # dispersion parameter of the negative binomial model
  S0 <- theta[5] # initial number of susceptible individuals
  E0 <- theta[6] # initial number if exposed individuals
  I0 <- theta[7] # initial number of infected individuals
  
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
  
  x_0 <- c(S = S0, E = E0, I = I0, R = N - (S0 + E0 + I0))
  out <-  rk4(x_0, times, SEIR_model, params2)
  
  I_values <- out[, "I"]
  mu_t <- I_values
  
  # Calculate the log-likelihood
  Loglikelihood <- sum(log(dnbinom(data, size = phi, mu = mu_t)))
  
  return(Loglikelihood)
}

prior <- function(theta){
  beta <- theta[1] # transmission rate
  delta <- theta[2] # incubation rate
  gamma <- theta[3] # recovery rate
  phi <- theta[4] # dispersion parameter of the negative binomial model
  S0 <- theta[5] # initial number of susceptible individuals
  E0 <- theta[6] # initial number if exposed individuals
  I0 <- theta[7] # initial number of infected individuals
  
  log_prior_beta <- dlnorm(beta, log(0.2976), 0.2, log = TRUE)  
  log_prior_delta <- dlnorm(beta, log(1/3.2), 0.2, log = TRUE)  
  log_prior_gamma <- dlnorm(gamma, log(1/8.5), 0.2, log = TRUE)
  log_prior_phi <- dgamma(phi, shape = 20000, rate = 200, log = TRUE)
  log_prior_S0 <- dnorm(S0, 986, 1, log = TRUE)
  log_prior_E0 <- dnorm(E0, 8, 1, log = TRUE)
  log_prior_I0 <- dnorm(I0, 6, 1, log = TRUE)
  
  return(log_prior_beta + log_prior_delta + log_prior_gamma + log_prior_phi + 
           log_prior_S0  + log_prior_E0 + log_prior_I0)
}

reference <- prior

K <-5000
d <- 7
init <- function(){c(rlnorm(1,log(0.2976), 0.2),
                     rlnorm(1,log(1/3.2), 0.2),
                     rlnorm(1,log(1/8.5), 0.2), 
                     rgamma(1, 20000,rate = 200),
                     rnorm(1, 986, 1),
                     rnorm(1, 8, 1),
                     rnorm(1, 6, 1))}

b <- 0.95
tuning_param <- list(eps = 0.5, eta = 0.99)
asmc_output <- asmc(K, theta_cov, tuning_param, init, data, likelihood, prior, reference)

R1 <- length(asmc_output$particles)

theta_1 <- t(sapply(asmc_output$particles[[1]], `[[`, "theta"))
theta_5 <- t(sapply(asmc_output$particles[[5]], `[[`, "theta"))
theta_10 <- t(sapply(asmc_output$particles[[10]], `[[`, "theta"))
theta_15 <- t(sapply(asmc_output$particles[[15]], `[[`, "theta"))
theta_R1 <- t(sapply(asmc_output$particles[[R1]], `[[`, "theta"))

W_R1 <- asmc_output$W[[R1]]

# post-run resampling
resampled_indices <- sample(x = 1:K, size = K, replace = TRUE, 
                            prob = W_R1)

resampled_particles <- theta_R1[resampled_indices,]

AASMC_rsp <- unlist(resampled_particles)
alpha_values <- unlist(asmc_output$alpha)
ess_values <- unlist(asmc_output$ESS)
W_values <-unlist(asmc_output$W)
logZ_values <- unlist(asmc_output$logZ)

AASMC_ESS <- data.frame(ESS = ess_values)
AASMC_alpha <- data.frame(alpha_val = alpha_values)
AASMC_W <- data.frame(W = W_values)
AASMC_logZ <- data.frame(logZ = logZ_values)

write.table(AASMC_ESS, file = "/scratch/theniw8/AASMC_Finale_1000_0.99_ESS.csv", sep = "," ,row.names = FALSE)
write.table(AASMC_alpha, file = "/scratch/theniw8/AASMC_Finale_1000_0.99_alpha.csv", sep = "," ,row.names = FALSE)
write.table(AASMC_W, file = "/scratch/theniw8/AASMC_Finale_1000_0.99_W.csv", sep = "," ,row.names = FALSE)
write.table(AASMC_logZ, file = "/scratch/theniw8/AASMC_Finale_1000_0.99_logZ.csv", sep = "," ,row.names = FALSE)
write.table(AASMC_rsp, file = "/scratch/theniw8/AASMC_Finale_1000_0.99_rsp.csv", sep = "," )
