set.seed(444)

# Data Simulation

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
x0 <- c(S = N-6, I = 6, R = 0)

# vector of time steps
times <- 0:59

# vector of parameters used in SIR model
params <- c(beta = 0.2976, gamma = 1/8.5)

SIR <- function(t, x, params) {
  with(as.list(c(params, x)), {
    dS <- -beta * S * I / N 
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I 
    list(c(dS, dI, dR))
  })
}

r <-  rk4(x0, times, SIR, params)

I_values <- r[, "I"]

mu_t <- I_values

h <- 100
p1 <- h/(h + mu_t)

data <- rnbinom(60, size = h, mu = mu_t)

data_table <- data.frame(Time = times, I_Values = round(I_values),Y = data)

# Adaptive Metropolis Algorithm
AM_SIR <- function(n_iter, n_iter0, p, init_val, b, d){
  
  # Initialize empty matrix to store samples
  samples_matrix <- matrix(0,nrow = n_iter, ncol = d)
  samples_matrix[1, ] <- init_val
  Sigma <- diag(d)
  
  # Define the target distribution as a function of all unknown parameters
  target_distribution <- function(params) {
    beta <- params[1] # transmission rate
    gamma <- params[2] # recovery rate
    phi <- params[3] # dispersion parameter of the negative binomial model
    S0 <- params[4] # initial number of susceptible individuals
    I0 <- params[5] # initial number of infected individuals
    
    params2 <- c(beta, gamma)
    
    SIR_model <- function(t, x, params2) {
      with(as.list(c(params2, x)), {
        dS <- -beta * S * I / N 
        dI <- beta * S * I / N - gamma * I 
        dR <- gamma * I 
        list(c(dS, dI, dR))
      })
    }
    
    x_0 <- c(S = S0, I = I0, R = N - (S0 + I0))
    out <-  rk4(x_0, times, SIR_model, params2)
    
    I_values <- out[, "I"]
    mu_t <- I_values
    
    # Calculate the log-likelihood
    loglikelihood <- sum(log(dnbinom(data, size = phi, mu = mu_t)))
    
    # Calculate the prior distributions for each parameter
    log_prior_beta <- dlnorm(beta, log(0.2976), 0.2, log = TRUE)    
    log_prior_gamma <- dlnorm(gamma, log(1/8.5), 0.2, log = TRUE)
    log_prior_phi <- dgamma(phi, shape = 20000, rate = 200, log = TRUE)
    log_prior_S0 <- dnorm(S0, 994, 1, log = TRUE)
    log_prior_I0 <- dnorm(I0, 6, 1, log = TRUE)
      
    # Return the product of the mixture and prior distributions
    
    return(loglikelihood + log_prior_beta + log_prior_gamma + log_prior_phi + 
             log_prior_S0  + log_prior_I0)
  }
  
  for (i in 2:n_iter) {
    current_params <- samples_matrix[(i-1), ]
    if (i <= n_iter0){
      proposed_params <- rtmvnorm(1, current_params, (0.1^2 * diag(d)) / d, lower = rep(0, d))
    }
    else{
      if (runif(1) < b){
        proposed_params <- rtmvnorm(1, current_params, (2.38^2 * Sigma) / d, lower = rep(0, d))
        
      }
      else{
        proposed_params <- rtmvnorm(1, current_params, (0.1^2 * diag(d)) / d, lower = rep(0, d))
      }
    }
    
    m1 <- target_distribution(proposed_params)
    m2 <- target_distribution(current_params)
    (a <- m1-m2)
    
    if (a > log(runif(1, min = 0, max =1))) {
      current_params <- proposed_params  
    }
    
    samples_matrix[i, ] <- current_params
    
    if (i %% p == 0 && i >= n_iter0) {
      Sigma <- cov(samples_matrix[1:(i - 1), ])
    }
  } 
  return(samples_matrix)
}

# Results
init1 <- c(runif(1,0,0.5), runif(1,0,0.3), runif(1,95,105), runif(1,980, 1000), runif(1, 2, 15)); 
init2 <- c(runif(1,0,0.5), runif(1,0,0.3), runif(1,95,105), runif(1,980, 1000), runif(1, 2, 15));
init3 <- c(runif(1,0,0.5), runif(1,0,0.3), runif(1,95,105), runif(1,980, 1000), runif(1, 2, 15));
init4 <- c(runif(1,0,0.5), runif(1,0,0.3), runif(1,95,105), runif(1,980, 1000), runif(1, 2, 15));
init5 <- c(runif(1,0,0.5), runif(1,0,0.3), runif(1,95,105), runif(1,980, 1000), runif(1, 2, 15));

chain1 <- AM_SIR(62000, 5000, 100, init1, 0.95, 5);
chain2 <- AM_SIR(62000, 5000, 100, init2, 0.95, 5);
chain3 <- AM_SIR(62000, 5000, 100, init3, 0.95, 5);
chain4 <- AM_SIR(62000, 5000, 100, init4, 0.95, 5);
chain5 <- AM_SIR(62000, 5000, 100, init5, 0.95, 5);

# Discarding first 10000 draws as burn-ins
chain11 <- chain1[10001: 62000, ]
chain22 <- chain2[10001: 62000, ]
chain33 <- chain3[10001: 62000, ]
chain44 <- chain4[10001: 62000, ]
chain55 <- chain5[10001: 62000, ]

mcmc_chain1 <- mcmc(chain11)
mcmc_chain2 <- mcmc(chain22)
mcmc_chain3 <- mcmc(chain33)
mcmc_chain4 <- mcmc(chain44)
mcmc_chain5 <- mcmc(chain55)

write.table(chain1, file = "/scratch/theniw8/AMCMC_SIR_chain1.csv", sep = "," ,row.names = FALSE)
write.table(chain2, file = "/scratch/theniw8/AMCMC_SIR_chain2.csv", sep = "," ,row.names = FALSE)
write.table(chain3, file = "/scratch/theniw8/AMCMC_SIR_chain3.csv", sep = "," ,row.names = FALSE)
write.table(chain4, file = "/scratch/theniw8/AMCMC_SIR_chain4.csv", sep = "," ,row.names = FALSE)
write.table(chain5, file = "/scratch/theniw8/AMCMC_SIR_chain5.csv", sep = "," ,row.names = FALSE)
write.table(chain11, file = "/scratch/theniw8/AMCMC_SIR_chain11.csv", sep = "," ,row.names = FALSE)
write.table(chain22, file = "/scratch/theniw8/AMCMC_SIR_chain22.csv", sep = "," ,row.names = FALSE)
write.table(chain33, file = "/scratch/theniw8/AMCMC_SIR_chain33.csv", sep = "," ,row.names = FALSE)
write.table(chain44, file = "/scratch/theniw8/AMCMC_SIR_chain44.csv", sep = "," ,row.names = FALSE)
write.table(chain55, file = "/scratch/theniw8/AMCMC_SIR_chain55.csv", sep = "," ,row.names = FALSE)