set.seed(444)
# Solving I(t) DE trajectory using ODE solver
library(deSolve)
library(mvtnorm)
library(MASS)
library(coda)
library(tictoc)
library(MCMCpack)
library(truncnorm)

# initial (state) values for SIR model
N <- 1000
x0 <- c(S = N-6, I = 6, R = 0)

# vector of time steps
times <- 0:59

# vector of parameters used in SIR model
params <- c(beta = 0.3, gamma = 0.1)

SIR <- function(t, x, params) {
  with(as.list(c(params, x)), {
    dS <- -beta * S * I / N 
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I 
    list(c(dS, dI, dR))
  })
}

r <-  rk4(x0, times, SIR, params)
plot(r)

I_values <- r[, "I"]

(I_max <- max(I_values))
mu_t <- 0.85 * I_values

h <- 100.5
p1 <- h/(h + mu_t)
y <- rnbinom(60, size = h, prob = p1)

data_table <- data.frame(Time = times, I_Values = round(I_values),Y = y)

plot(times,y, main = 'Observations', xlab = 'time')

(y_max <- max(y))

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
  
  # Calculate the log-likelihood
  loglikelihood <- sum(log(dnbinom(y, size = phi, prob = phi / (phi + mu_t))))
  
  # Calculate the prior distributions for each parameter
  log_prior_beta <- dlnorm(beta, log(0.3^2/sqrt(0.3^2 + 0.05^2)), sqrt(log(1 + (0.05^2/0.3^2))), log = TRUE)  
  log_prior_gamma <- dlnorm(gamma, log(0.1^2/sqrt(0.1^2 + 0.1^2)), sqrt(log(1 + (0.1^2/0.1^2))), log = TRUE)
  log_prior_phi <- dgamma(phi, shape = 20000, rate = 200, log = TRUE)
  log_prior_S0 <- dnorm(S0, 994, 1, log = TRUE)
  log_prior_I0 <- dnorm(I0, 6, 1, log = TRUE)
  
  # Return the product of the mixture and prior distributions
  
  return(loglikelihood + log_prior_beta + log_prior_gamma + log_prior_phi + 
           log_prior_S0  + log_prior_I0)
}

n_chains <- 4

sample_matrices <- list()

for (j in 1:n_chains){
  # Initialize starting point
  current_params <- c(runif(1,0.1,0.5), runif(1,0,0.4), runif(1,90,110), 
                      runif(1,990, 996), runif(1, 4, 9))
  
  n_iter <- 2500 #G
  n_iter0 <- 500 #G0
  p <- 100
  b <- 0.05
  b1 <- 1 - b
  d <- 5
  
  # Initialize empty matrix to store samples
  samples_matrix <- matrix(0,nrow = n_iter, ncol = d)
  Sigma <- diag(d)
  # Initialize counters for accepted and rejected proposals
  accepted_count <- 0
  
  # Run the Metropolis-Hastings algorithm
  for (i in 1:n_iter) {
    if (i < n_iter0){
      proposed_params <- rnorm(d, current_params, 0.1/d)
    }
    else{
      if (runif(1) < b){
        proposed_params <- rnorm(d, current_params, 0.1/ d)
      }
      else{
        proposed_params <- mvrnorm(1, current_params, (2.38^2 / d) * Sigma)
      }
    }
    
    # Step 4: Compute acceptance probability
    
    m1 <- target_distribution(proposed_params)
    m2 <- target_distribution(current_params)
    (r <- exp(m1-m2))
    
    a <- min(1, r)
    
    # Step 5: Accept or reject the proposal
    if (!is.na(a) && runif(1) < a) {
      current_params <- proposed_params  # Accept the proposal
      accepted_count <- accepted_count + 1
    }
    
    samples_matrix[i, ] <- current_params
    
    if (i %% p == 0 && i >= n_iter0) {
      Sigma <- cov(samples_matrix[1:i, ])
    }
  }
  
  # Calculate the mean acceptance rate
  mean_acceptance_rate <- accepted_count / n_iter
  cat("Mean Acceptance Rate:", mean_acceptance_rate, "\n")
  
  # Store the samples for this chain in the list
  sample_matrices[[j]] <- mcmc(samples_matrix)
}

# Convert the list of sample matrices into a mcmc.list object
mcmc_list <- mcmc.list(sample_matrices)

# Compute the Gelman-Rubin diagnostic (R-hat) for each parameter in each chain
(rhats <- gelman.diag(mcmc_list))

par(mfrow = c(3, 2)) 

# beta values
samples_beta <- lapply(sample_matrices, function(mat) mat[, 1])

# gamma values
samples_gamma <- lapply(sample_matrices, function(mat) mat[, 2])

# phi values
samples_phi <- lapply(sample_matrices, function(mat) mat[, 3])

# S0 values
samples_S0 <- lapply(sample_matrices, function(mat) mat[, 4])

# I0 values
samples_I0 <- lapply(sample_matrices, function(mat) mat[, 5])


# Create a vector of colors for the lines
line_colors <- c("red", "blue", "green", "orange")


# Plot for Beta
plot(1, type = "l", xlim = c(1, n_iter), ylim = range(unlist(samples_beta)), 
     xlab = "Iterations", ylab = "Beta", main = "Trace Plots for Beta")

for (i in 1:length(samples_beta)) {
  lines(samples_beta[[i]], col = line_colors[i])
}

legend("bottomright", legend = paste("C", 1:length(samples_beta)), 
       col = line_colors, lty = 1, ncol = length(samples_beta), cex = 0.8)


# Plot for Gamma
plot(1, type = "l", xlim = c(1, n_iter), ylim = range(unlist(samples_gamma)), 
     xlab = "Iterations", ylab = "Gamma", main = "Trace Plots for Gamma")

for (i in 1:length(samples_gamma)) {
  lines(samples_gamma[[i]], col = line_colors[i])
}

legend("topright", legend = paste("C", 1:length(samples_gamma)), 
       col = line_colors, lty = 1, ncol = length(samples_gamma), cex = 0.8)


# Plot for Phi
plot(1, type = "l", xlim = c(1, n_iter), ylim = range(unlist(samples_phi)), 
     xlab = "Iterations", ylab = "Phi", main = "Trace Plots for Phi")

for (i in 1:length(samples_phi)) {
  lines(samples_phi[[i]], col = line_colors[i])
}

legend("bottomright", legend = paste("C", 1:length(samples_phi)), 
       col = line_colors, lty = 1, ncol = length(samples_phi), cex = 0.8)

# Plot for S0
plot(1, type = "l", xlim = c(1, n_iter), ylim = range(unlist(samples_S0)), 
     xlab = "Iterations", ylab = "S0", main = "Trace Plots for S0")

for (i in 1:length(samples_S0)) {
  lines(samples_S0[[i]], col = line_colors[i])
}

legend("bottomright", legend = paste("C", 1:length(samples_S0)), 
       col = line_colors, lty = 1, ncol=length(samples_S0), cex = 0.8)


# Plot for I0
plot(1, type = "l", xlim = c(1, n_iter), ylim = range(unlist(samples_I0)), 
     xlab = "Iterations", ylab = "I0", main = "Trace Plots for I0")

for (i in 1:length(samples_I0)) {
  lines(samples_I0[[i]], col = line_colors[i])
}

legend("topright", legend = paste("C", 1:length(samples_I0)), 
       col = line_colors, lty = 1, ncol = length(samples_I0), cex = 0.8)


# Calculate Mean Parameter Values for Each Chain
(mean_beta <- lapply(sample_matrices, function(mat) mean(mat[, 1])))
(mean_gamma <- lapply(sample_matrices, function(mat) mean(mat[, 2])))
(mean_phi <- lapply(sample_matrices, function(mat) mean(mat[, 3])))
(mean_S0 <- lapply(sample_matrices, function(mat) mean(mat[, 4])))
(mean_I0 <- lapply(sample_matrices, function(mat) mean(mat[, 5])))

# Compute Credible Interval (95% CI) for Parameters
(quantiles_beta <- quantile(unlist(samples_beta), probs = c(0.025, 0.975)))
(quantiles_gamma <- quantile(unlist(samples_gamma), probs = c(0.025, 0.975)))
(quantiles_phi <- quantile(unlist(samples_phi), probs = c(0.025, 0.975)))
(quantiles_S0 <- quantile(unlist(samples_S0), probs = c(0.025, 0.975)))
(quantiles_I0 <- quantile(unlist(samples_I0), probs = c(0.025, 0.975)))
