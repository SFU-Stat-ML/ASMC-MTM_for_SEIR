getwd()
setwd("C:\\SFU\\Research\\MTM")
set.seed(444)

# Solving I(t) DE trajectory using ODE solver
library(deSolve)
library(mvtnorm)
library(MASS)
library(coda)
library(tictoc)
library(MCMCpack)
library(truncnorm)
library(numDeriv)

# initial (state) values for SIR model
N <- 1000
z0 <- c(S = N-6, I = 6, R = 0)

# vector of time steps
times <- 0:59

# vector of parameters used in SIR model
params <- c(beta = 0.3, gamma = 0.15)

SIR <- function(t, z, params) {
  with(as.list(c(params, z)), {
    dS <- -beta * S * I / N 
    dI <- beta * S * I / N - gamma * I
    dR <- gamma * I 
    list(c(dS, dI, dR))
  })
}

r <-  rk4(z0, times, SIR, params)
plot(r)

I_values <- round(r[, "I"])

(I_max <- max(I_values))
mu_t <- 0.85 * I_values

h <- 100.5
p1 <- h/(h + mu_t)
data <- rnbinom(60, size = h, prob = p1)

data_table <- data.frame(Time = times, I_Values = I_values, Y = data)

plot(times,data, main = 'Observations', xlab = 'time')

(data_max <- max(data))

# Define the target distribution as a function of all unknown parameters
td <- function(x) {
  beta <- x[1] # transmission rate
  gamma <- x[2] # recovery rate
  phi <- x[3] # dispersion parameter of the negative binomial model
  S0 <- x[4] # initial number of susceptible individuals
  I0 <- x[5] # initial number of infected individuals
  
  x <- c(beta, gamma, phi, S0, I0)
  
  params2 <- c(beta, gamma)
  
  SIR_model <- function(t, z, params2) {
    with(as.list(c(params2, z)), {
      dS <- -beta * S * I / N 
      dI <- beta * S * I / N - gamma * I
      dR <- gamma * I 
      list(c(dS, dI, dR))
    })
  }
  
  z_0 <- c(S = S0, I = I0, R = N - (S0 + I0))
  out <-  rk4(z_0, times, SIR_model, params2)
  
  I_values <- round(out[, "I"])
  
  # Calculate the log-likelihood
  loglikelihood <- sum(log(dnbinom(data, size = phi, prob = phi / (phi + mu_t))))
  
  # Calculate the prior distributions for each parameter
  log_prior_beta <- dlnorm(beta, log(0.3^2/sqrt(0.3^2 + 0.05^2)), sqrt(log(1 + (0.05^2/0.3^2))), log = TRUE)  
  log_prior_gamma <- dlnorm(gamma, log(0.12^2/sqrt(0.12^2 + 0.025^2)), sqrt(log(1 + (0.025^2/0.12^2))), log = TRUE)
  log_prior_phi <- dgamma(phi, shape = 20000, rate = 200, log = TRUE)
  log_prior_S0 <- dnorm(S0, 994, 1, log = TRUE)
  log_prior_I0 <- dnorm(I0, 6, 1, log = TRUE)
  
  # Return the product of the mixture and prior distributions
  
  return(loglikelihood + log_prior_beta + log_prior_gamma + log_prior_phi  + 
           log_prior_S0  + log_prior_I0)
}

n_chains <- 4

sample_matrices <- list()
for(l in 1:n_chains){
  # Initialize starting point
  xt <- c(runif(1,0.1,0.5), runif(1,0,0.4), runif(1,90,110), 
                      runif(1,990, 996), runif(1, 4, 9))
  
  n_iter <- 2500 #G
  n_iter0 <- 500 #G0
  pp <- 100
  b <- 0.05
  b1 <- 1 - b
  d <- 5
  k <- 5
  # Initialize empty matrix to store samples
  samples_matrix <- matrix(0,nrow = n_iter, ncol = d)
  Sigma <- diag(d)
  v_i <- c()
  # Initialize counters for accepted and rejected proposals
  accepted_count <- 0
  
# Run the Metropolis-Hastings algorithm
for (i in 1:n_iter) {
  
  # Propose sets of new values
  y <- matrix(NA, k, d)
 # pd <- c()
  for (j in 1:k){
    if (i < n_iter0){
      v0 <- abs(xt)
      params_sd <- 0.1 * v0
      y[j,] <- rnorm(d, xt, params_sd/d)
    }
    else{
      v_i <- colMeans(samples_matrix[1:i, ])
      if (runif(1) < b){
        y[j,] <- rnorm(d, xt, 0.1 * abs(v_i)/ d)
      }
      else{
        y[j,] <- rmvnorm(1, xt, (2.38^2 / d) * Sigma)
      }
    }
  } 
  
  # Compute the weights w(yj, x_t) for each trial proposal yj
  logp <- c()
  for (j in 1:k){
    logp[j] <- td(y[j,])
  }
  
  p <- exp(logp)
  
  # Step 2: Select Y among the trial set with probability proportional to weights
  selected_index <- sample(1:k, 1, prob = p)  # Select index based on weights
  
  Y <- y[selected_index,]  # Select Y based on the selected index
  
  # Draw x1*, ..., x(k-1)* from the distribution T(Y, .)
  X_star <- matrix(NA, k-1, d)
  for (j in 1:k - 1){
    if (i < n_iter0){
      v0 <- abs(Y)
      params_sd <- 0.1 * v0
      X_star[j,] <- rnorm(d, Y, params_sd/d)
    }
    else{
      v_i <- colMeans(samples_matrix[1:i, ])
      if (runif(1) < b){
        X_star[j,] <- rnorm(d, Y, 0.1 * v_i/ d)
      }
      else{
        X_star[j,] <- rmvnorm(1, Y, (2.38^2 / d) * Sigma)
      }
    }
  } 
  
  x_star <- rbind(X_star, xt)
  
  # Compute the weights w(x_starj,y) for each x_starj
  logpx <- c()
  for (j in 1:k){
    logpx[j] <- td(x_star[j,]) 
  }
  
  px <- exp(logpx)
  
  # calculating the acceptance probability
  a <- sum(p) / sum(px)
  rg <- min(1, a)
  
  
  # Decide whether to accept or reject the proposal
  if (!is.na(rg) && runif(1) < rg) {
    xt <- Y
    accepted_count <- accepted_count + 1
  }
  
  # Store the current parameters
  samples_matrix[i, ] <- xt
  
  if (i %% pp == 0 && i >= n_iter0) {
    Sigma <- cov(samples_matrix[1:i, ])
  }
} 
  # Calculate the mean acceptance rate
  mean_acceptance_rate <- accepted_count / n_iter
  cat("Mean Acceptance Rate:", mean_acceptance_rate, "\n")
  
  # Store the samples for this chain in the list
  sample_matrices[[l]] <- mcmc(samples_matrix)
}

# Convert the list of sample matrices into a mcmc.list object
mcmc_list <- mcmc.list(sample_matrices)

# Compute the Gelman-Rubin diagnostic (R-hat) for each parameter in each chain
(rhats <- gelman.diag(mcmc_list))

par(mfrow = c(1, 1)) 

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
