using Distributions, Random, LinearAlgebra, RCall

function mcmcmove(particle, data, likelihood, prior, alpha, r, d, b, weighted_cov)
    theta_old = Float64.(particle.theta)
    lower = zeros(d)
    cov_base = (0.1^2 / d) * I(d)
    Σ = (weighted_cov + weighted_cov') / 2
    Σ += 1e-6 * I(d)
    cov_adaptive = (2.38^2 / d) * Σ

    # Robust Adaptive MH covariance logic
    if r <= 2 * d || rand() >= b
        theta_cov = cov_base
    else
        theta_cov = cov_adaptive
    end

    sigma = Matrix(theta_cov)
    # Sample proposal from truncated MVN
    theta_new = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $theta_old, sigma = $sigma, lower = $lower))""")

    # Forward density: proposed state given current state.
    log_q_forward = rcopy(Float64, R"""as.numeric(tmvtnorm::dtmvnorm($theta_new, mean = $theta_old, sigma = $sigma, lower = $lower, log = TRUE))""")

    # Reverse density: current state given proposed state.
    log_q_reverse = rcopy(Float64, R"""as.numeric(tmvtnorm::dtmvnorm($theta_old, mean = $theta_new, sigma = $sigma, lower = $lower, log = TRUE))""")

    # compute MH ratio
    log_target_new = likelihood(data, theta_new) * alpha + prior(theta_new)
    log_target_old = likelihood(data, theta_old) * alpha + prior(theta_old)
    ratio = (log_target_new - log_target_old) + (log_q_reverse - log_q_forward)

    if log(rand()) < ratio
        return (theta = theta_new, accept = true)
    else
        return (theta = theta_old, accept = false)
    end
end
