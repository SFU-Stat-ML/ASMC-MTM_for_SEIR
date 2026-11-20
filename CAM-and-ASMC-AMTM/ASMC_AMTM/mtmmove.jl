using Distributions, Random, LinearAlgebra, RCall

function mtmmove(particle, data, likelihood, prior, alpha, r, d, b, weighted_cov, J)
    theta_old = Float64.(particle.theta)
    lower = zeros(d)
    cov_base = (0.1^2 / d) * I(d)
    Σ = (weighted_cov + weighted_cov') / 2
    Σ += 1e-6 * I(d)
    cov_adaptive = (2.38^2 / d) * Σ

    # Propose J candidates from T(theta_old | .)
    theta_proposals = Matrix{Float64}(undef, J, d)
    
    for j in 1:J
        if r <= 2 * d || rand() >= b
            theta_cov = cov_base
        else
            theta_cov = cov_adaptive
        end
            # Sample proposal from truncated MVN
            theta_proposals[j, :] = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $theta_old, sigma = $theta_cov, lower = $lower))""")
    end
   
    # computing T(y_j, x_t)
    logTyx = zeros(J)
    for j in 1:J
        if r <= 2 * d
            logTyx[j] = rcopy(Float64, R"as.numeric(tmvtnorm::dtmvnorm($theta_old, mean = $(theta_proposals[j, :]), sigma = $cov_base, lower = $lower, log = TRUE))")
        else
            logTyx_adaptive = rcopy(Float64, R"as.numeric(tmvtnorm::dtmvnorm($theta_old, mean = $(theta_proposals[j, :]), sigma = $cov_adaptive, lower = $lower, log = TRUE))")
            logTyx_base = rcopy(Float64, R"as.numeric(tmvtnorm::dtmvnorm($theta_old, mean = $(theta_proposals[j, :]), sigma = $cov_base, lower = $lower, log = TRUE))")
    
            logTyx[j] = log_mixture_density(logTyx_adaptive, logTyx_base, b)
        end
    end

    # Compute log weights and sample one proposal
    likelihoods = [alpha[r] * likelihood(data, theta) for theta in eachrow(theta_proposals)]
    priors = [prior(theta) for theta in eachrow(theta_proposals)]
    logwyx = likelihoods .+ priors .+ logTyx 

    # If every candidate has zero weight, stay at the current particle
    if any(isnan, logwyx) || any(==(Inf), logwyx)
        error("Invalid forward log weight: check the likelihood, prior, and proposal density")
    end

    if all(==(-Inf), logwyx)
        return (theta = theta_old, accept = false)
    end

    log_probs = logwyx .- logsum(logwyx)

    idx = sample(1:J, Weights(exp.(log_probs)))
    theta_new = theta_proposals[idx, :]
    
    # Generate backward proposals and compute reverse weights
    new_proposals = Matrix{Float64}(undef, J - 1, d)

    for j in 1: (J - 1)
        if r <= 2 * d || rand() >= b
            theta_cov = cov_base
        else
            theta_cov = cov_adaptive
        end
            
            new_proposals[j, :] = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $theta_new, sigma = $theta_cov, lower = $lower))""")
    end

    New_proposals = vcat(new_proposals, reshape(theta_old, 1, :))

    # computing T(x*_j, y)
    logTxy = zeros(J)
    for j in 1:J
        if r <= 2 * d
            logTxy[j] = rcopy(Float64, R"as.numeric(tmvtnorm::dtmvnorm($theta_new, mean = $(New_proposals[j, :]), sigma = $cov_base, lower = $lower, log = TRUE))")
        else
            logTxy_adaptive = rcopy(Float64, R"as.numeric(tmvtnorm::dtmvnorm($theta_new, mean = $(New_proposals[j, :]), sigma = $cov_adaptive, lower = $lower, log = TRUE))")
            logTxy_base = rcopy(Float64, R"as.numeric(tmvtnorm::dtmvnorm($theta_new, mean = $(New_proposals[j, :]), sigma = $cov_base, lower = $lower, log = TRUE))")
    
            logTxy[j] = log_mixture_density(logTxy_adaptive, logTxy_base, b)
        end
    end

    new_likelihoods = [alpha[r] * likelihood(data, theta) for theta in eachrow(New_proposals)]
    new_priors = [prior(theta) for theta in eachrow(New_proposals)]
    logwxy = new_likelihoods .+ new_priors .+ logTxy

    # Accept/reject 
    ratio = logsum(logwyx) - logsum(logwxy)
    if log(rand()) < ratio
        return (theta = theta_new, accept = true)
    else
        return (theta = theta_old, accept = false)
    end
end
