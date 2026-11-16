using Distributions, Random, LinearAlgebra

function mtmmove(particle, data, likelihood, prior, alpha, r, d, b, weighted_cov, J)
    theta_old = Float64.(particle.theta)
    lower = zeros(d)
    # Propose J candidates from T(theta_old | .)
    theta_proposals = Matrix{Float64}(undef, J, d)

    for j in 1:J
        if r <= 2 * d || rand() >= b
            theta_cov = (0.1^2 / d) * I(d)
        else
            Σ = (weighted_cov + weighted_cov') / 2         # force symmetry
            Σ += 1e-6 * I(d)                                # add jitter for PD
            theta_cov = (2.38^2 / d) * Σ
        end
            # Sample proposal from truncated MVN
            sampler = TruncatedMVNSampler(MvNormal(theta_old, theta_cov), lower)
            theta_proposals[j, :] = rand(sampler)
    end
   
    # Compute log weights and sample one proposal
    likelihoods = [alpha[r] * likelihood(data, theta) for theta in eachrow(theta_proposals)]
    priors = [prior(theta) for theta in eachrow(theta_proposals)]
    logwyx = likelihoods .+ priors
    log_probs = logwyx .- logsum(logwyx)

    idx = sample(1:J, Weights(exp.(log_probs)))
    theta_new = theta_proposals[idx, :]
    
    # Generate backward proposals and compute reverse weights
    new_proposals = Matrix{Float64}(undef, J - 1, d)

    for j in 1: (J - 1)
        if r <= 2 * d || rand() >= b
            theta_cov = (0.1^2 / d) * I(d)
        else
            Σ = (weighted_cov + weighted_cov') / 2         # force symmetry
            Σ += 1e-6 * I(d)                                # add jitter for PD
            theta_cov = (2.38^2 / d) * Σ
        end
            # Sample proposal from truncated MVN
            sampler = TruncatedMVNSampler(MvNormal(theta_new, theta_cov), lower)
            new_proposals[j, :] = rand(sampler)
    end

    New_proposals = vcat(new_proposals, reshape(theta_old, 1, :))
    new_likelihoods = [alpha[r] * likelihood(data, theta) for theta in eachrow(New_proposals)]
    new_priors = [prior(theta) for theta in eachrow(New_proposals)]
    logwxy = new_likelihoods .+ new_priors

    # Accept/reject
    ratio = logsum(logwyx) - logsum(logwxy)
    if log(rand()) < ratio
        return (theta = theta_new, accept = true)
    else
        return (theta = theta_old, accept = false)
    end
end
