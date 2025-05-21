using Distributions, Random, LinearAlgebra

function mcmcmove(particle, data, likelihood, prior, alpha, r, d, b, weighted_cov)
    theta_old = Float64.(particle.theta)

    #@show size(weighted_cov)
    #@show weighted_cov

    # Robust Adaptive MH covariance logic
    if r <= 2 * d || rand() >= b
        theta_cov = (0.1^2 / d) * I(d)
    else
        Σ = (weighted_cov + weighted_cov') / 2         # force symmetry
        Σ += 1e-6 * I(d)                                # add jitter for PD
        theta_cov = (2.38^2 / d) * Σ
    end

    # Sample proposal from truncated MVN
    mvn = MvNormal(theta_old, theta_cov)
    lower = zeros(d)
    sampler = TruncatedMVNSampler(mvn, lower)
    theta_new = rand(sampler)
    
    # compute MH ratio
    num = likelihood(data, theta_new) * alpha + prior(theta_new)
    den = likelihood(data, theta_old) * alpha + prior(theta_old)
    ratio = num - den

    if log(rand()) < ratio
        return (theta = theta_new, accept = true)
    else
        return (theta = theta_old, accept = false)
    end
end
