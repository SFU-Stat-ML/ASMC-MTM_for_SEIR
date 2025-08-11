using Distributions, Random, LinearAlgebra

function mcmcmove(particle, data, likelihood, prior, alpha)
    theta_old = Float64.(particle.theta)
    d = length(theta_old)

    # Fixed proposal covariance
    theta_cov = Diagonal([0.03, 0.005, 0.015, 0.1])

    # Sample proposal
    mvn = MvNormal(theta_old, theta_cov)
    lower = zeros(d)
    sampler = TruncatedMVNSampler(mvn, lower)
    theta_new = rand(sampler)

    # MH ratio (log-domain)
    num = likelihood(data, theta_new) * alpha + prior(theta_new)
    den = likelihood(data, theta_old) * alpha + prior(theta_old)
    ratio = num - den

    if log(rand()) < ratio
        return (theta = theta_new, accept = true)
    else
        return (theta = theta_old, accept = false)
    end
end
