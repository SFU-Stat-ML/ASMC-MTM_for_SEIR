using LinearAlgebra, Statistics

"""
    local_cov(sample_matrix::Matrix{Float64}, d::Int, b::Float64, r::Int)

Given past MCMC samples (each row is a sample), compute the adaptive
Metropolis local proposal covariance matrix.
"""
function local_cov(sample_matrix::Matrix{Float64}, d::Int, b::Float64, r::Int)
    if r <= 2d || rand() >= b
        return (0.1^2 / d) * I(d)
    else
        # Compute empirical covariance from previous samples (1 to r-1)
        past_samples = sample_matrix[1:r-1, :]  # (r-1 × d)
        Σ = cov(past_samples, dims=1)           # Sample covariance along columns
        adaptive_cov = (2.38^2 / d) * Symmetric(Σ)
        return adaptive_cov
    end
end
