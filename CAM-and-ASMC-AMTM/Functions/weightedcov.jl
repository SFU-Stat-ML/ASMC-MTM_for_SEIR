using LinearAlgebra, Statistics, StatsBase

function weighted_covariance_matrix(X::Matrix{Float64}, weights::Vector{Float64})
    weights = weights / sum(weights)                   # ensure weights sum to 1
    μ = vec(sum(X .* weights, dims=1))                 # weighted mean
    X_centered = X .- μ'                               # center the particles
    return X_centered' * Diagonal(weights) * X_centered
end
