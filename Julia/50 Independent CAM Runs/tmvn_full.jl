# tmvn.jl - Truncated Multivariate Normal Support

using Distributions, Random, LinearAlgebra, Statistics, Base

# --- Sampler Definition ---

struct TruncatedMVNSampler <: Sampleable{Multivariate, Continuous}
    mvn::MvNormal
    lower::Vector{Float64}
end

Base.length(s::TruncatedMVNSampler) = length(s.lower)

function Distributions._rand!(rng::AbstractRNG, s::TruncatedMVNSampler, x::AbstractVector{T}) where T <: Real
    while true
        rand!(rng, s.mvn, x)
        if all(x .>= s.lower)
            return x
        end
    end
end

# --- PDF Estimation for Truncated MVN ---

"""
    tmvn_pdf(x::Vector{Float64}, μ::Vector{Float64}, Σ, lower::Vector{Float64}; N=10_000)

Estimate the PDF of a lower-truncated multivariate normal at point `x`.
Returns 0.0 if `x` is outside the truncation region.
Supports `Σ` as Matrix or PDMat.
"""
function tmvn_pdf(x::Vector{Float64}, μ::Vector{Float64}, Σ, lower::Vector{Float64}; N::Int=10_000)
    if any(x .< lower)
        return 0.0
    end

    Σmat = Matrix(Σ)  # ensure PDMat is converted
    mvn = MvNormal(μ, Σmat)
    numerator = pdf(mvn, x)

    # Monte Carlo estimate of the normalization constant
    count = 0
    buffer = zeros(length(μ))
    for _ in 1:N
        rand!(mvn, buffer)
        if all(buffer .>= lower)
            count += 1
        end
    end

    normalization_constant = count / N
    return numerator / normalization_constant
end

"""
    tmvn_pdf(x::Vector{Float64}, sampler::TruncatedMVNSampler; N=10_000)

Estimate the PDF using a `TruncatedMVNSampler` object.
"""
function tmvn_pdf(x::Vector{Float64}, sampler::TruncatedMVNSampler; N::Int=10_000)
    tmvn_pdf(x, mean(sampler.mvn), cov(sampler.mvn), sampler.lower; N=N)
end
