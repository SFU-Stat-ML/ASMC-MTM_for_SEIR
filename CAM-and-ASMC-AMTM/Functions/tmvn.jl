# Truncated Multivariate Normal Distribution

using Distributions
using Random
using LinearAlgebra

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
