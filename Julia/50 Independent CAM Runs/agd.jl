# ASMC-post AGD

using Distributions, LinearAlgebra, StatsBase

function agd(asmc_output, I::Int)
    R1 = length(asmc_output.particles)
    theta_R1 = hcat([p.theta for p in asmc_output.particles[R1]]...)'
    W_R1 = asmc_output.W[R1]

    z_list = []  

    for _ in 1:I
        idx = sample(1:size(theta_R1, 1), Weights(W_R1))
        push!(z_list, theta_R1[idx, :])
    end

    return z_list  
end
