# ASMC-post AGD

using Distributions, LinearAlgebra, StatsBase

function agd_all(asmc_output, I::Int)
    R_total = length(asmc_output.particles)
    z_list = []  

    for _ in 1:I
        # 1. Uniformly select an intermediate distribution (from 1 to R_total)
        r = rand(1:R_total)
        
        # 2. Extract particles and weights for the chosen iteration
        theta_r = Matrix(hcat([p.theta for p in asmc_output.particles[r]]...)')
        W_r = asmc_output.W[r]
        
        # 3. Sample a particle based on its normalized weight for that iteration
        idx = sample(1:size(theta_r, 1), Weights(W_r))
        push!(z_list, theta_r[idx, :])
    end

    return z_list  
end
