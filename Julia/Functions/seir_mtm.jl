using Pkg
Pkg.activate("/scratch/theniw8/")

# SEIR model
include("bisection.jl")
include("additionalfns.jl")
include("mtmmove.jl")
include("tmvn.jl")
include("asmc_mtm.jl")

using DifferentialEquations, Distributions, StatsBase, Random, CSV, DataFrames

# initial (state) values for the SEIR model
N = 22e6
i0 = 2.0
x0 = [N - i0, 0.4 * i0, 0.6 * i0, 0.0]

t = 0:1:136

params = (beta = 0.1447, delta = 1/3.2, gamma = 1/8.5, logphi = log(10))

# SEIR function
function seir!(du, u, p, t)
    S, E, I, R = max.(u, 0.0)
    β, δ, γ = p.beta, p.delta, p.gamma

    du[1] = - β * S * I / N
    du[2] = β * S * I / N - δ * E
    du[3] = δ * E - γ * I
    du[4] = γ * I
end

# Data simulation
function likelihood_NB_dfun(data, fit, ODEpars, params; log=true)
    phi = exp(params.logphi)
    mu = fit[:, "I"]  

    # Clamp probabilities to valid range (0, 1]
    p = clamp.(phi ./ (phi .+ mu), eps(), 1.0)
    dist = NegativeBinomial.(phi, p)

    probs = logpdf.(dist, data[:y]) 
    return sum(probs[.!isnan.(probs)])
end

function likelihood_NB_getMean(fit, params)
    mu = [u[3] for u in fit.u]  # I values
    t  = fit.t                  # time points
    return (t = t, mu = mu)
end

# single simulation function
function sim_data(n, likelihood, ODEmodel, initial_states, times, params)
    # Solve the ODE
    prob = ODEProblem(ODEmodel, initial_states, (times[1], times[end]), params)
    sol = solve(prob, Tsit5(), saveat=times)

    # Extract model means (e.g. I values over time)
    means = likelihood_NB_getMean(sol, params)  # returns (t = ..., mu = ...)
    t_vals = means.t
    mu_vals = means.mu

    # Sample n evenly spaced indices
    sampled_indices = round.(Int, range(1, length(mu_vals), length=n))

    # Compute dispersion parameter
    phi = exp(params.logphi)

    # Simulate observations
    y = [rand(NegativeBinomial(phi, clamp(phi / (phi + mu_vals[i]), eps(), 1.0))) for i in sampled_indices]

    return (
        y = y,
        t = t_vals[sampled_indices],
        n_DE = length(initial_states),
        DE_val_names = propertynames(initial_states),
        initial_state = initial_states,
        partial_likelihood = false,
        fit = sol  # the full ODE solution
    )
end

# Define the likelihood
function likelihood(data, theta)
    beta, delta, gamma, phi = theta

    S0 = N - i0
    E0 = 0.4 * i0
    I0 = 0.6 * i0
    R0 = N - (E0 + I0 + S0)

    x_0 = [S0, E0, I0, R0]

    # Define parameter tuple for ODE
    params2 = (beta = beta, delta = delta, gamma = gamma)

    # Solve the ODE
    prob = ODEProblem(seir!, x_0, (minimum(t), maximum(t)), params2)
    sol = solve(prob, Tsit5(), saveat=t)

    # Extract predicted I(t)
    mu = [u[3] for u in sol.u]  # I is the 3rd state

    # Compute NB log-likelihood
    p = clamp.(phi ./ (phi .+ mu), eps(), 1.0)
    dist = NegativeBinomial.(phi, p)
    Loglikelihood = sum(logpdf.(dist, data))

    return Loglikelihood 
end

# Prior Distribution
function prior(theta)
    beta, delta, gamma, phi = theta

    log_prior_beta  = logpdf(LogNormal(log(0.1447), 0.1), beta)
    log_prior_delta = logpdf(LogNormal(log(1/3.2), 0.07), delta)
    log_prior_gamma = logpdf(LogNormal(log(1/8.5), 0.05), gamma)
    log_prior_phi   = logpdf(Gamma(50, 0.2), phi)
    return log_prior_beta + log_prior_delta + log_prior_gamma + log_prior_phi
end

# Reference Distribution
reference = prior

# Initial theta sampler (matching rlnorm and rgamma)
init() = [
    rand(LogNormal(log(0.1447), 0.1)),
    rand(LogNormal(log(1 / 3.2), 0.07)),
    rand(LogNormal(log(1 / 8.5), 0.05)),
    rand(Gamma(50, 0.2))
]

# Number of particles and dimension
K = 55
d = 4
J = 10

# Tuning parameters
b = 0.95
tuning_param = Dict(:eps => 0.5, :eta => 0.95)

H = 50

Threads.@threads for i in 1:H
    println("Running dataset $i")
    Random.seed!(1000 + i)  # Reproducible results per dataset

    simulated_data = sim_data(137, likelihood_NB_dfun, seir!, x0, t, params)
    data = simulated_data.y

    asmc_output = asmc(K, d, b, tuning_param, init, data, likelihood, prior, reference, J)

    R1 = length(asmc_output.particles)
    theta_R1 = reduce(vcat, [reshape(p.theta, 1, :) for p in asmc_output.particles[R1]])

    W_R1 = asmc_output.W[R1]
    resampled_indices = sample(1:K, Weights(W_R1), K, replace=true)
    resampled_particles = theta_R1[resampled_indices, :]

    #println("resampled particle example:", resampled_particles)
    
    df_resampled = DataFrame(resampled_particles, :auto)
    rename!(df_resampled, ["Beta", "Delta", "Gamma", "Phi"])
    CSV.write("ASMC-AMTM-SCB-rsp-$(i).csv", df_resampled)

    # Save resampled particles (flattened)
    #ASMC_AMTM_rsp = vec(resampled_particles)  # Flattened like R's unlist()
    #CSV.write("ASMC-AMTM-SCB-rsp-$(i).csv", DataFrame(rsp = ASMC_AMTM_rsp))


    # Save alpha, ESS, logZ (all unlisted vectors)
    CSV.write("ASMC-AMTM-SCB-alpha-$(i).csv", DataFrame(alpha_val = asmc_output.alpha))
    CSV.write("ASMC-AMTM-SCB-ESS-$(i).csv", DataFrame(ESS = asmc_output.ESS))
    CSV.write("ASMC-AMTM-SCB-logZ-$(i).csv", DataFrame(logZ = asmc_output.logZ))
    CSV.write("ASMC-AMTM-SCB-W-$(i).csv", DataFrame(W = reduce(vcat, asmc_output.W)))

end
