using Pkg
Pkg.activate("/scratch/theniw8/")

# SEIR model
include("bisection.jl")
include("additionalfns.jl")
include("tmvn_full.jl")
include("weightedcov.jl")
include("asmc.jl")
include("cam_2.jl")
include("agd.jl")
include("local_cov.jl")

using DifferentialEquations, Distributions, StatsBase, Random, CSV, DataFrames 

# initial (state) values for the SEIR model
N = 22e6
i0 = 2.0
x0 = [N - i0, 0.4 * i0, 0.6 * i0, 0.0]

t = 0:1:136

params = (beta = 0.1447, delta = 1/3.2, gamma = 1/8.5, logphi = log(10)) #0.3415802 #0.3824282

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

    #println("mu values: ", mu_vals')

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
    sol = solve(prob, Rodas5(), saveat=t)

    # Extract predicted I(t)
    mu = [u[3] for u in sol.u]  # I is the 3rd state

    #if any(mu .<= eps())
    #println("⚠️ Small or zero mu detected for θ = ", theta)
    #println("min(mu) = ", minimum(mu), ", max(mu) = ", maximum(mu))
    #return -Inf
   # end


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

function target_dis(data, theta)
    return likelihood(data, theta) + prior(theta)
end

# Reference Distribution
reference = prior

# Initial theta sampler (matching rlnorm and rgamma)
init() = [
    rand(LogNormal(log(0.1447), 0.1)),
    rand(LogNormal(log(1 / 3.2), 0.07)),
    rand(LogNormal(log(1 / 8.5), 0.05)),
    rand(Gamma(50, 0.2)) #10.58, 0.2174
]

# Number of particles and dimension
K = 250

d = 4
M = 100
Mx = 50
n = 500
#burn_in = 100

# Tuning parameters
b = 0.95
tuning_param = Dict(:eps => 0.5, :eta => 0.95)

H = 50

Threads.@threads for i in 1 : H
    println("Running dataset $i")
    Random.seed!(1000 + i)

    simulated_data = sim_data(137, likelihood_NB_dfun, seir!, x0, t, params)
    # Extract simulated values
    data = simulated_data.y

    asmc_output = asmc(K, d, b, tuning_param, init, data, likelihood, prior, reference)
    theta_init = init()
    #cov_local = Diagonal(0.01 .* ones(d))
    cov_local = Diagonal([0.01, 0.0015, 0.002, 0.05])

    samples = cam(data, target_dis, theta_init, n, Mx, M, cov_local, asmc_output)

    # Convert samples (matrix) to DataFrame
    df = DataFrame(samples, :auto)  # :auto names the columns like x1, x2, x3, x4

    # Optional: Rename columns to parameter names
    rename!(df, ["Beta", "Delta", "Gamma", "Phi"])

    # Save to CSV
    CSV.write("cam_samples_dataset_$(i).csv", df)
end

