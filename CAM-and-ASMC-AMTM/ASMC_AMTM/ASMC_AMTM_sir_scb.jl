# SIR model
include("asmc_mtm.jl")

using DifferentialEquations, Distributions, StatsBase, Random, CSV, DataFrames 

# initial (state) values for the SIR model
N = 22e6
i0 = 2.0
x0 = [N - i0, i0, 0.0]

t = 0:1:136

params = (beta = 0.1447, gamma = 1/8.5, logphi = log(1))

# SIR function
function sir!(du, u, p, t)
    S, I, R = max.(u, 0.0)
    β, γ = p.beta, p.gamma

    du[1] = - β * S * I / N
    du[2] = β * S * I / N - γ * I
    du[3] = γ * I 
end

# Data simulation
function likelihood_NB_dfun(data, fit, ODEpars, params; log=true)
    phi = exp(params.logphi)
    mu = fit[:, "I"]  

    p = clamp.(phi ./ (phi .+ mu), eps(), 1.0)
    dist = NegativeBinomial.(phi, p)

    probs = logpdf.(dist, data[:y]) 
    return sum(probs[.!isnan.(probs)])
end

function likelihood_NB_getMean(fit, params)
    mu = [u[2] for u in fit.u]  # I values
    t  = fit.t                  # time points
    return (t = t, mu = mu)
end

# single simulation function
function sim_data(n, likelihood, ODEmodel, initial_states, times, params)
    # Solve the ODE
    prob = ODEProblem(ODEmodel, initial_states, (times[1], times[end]), params)
    sol = solve(prob, Tsit5(), saveat=times)

    # Extract model means 
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
    beta, gamma, phi = theta

    S0 = N - i0
    I0 = i0
    R0 = N - (I0 + S0)

    x_0 = [S0, I0, R0]

    # Define parameter tuple for ODE
    params2 = (beta = beta, gamma = gamma)

    # Solve the ODE
    prob = ODEProblem(sir!, x_0, (minimum(t), maximum(t)), params2)
    sol = solve(prob, Tsit5(), saveat=t)

    # Extract predicted I(t)
    mu = [u[2] for u in sol.u]  

    # Compute NB log-likelihood
    p = clamp.(phi ./ (phi .+ mu), eps(), 1.0)
    dist = NegativeBinomial.(phi, p)
    Loglikelihood = sum(logpdf.(dist, data))

    return Loglikelihood 
end

# Prior Distribution
function prior(theta)
    beta, gamma, phi = theta

    log_prior_beta  = logpdf(LogNormal(log(0.1447), 0.35), beta)
    log_prior_gamma = logpdf(LogNormal(log(1/8.5), 0.15), gamma)
    log_prior_phi   = logpdf(Gamma(2.0, 1.0), phi)
    return log_prior_beta + log_prior_gamma + log_prior_phi
end

function target_dis(data, theta)
    return likelihood(data, theta) + prior(theta)
end

# Reference Distribution
reference = prior

# Initial theta sampler (matching rlnorm and rgamma)
init() = [
    rand(LogNormal(log(0.1447), 0.35)),
    rand(LogNormal(log(1 / 8.5), 0.15)),
    rand(Gamma(2.0, 1.0)) 
]


# Number of particles and dimension
K = 639
d = 3
J = 10

# Tuning parameters
b = 0.95
tuning_param = Dict(:eps => 0.5, :eta => 0.9)

i = parse(Int, ARGS[1]) 

data_path = joinpath(@__DIR__, "simulated_sir_data.csv")
data = CSV.read(data_path, DataFrame; types=Dict(:y => Int)).y

@assert length(data) == length(t) "Expected one count per time point (137 rows)."

seed = 1000 + i
Random.seed!(seed)

r_seed = 2000 + i
R"set.seed($r_seed)"

println("Running ASMC_AMTM chain $i with Julia seed $seed and R seed $r_seed")

asmc_output = asmc(K, d, b, tuning_param, init, data, likelihood, prior, reference, J)
R1 = length(asmc_output.particles)
theta_R1 = reduce(vcat, [reshape(p.theta, 1, :) for p in asmc_output.particles[R1]])

W_R1 = asmc_output.W[R1]

resampled_indices = sample(1:K, Weights(W_R1), K, replace=true)
resampled_particles = theta_R1[resampled_indices, :]

df_resampled = DataFrame(resampled_particles, :auto)
rename!(df_resampled, ["Beta", "Gamma", "Phi"])
CSV.write("ASMC_AMTM_SIR_SCB_samples$(i).csv", df_resampled)

