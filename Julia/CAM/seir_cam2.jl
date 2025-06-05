# SEIR model
include("bisection.jl")
include("additionalfns.jl")
include("tmvn_full.jl")
include("weightedcov.jl")
include("asmc.jl")

using DifferentialEquations, Distributions, StatsBase, Random

Random.seed!(444)

# initial (state) values for the SEIR model
N = 22e6
i0 = 2.0
x0 = [N - i0, 0.4 * i0, 0.6 * i0, 0.0]

t = 0:1:136

params = (beta = 0.1447, delta = 1/3.2, gamma = 1/8.5, logphi = log(0.3415802))

# SEIR function
function seir!(du, u, p, t)
    S, E, I, R = u
    β, δ, γ = p.beta, p.delta, p.gamma

    du[1] = - β * u[1] * u[3] / N
    du[2] = β * (u[1] * u[3] / N) - δ * u[2]
    du[3] = δ * u[2] - γ * u[3]
    du[4] = γ * u[3]
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

# simulated data
simulated_data = sim_data(
    137,            # n
    likelihood_NB_dfun,  # likelihood function/module
    seir!,          # your ODE model function
    x0,             # initial_states
    t,          # times
    params          # parameters
)

# Extract simulated values
simulated_y = simulated_data.y
simulated_t = simulated_data.t

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

    log_prior_beta  = logpdf(LogNormal(log(0.1447), 0.25), beta)
    log_prior_delta = logpdf(LogNormal(log(1/3.2), 0.2), delta)
    log_prior_gamma = logpdf(LogNormal(log(1/8.5), 0.25), gamma)
    log_prior_phi   = logpdf(Gamma(382.4282, 1/1000), phi)
    return log_prior_beta + log_prior_delta + log_prior_gamma + log_prior_phi
end

function target_dis(data, theta)
    return likelihood(data, theta) + prior(theta)
end

# Reference Distribution
reference = prior

# Number of particles and dimension
K = 2500

d = 4
M = 6
Mx = 3
n = 1000

# Initial theta sampler (matching rlnorm and rgamma)
init() = [
    rand(LogNormal(log(0.1447), 0.25)),
    rand(LogNormal(log(1 / 3.2), 0.2)),
    rand(LogNormal(log(1 / 8.5), 0.25)),
    rand(Gamma(382.4282, 1 / 1000))
]

# Tuning parameters
b = 0.95
tuning_param = Dict(:eps => 0.5, :eta => 0.95)

data = simulated_y

# Call your ASMC function (assumes it's defined in asmc.jl)
asmc_output = asmc(K, d, b, tuning_param, init, data, likelihood, prior, reference);

include("cam_2.jl")
include("agd.jl")
include("local_cov.jl")

# Step 1: Initial theta and cov
theta_init = init()

# Step 2: Estimate local covariance (diagonal for simplicity)
cov_local = Diagonal(0.01 .* ones(d))

# Step 4: Run CAM
samples = cam(data, target_dis, theta_init, n, Mx, M, cov_local, asmc_output)

include("analysis_cam.jl")