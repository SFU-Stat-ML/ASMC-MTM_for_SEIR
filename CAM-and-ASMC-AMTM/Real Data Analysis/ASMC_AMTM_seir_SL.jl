# SEIR model
include("asmc_mtm.jl")

using DifferentialEquations, Distributions, StatsBase, Random, CSV, DataFrames, LinearAlgebra, RCall

# initial (state) values for the SEIR model
N = 22e6
i0 = 2.0
x0 = [N - i0, 0.4 * i0, 0.6 * i0, 0.0]

# Load SL data
df = CSV.read("SL_COVID_data.csv", DataFrame)
data = df.New_Confirmed
t = 0:(length(data)-1)

# SEIR function
function seir!(du, u, p, t)
    S, E, I, R = max.(u, 0.0)
    β, δ, γ = p.beta, p.delta, p.gamma

    du[1] = - β * S * I / N
    du[2] = β * S * I / N - δ * E
    du[3] = δ * E - γ * I
    du[4] = γ * I
end

# Define the likelihood
function likelihood(data, theta)
    beta, delta, gamma, phi = theta
    S0 = N - i0
    E0 = 0.4 * i0
    I0 = 0.6 * i0
    R0 = N - (E0 + I0 + S0)
    x_0 = [S0, E0, I0, R0]
    params2 = (beta = beta, delta = delta, gamma = gamma)
    prob = ODEProblem(seir!, x_0, (minimum(t), maximum(t)), params2)
    sol = solve(prob, Tsit5(), saveat=t)
    mu = [u[3] for u in sol.u]  
    p = clamp.(phi ./ (phi .+ mu), eps(), 1.0)
    dist = NegativeBinomial.(phi, p)
    Loglikelihood = sum(logpdf.(dist, data))
    return Loglikelihood 
end

# Prior Distribution
function prior(theta)
    beta, delta, gamma, phi = theta

    log_prior_beta  = logpdf(LogNormal(log(0.1447), 0.35), beta)
    log_prior_delta = logpdf(LogNormal(log(1/3.2), 0.25), delta)
    log_prior_gamma = logpdf(LogNormal(log(1/8.5), 0.15), gamma)
    log_prior_phi   = logpdf(Gamma(4.0, 0.5), phi)
    return log_prior_beta + log_prior_delta + log_prior_gamma + log_prior_phi
end

function target_dis(data, theta)
    return likelihood(data, theta) + prior(theta)
end

# Reference Distribution
reference = prior

# Initial theta sampler 
init() = [
    rand(LogNormal(log(0.1447), 0.35)),
    rand(LogNormal(log(1 / 3.2), 0.25)),
    rand(LogNormal(log(1 / 8.5), 0.15)),
    rand(Gamma(4.0, 0.5)) 
]

# Number of particles and dimension
K = 854
d = 4
J = 25

# Tuning parameters
b = 0.95

tuning_param = Dict(:eps => 0.5, :eta => 0.95)

i = parse(Int, ARGS[1])

seed = 1000 + i
Random.seed!(seed)

r_seed = 2000 + i
R"set.seed($r_seed)"

println("Running CAM chain $i with Julia seed $seed and R seed $r_seed")

asmc_output = asmc(K, d, b, tuning_param, init, data, likelihood, prior, reference, J)
R1 = length(asmc_output.particles)
theta_R1 = reduce(vcat, [reshape(p.theta, 1, :) for p in asmc_output.particles[R1]])
W_R1 = asmc_output.W[R1]

resampled_indices = sample(1:K, Weights(W_R1), K, replace=true)
resampled_particles = theta_R1[resampled_indices, :]

df_resampled = DataFrame(resampled_particles, :auto)
rename!(df_resampled, ["Beta", "Delta", "Gamma", "Phi"])

CSV.write("ASMC-AMTM-alpha-SL$(i).csv", DataFrame(alpha_val = asmc_output.alpha))
CSV.write("ASMC-AMTM-ESS-SL$(i).csv", DataFrame(ESS = asmc_output.ESS))
CSV.write("ASMC-AMTM-logZ-SL$(i).csv", DataFrame(logZ = asmc_output.logZ))
CSV.write("ASMC-AMTM-W-SL$(i).csv", DataFrame(W = reduce(vcat, asmc_output.W)))
CSV.write("ASMC_AMTM_SEIR_SCB_samples$(i).csv", df_resampled)
