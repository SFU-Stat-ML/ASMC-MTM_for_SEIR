include("bisection.jl")
include("additionalfns.jl")
include("tmvn_full.jl")
include("asmc_sd.jl")
include("cam_ad.jl")
include("agd.jl")

using DifferentialEquations, Distributions, StatsBase, Random, CSV, DataFrames

N = 22e6
i0 = 2.0
x0 = [N - i0, 0.4 * i0, 0.6 * i0, 0.0]
t = 0:1:136
params = (beta = 0.1447, delta = 1/3.2, gamma = 1/8.5, logphi = log(10))

function seir!(du, u, p, t)
    S, E, I, R = max.(u, 0.0)
    β, δ, γ = p.beta, p.delta, p.gamma
    du[1] = - β * S * I / N
    du[2] = β * S * I / N - δ * E
    du[3] = δ * E - γ * I
    du[4] = γ * I
end

function likelihood_NB_getMean(fit, params)
    mu = [u[3] for u in fit.u]
    t  = fit.t
    return (t = t, mu = mu)
end

function sim_data(n, likelihood, ODEmodel, initial_states, times, params)
    prob = ODEProblem(ODEmodel, initial_states, (times[1], times[end]), params)
    sol = solve(prob, Tsit5(), saveat=times)
    means = likelihood_NB_getMean(sol, params)
    t_vals = means.t
    mu_vals = means.mu
    sampled_indices = round.(Int, range(1, length(mu_vals), length=n))
    phi = exp(params.logphi)
    y = [rand(NegativeBinomial(phi, clamp(phi / (phi + mu_vals[i]), eps(), 1.0))) for i in sampled_indices]
    return (y = y, t = t_vals[sampled_indices])
end

function likelihood(data, theta)
    beta, delta, gamma, phi = theta
    S0 = N - i0; E0 = 0.4 * i0; I0 = 0.6 * i0; R0 = N - (E0 + I0 + S0)
    x_0 = [S0, E0, I0, R0]
    params2 = (beta = beta, delta = delta, gamma = gamma)
    prob = ODEProblem(seir!, x_0, (minimum(t), maximum(t)), params2)
    sol = solve(prob, Tsit5(), saveat=t)
    mu = [u[3] for u in sol.u]
    p = clamp.(phi ./ (phi .+ mu), eps(), 1.0)
    dist = NegativeBinomial.(phi, p)
    return sum(logpdf.(dist, data))
end

function prior(theta)
    beta, delta, gamma, phi = theta
    return logpdf(LogNormal(log(0.1447), 0.1), beta) +
           logpdf(LogNormal(log(1/3.2), 0.07), delta) +
           logpdf(LogNormal(log(1/8.5), 0.05), gamma) +
           logpdf(Gamma(50, 0.2), phi)
end

function target_dis(data, theta)
    return likelihood(data, theta) + prior(theta)
end

reference = prior
init() = [rand(LogNormal(log(0.1447), 0.1)),
          rand(LogNormal(log(1/3.2), 0.07)),
          rand(LogNormal(log(1/8.5), 0.05)),
          rand(Gamma(50, 0.2))]

K = 500
d = 4
M = 10
Mx = 5
n = 10000
b = 0.95
tuning_param = Dict(:eps => 0.5, :eta => 0.95)

Random.seed!(1234)
simulated = sim_data(137, likelihood_NB_getMean, seir!, x0, t, params)
data = simulated.y
CSV.write("simulated_data.csv", DataFrame(y = data))

i = parse(Int, ARGS[1])
seeds = 1001:1050
seed = seeds[i]
println("Running chain $i with seed $seed")
Random.seed!(seed)

asmc_output = asmc(K, d, tuning_param, init, data, likelihood, prior, reference)
theta_init = init()

# Save ASMC diagnostics
R1 = length(asmc_output.particles)
theta_R1 = reduce(vcat, [reshape(p.theta, 1, :) for p in asmc_output.particles[R1]])
W_R1 = asmc_output.W[R1]

#resampled_indices = sample(1:K, Weights(W_R1), K, replace=true)
#resampled_particles = theta_R1[resampled_indices, :]
#df_resampled = DataFrame(resampled_particles, :auto)
#rename!(df_resampled, ["Beta", "Delta", "Gamma", "Phi"])
#CSV.write("ASMC-CAM-rsp-$(i).csv", df_resampled)

#CSV.write("ASMC-CAM-alpha-$(i).csv", DataFrame(alpha_val = asmc_output.alpha))
#CSV.write("ASMC-CAM-ESS-$(i).csv", DataFrame(ESS = asmc_output.ESS))
#CSV.write("ASMC-CAM-logZ-$(i).csv", DataFrame(logZ = asmc_output.logZ))
#CSV.write("ASMC-CAM-W-$(i).csv", DataFrame(W = reduce(vcat, asmc_output.W)))

samples, move_type = cam(data, target_dis, theta_init, n, Mx, M, asmc_output, b)

df = DataFrame(samples, :auto)
rename!(df, ["Beta", "Delta", "Gamma", "Phi"])
df.move_type = move_type
CSV.write("cam_rhat_chain_$(i).csv", df)
