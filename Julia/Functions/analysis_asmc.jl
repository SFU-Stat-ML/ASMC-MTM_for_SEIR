using StatsBase, DataFrames, Plots

# Existing diagnostics extraction...
R1 = length(asmc_output.particles)
function get_theta_matrix(particles_at_r)
    return hcat([p.theta for p in particles_at_r]...)'
end

theta_R1  = get_theta_matrix(asmc_output.particles[R1])

W_R1 = asmc_output.W[R1]
resampled_indices = sample(1:size(theta_R1, 1), Weights(W_R1), size(theta_R1, 1), replace=true)
resampled_particles = theta_R1[resampled_indices, :]

# Parameter names
param_names = ["β", "δ", "γ", "φ"]

# Summary DataFrame
summary_stats = DataFrame(
    Parameter = String[],
    Mean = Float64[],
    CI_lower = Float64[],
    CI_upper = Float64[]
)


# Compute mean and 95% CI for each parameter column
for j in 1:size(resampled_particles, 2)
    samples = resampled_particles[:, j]
    μ = mean(samples)
    ci = quantile(samples, [0.025, 0.975])
    push!(summary_stats, (param_names[j], μ, ci[1], ci[2]))
end

# View summary
println(summary_stats)

AASMC_alpha = DataFrame(alpha_val = asmc_output.alpha)
AASMC_ESS   = DataFrame(ESS = asmc_output.ESS)
AASMC_logZ  = DataFrame(logZ = asmc_output.logZ)
AASMC_W     = DataFrame(W = reduce(vcat, asmc_output.W))

# ---- Plotting diagnostics ----
alpha = AASMC_alpha.alpha_val
ESS = AASMC_ESS.ESS

iters = 1:length(alpha)

# Plot ESS over iterations
plot_ess = plot(
    iters,
    ESS,
    seriestype = :scatter,
    marker = (:circle, 4),
    line = (:dot, :black),
    xlabel = "Iteration",
    ylabel = "ESS",
    title = "ESS Over Iterations",
    legend = false
)

# Plot alpha over iterations
plot_alpha = plot(
    iters,
    alpha,
    seriestype = :scatter,
    marker = (:circle, 4),
    line = (:dot, :blue),
    xlabel = "Iteration",
    ylabel = "Alpha",
    title = "Alpha Over Iterations",
    legend = false
)

# Display side by side
plot(plot_ess, plot_alpha, layout = (1, 2), size=(800, 600))
savefig("asmc_diagnostics.png")
