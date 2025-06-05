# Analysis using CAM
using StatsBase, DataFrames, Plots

# Discard burn-in
burnin = 200
posterior_samples = samples[burnin+1:end, :]

# Compute summaries
means = mean(posterior_samples, dims=1)
ci_lower = mapslices(x -> quantile(x, 0.025), posterior_samples, dims=1)
ci_upper = mapslices(x -> quantile(x, 0.975), posterior_samples, dims=1)

println("Posterior means: ", means)
println("95% Credible Intervals:")
for i in 1:d
    println("Param $i: [", ci_lower[i], ", ", ci_upper[i], "]")
end

using Plots

param_names = ["β", "δ", "γ", "ϕ"]
for i in 1:d
    plt = plot(samples[:, i], label=param_names[i], xlabel="Iteration", ylabel="Value", title="Trace plot of $(param_names[i])")
    savefig(plt, "traceplot_param_$(param_names[i]).png")
end
