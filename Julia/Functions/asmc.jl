include("bisection.jl")
include("mcmcmove.jl")
include("additionalfns.jl")
include("weightedcov.jl")

using StatsBase, Statistics

# ASMC Function
function asmc(K, d, b, tuning_param, init, data, likelihood, prior, reference)

    # initialize storage lists
    alpha = Vector{Float64}()         
    ESS = Vector{Float64}()            
    logZ = Vector{Float64}()           
    particles = Vector{Vector{Any}}()  
    W = Vector{Vector{Float64}}()      

    # check if adaptive annealing scheme is being used
    adaptive_alpha = haskey(tuning_param, :eta)

    # initialize the iterations
    r = 1

    # initialize the log-normalizing constant
    push!(logZ, 0.0)
    

    # initialize the annealing parameter
    if adaptive_alpha
        push!(alpha, 0.0) # ✔ 
        alphaDiff = 0
    else
        alpha = copy(tuning_param[:alpha])
    end

    # initialize particles
    push!(particles, [ (theta = init(), accept = false) for k in 1:K ])

    # initialize the weighted covariance matrix
    weighted_cov = Matrix{Float64}(I, d, d)

    # normalized weights
    push!(W, fill(1 / K, K))   
    logW = log.(W[r])         

    # unnormalized weights
    w = fill(1.0, K)          
    logw = fill(0.0, K) # debugged upto here

    # begin iterative process
    while alpha[r] < 1
    
        print("iteration: ", r, "\n")

        r += 1 # ✔

        # evaluate loglikelihood for updating alpha
        u = fill(0.0, K)

        # evaluate the loglikelihood for each particle using the current and previous data
        u = [likelihood(data, particles[r - 1][k].theta) for k in 1:K] # ✔

       if adaptive_alpha
        alphaDiff = bisection(0.0, 1.0, W[r - 1], u, tuning_param[ :eta])
        push!(alpha, alpha[r - 1] + alphaDiff)
       else
        alphaDiff = alpha[r] - alpha[r - 1] # ✔
       end

        #if r > length(tuning_param[:alpha])
        #error("Index $r exceeds the length of the fixed alpha schedule. Make sure tuning_param[:alpha] has enough entries.")
        #end
        #push!(alpha, tuning_param[:alpha][r])
        
       print("annealing parameter: ", alpha[r], "\n")

       if alpha[r] > 1
        alpha[r] = 1
        alphaDiff = 1 - alpha[r - 1]
       end

       # updated MCMCmove function
       push!(particles, [
        mcmcmove(p, data, likelihood, prior, alpha[r], r, d, b, weighted_cov)
        for p in particles[r - 1]
        ]) # ✔

       # converting particles[r] list to a matrix
       theta_values = [p.theta for p in particles[r - 1]]
       #theta_matrix = reduce(vcat, [theta' for theta in theta_values])
       theta_matrix = Matrix(hcat(theta_values...)')  # always gives K × d matrix
    
       # compute ESS
       log_incremental_w = alphaDiff .* u
       logw .= log_incremental_w .+ logw

       logmax = maximum(logw)
       push!(logZ, logZ[r - 1] + logsum(log_incremental_w .+ logW)) # ✔

       push!(W, exp.(logw .- logmax))
       W[r] ./= sum(W[r])         
       logW = log.(W[r])
       push!(ESS, rESS(logW))

       # resample if ESS is below the threshold
       if ESS[r - 1] < tuning_param[:eps]
        println("Resample: ESS = ", ESS[r - 1], "\n") # ✔
    
        ancestors = systematic_resample(W[r])
        particles[r] = particles[r][ancestors]
    
        W[r] = fill(1/K, K)
        logW = log.(W[r])
        w = fill(1.0, K)
        logw = fill(0.0, K)
       end
    
       # calculate the weighted covariance matrix
       weighted_cov = weighted_covariance_matrix(theta_matrix, W[r - 1])
 
    end

    return (
        particles = particles,
        alpha = alpha,
        ESS = ESS,
        logZ = logZ,
        W = W
    )

end
