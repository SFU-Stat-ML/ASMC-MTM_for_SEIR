# log-sum-exponential evaluation of log(sum(w)) to avoid
# numerical overflow and underflow

function logsum(logw)
    logmax = maximum(logw)
    log(sum(exp.(logw.-logmax))) + logmax    
end

# relative conditional effective sample size
function rCESS(W, u, a, eta)
    logw = a .* u
    exp(2 * logsum(log.(W) .+ logw) - logsum(log.(W) .+ 2 .* logw)) - eta
end

# Effective sample size
function rESS(logW)
    K = length(logW)
    logWmax = maximum(logW)
    logrESS = -(2 * logWmax + log(sum(exp.(2 .* logW .- 2 * logWmax)))) - log(K)
    return exp(logrESS)
end

# systematic resampling
function systematic_resample(W)
    K = length(W)
    U = rand() / K .+ (0:(K-1)) ./ K
    Wsum = cumsum(W)
    N1 = zeros(Int, K)
    j = 1

    # re-sample loop
    for i in 1:K
        found = false
        while !found
            if U[i] > Wsum[j]
                j = j + 1
            else
                found = true
            end
        end
        N1[i] = j
    end
    return N1
end
