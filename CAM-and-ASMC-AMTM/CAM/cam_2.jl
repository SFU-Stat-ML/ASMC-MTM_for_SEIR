# CAM Algorithm correct
include("agd.jl")
include("local_cov.jl")
include("weightedcov.jl")

using LinearAlgebra, StatsBase, Distributions, Random, RCall

function make_posdef_cov(Σ; jitter=1e-6)
    Σ = Matrix(Symmetric(Σ))
    d = size(Σ, 1)

    for k in 0:8
        Σj = Symmetric(Σ + (jitter * 10.0^k) * I(d))
        try
            cholesky(Σj)
            return Σj
        catch e
            if !(e isa PosDefException)
                rethrow(e)
            end
        end
    end

    error("Could not make covariance positive definite")
end

function extract_auxiliary_info(asmc_output)
    R1 = length(asmc_output.particles)
    # Use the final target distribution to build the auxiliary covariance matrix
    theta_mat = Matrix(hcat([p.theta for p in asmc_output.particles[R1]]...)')
    weights = asmc_output.W[R1]
    cov_aux = weighted_covariance_matrix(theta_mat, weights)
    
    # NEW: Use agd_all to draw the auxiliary variable from ALL intermediate distributions
    z = agd_all(asmc_output, 1)[1]
    
    return z, cov_aux
end

function cam(data, target_dis, theta_old, n, Mx, M, asmc_output)

    d = length(theta_old)
    lower = zeros(d)

    sample_matrix = Matrix{Float64}(undef, n, d)
    sample_matrix[1, :] = theta_old
    
    acceptances = 0
    
    for r in 2 : n
        # 1. Update Adaptive Covariance dynamically
        cov_local = make_posdef_cov(local_cov(sample_matrix, d, 0.95, r))

        # 2. Extract NEW auxiliary variable for THIS iteration
        z, cov_aux_raw = extract_auxiliary_info(asmc_output)
        cov_aux = make_posdef_cov(cov_aux_raw)

        # 3. Generate Candidates
        ym = Matrix{Float64}(undef, M, d)
        log_um = Vector{Float64}(undef, M)
        T_fwd_vec = Vector{Float64}(undef, M)

        for m in 1 : M
            if m <= Mx
                proposal_mean = theta_old
                sigma = Matrix(cov_local)
                y_m = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $proposal_mean, sigma = $sigma, lower = $lower))""")
                reverse_mean = y_m
            else
                proposal_mean = z
                sigma = Matrix(cov_aux)
                y_m = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $proposal_mean, sigma = $sigma, lower = $lower))""")
                reverse_mean = z
            end

            ym[m, :] = y_m

            T_fwd_vec[m] = rcopy(Float64, R"""as.numeric(tmvtnorm::dtmvnorm($y_m, mean = $proposal_mean, sigma = $sigma, lower = $lower, log = TRUE))""")

            log_q_back = rcopy(Float64, R"""as.numeric(tmvtnorm::dtmvnorm($theta_old, mean = $reverse_mean, sigma = $sigma, lower = $lower, log = TRUE))""")

            log_um[m] = target_dis(data, y_m) + log_q_back
        end

        # Normalize weights and sample J
        logsum_um = logsum(log_um)
        N_um = log_um .- logsum_um
        um = exp.(N_um)

        J = sample(1:M, Weights(um))
        theta_new = ym[J, :]
        T_fwd = T_fwd_vec[J]

        # 4. Generate Reverse Samples
        xm_star = Matrix{Float64}(undef, M, d)  
        log_um_star = Vector{Float64}(undef, M)    
        T_rev_vec = Vector{Float64}(undef, M)

        for m in 1:M
            if m <= Mx
                proposal_mean = theta_new
                sigma = Matrix(cov_local)

                if m == J
                    x_star = theta_old
                else
                    x_star = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $proposal_mean, sigma = $sigma, lower = $lower))""")
                end

                back_mean = x_star
            else
                proposal_mean = z
                sigma = Matrix(cov_aux)

                if m == J
                    x_star = theta_old
                else
                    x_star = rcopy(Vector{Float64}, R"""as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = $proposal_mean, sigma = $sigma, lower = $lower))""")
                end

                back_mean = z
            end

            xm_star[m, :] = x_star

            # Eq. 15: density of the reverse candidate, starting from theta_new.
            T_rev_vec[m] = rcopy(Float64, R"""as.numeric(tmvtnorm::dtmvnorm($x_star, mean = $proposal_mean, sigma = $sigma, lower = $lower, log = TRUE))""")

            # Eq. 14: density of returning to theta_new from this candidate.
            log_q_back = rcopy(Float64, R"""as.numeric(tmvtnorm::dtmvnorm($theta_new, mean = $back_mean, sigma = $sigma, lower = $lower, log = TRUE))""")

            log_um_star[m] = target_dis(data, x_star) + log_q_back
        end

        # Normalize reverse weights
        logsum_um_star = logsum(log_um_star)
        N_um_star = log_um_star .- logsum_um_star

        # 5. Compute the Acceptance Probability
        log_pi_y = target_dis(data, theta_new)
        log_pi_x = target_dis(data, theta_old)
        
        N_um_J = N_um[J]
        N_um_star_J = N_um_star[J]
        T_rev = T_rev_vec[J]

        num = N_um_star_J + log_pi_y + T_rev
        den = N_um_J + log_pi_x + T_fwd

        log_ratio = min(0.0, num - den) 

        # 6. Accept/Reject and UPDATE STATE
        if log(rand()) < log_ratio
            theta_old = theta_new   
            sample_matrix[r, :] = theta_new
            acceptances += 1
        else
            sample_matrix[r, :] = theta_old
        end
    end
    
    println("CAM Acceptance Rate: ", round(acceptances / (n - 1) * 100, digits=2), "%")
    return sample_matrix
end