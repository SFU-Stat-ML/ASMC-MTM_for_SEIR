# CAM Algorithm correct
include("agd.jl")
include("tmvn_full.jl")
#include("local_cov.jl")
#include("weightedcov.jl")


using LinearAlgebra, StatsBase, Distributions, Random

# target_dis = target Distribution
# theta_old = current state
# r = rth iteration
# n = number of iterations
# M = total number of proposals
# Mz = total number of auxiliary proposals
# Mx = total number of local proposals and M = Mz + Mx
# cov_local = local covariance matrix using Adaptive Metropolis
# cov_aux = auxiliary covariance matrix using weighted_cov from ASMC samples

function extract_auxiliary_info(asmc_output)
    R1 = length(asmc_output.particles)
    #theta_mat = hcat([p.theta for p in asmc_output.particles[R1]]...)'
    theta_mat = Matrix(hcat([p.theta for p in asmc_output.particles[R1]]...)')
    weights = asmc_output.W[R1]
    cov_aux = weighted_covariance_matrix(theta_mat, weights)
    z = agd(asmc_output, 1)[1]
    return z, cov_aux
end

function cam(data, target_dis, theta_old, n, Mx, M, asmc_output, b)

    d = length(theta_old)
    lower = zeros(d)
    #lower = 0.01 * ones(d)

    # 4. generate auxiliary variables
    # For simplicity, a fixed auxiliary variable z is used throughout all CAM iterations.
    # This is a simplified version of CAM where AGD is not refreshed every iteration.
    z, cov_aux = extract_auxiliary_info(asmc_output)
    cov_aux = Symmetric(cov_aux)

    sample_matrix = Matrix{Float64}(undef, n, d)
    sample_matrix[1, :] = theta_old

    move_type = Vector{Int}(undef, n) ##
    move_type[1] = 0  # First sample is just the initial state ##

    function empirical_covariance(samples::Matrix{Float64})
        μ = mean(samples, dims=1)
        centered = samples .- μ
        return centered' * centered / (size(samples, 1) - 1)
    end

    for r in 2 : n

        if r <= 2 * d || rand() >= b
            cov_local = (0.1^2 / d) * I(d)
        else
            Σ = empirical_covariance(sample_matrix[1:r-1, :])
            Σ = (Σ + Σ') / 2
            Σ += 1e-6 * I(d)
            cov_local = (2.38^2 / d) * Σ
        end

        # 5.
        ym = Matrix{Float64}(undef, M, d)
        log_um = Vector{Float64}(undef, M)
        N_um = Vector{Float64}(undef, M)
        T_fwd_vec = Vector{Float64}(undef, M)

        for m in 1 : M
            # Generate candidate ym
            if m <= Mx
                fwd_sampler = TruncatedMVNSampler(MvNormal(theta_old, cov_local), lower)
            else
                fwd_sampler = TruncatedMVNSampler(MvNormal(z, cov_aux), lower)
            end

            ym[m, :] = rand(fwd_sampler)  

            # compute weights um
            logp = target_dis(data, ym[m, :])
            #um[m] = sqrt(exp(logp))
            log_um[m] = 0.5 * logp
            T_fwd_vec[m] = tmvn_pdf(ym[m, :], fwd_sampler)
        end

        # normalized weights
        #um ./= sum(um)
        logsum_um = logsum(log_um)
        
        N_um = log_um .- logsum_um
        um = exp.(log_um .- logsum_um)

        # sample J
        J = sample(1:M, Weights(um))

        move_type[r] = J > Mx ? 1 : 0  # 1 = auxiliary, 0 = local ##


        # set y = yJ
        theta_new = ym[J, :]

        T_fwd = T_fwd_vec[J]


        # reverse sample generation
        xm_star = Matrix{Float64}(undef, M, d)  # reverse samples
        log_um_star = Vector{Float64}(undef, M)    # reverse weights
        N_um_star = Vector{Float64}(undef, M)    # reverse weights
        T_rev_vec = Vector{Float64}(undef, M)

        for m in 1:M
            if m == J
             xm_star[m, :] = theta_old
             T_rev_vec[m] = T_fwd
             logp_star = target_dis(data, xm_star[m, :])
             #um_star[m] = sqrt(exp(logp_star))
             log_um_star[m] = 0.5 * logp_star
            else
                if m <= Mx
                    rev_sampler = TruncatedMVNSampler(MvNormal(theta_new, cov_local), lower)
                else
                    rev_sampler = TruncatedMVNSampler(MvNormal(z, cov_aux), lower)
                end

                xm_star[m, :] = rand(rev_sampler)
                # reverse weights
                logp_star = target_dis(data, xm_star[m, :])
                #um_star[m] = sqrt(exp(logp_star))
                log_um_star[m] = 0.5 * logp_star
                T_rev_vec[m] = tmvn_pdf(xm_star[m, :], rev_sampler)
            end
        end

        # x_star_J = xm_star[J, :]
        T_rev = T_rev_vec[J]

        # normalized reverse weights
        #um_star ./= sum(um_star)
        logsum_um_star = logsum(log_um_star)
        #um_star = exp.(log_um_star .- logsum_um_star)
        N_um_star = log_um_star .- logsum_um_star

        # compute the acceptance probability
        log_pi_y = target_dis(data, theta_new)
        log_pi_x = target_dis(data, theta_old)
        N_um_J = N_um[J]
        N_um_star_J = N_um_star[J]

        num = 2 * N_um_star_J + log_pi_y + log(T_rev)
        den = 2 * N_um_J + log_pi_x + log(T_fwd)

        log_ratio =  min(0, num - den) 

        #ratio = min(1.0, num / den)

        # accept/reject
        decision = log(rand()) < log_ratio
        
        if decision
            sample_matrix[r, :] = theta_new
        else
            sample_matrix[r, :] = theta_old
        end
    end
    return sample_matrix, move_type ##
end
