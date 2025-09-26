""" Fit model for different regularization strength and functional choices """

# Import packages
using DataFrames, CSV, PdbTool

# Functions
# Optimize with linear distance specific regularization
function optimize_unif(λ1::Vector{Float64}, optx::Vector{OptX}, qform::Vector{QForm};
    λ2::Float64=1.0e-6,
    epsconv::Float64=1.0e-10,
    verbose::Bool=true,
    maxiter::Int64=5000)

    # Total length of coupling parameters
    tot_num_j = sum([optx[m].num_j for m in eachindex(optx)])

    # Initialize vectors to store results
    e_vec = zeros(length(λ1))
    loss = zeros(length(λ1))
    zf = zeros(length(λ1))

    # Loop over regularization strengths
    for l in eachindex(λ1)

        l1_reg = fill(λ1[l], tot_num_j)
        loss[l] = SC2Epistasis.learn(optx, qform, l1_reg; λ2=λ2, epsconv=epsconv, verbose=verbose, maxiter=maxiter)
        e_vec[l] = SC2Epistasis.energy(optx, qform)
        J = vcat([optx[m].J for m in 1:length(optx)]...)
        zf[l] = sum(abs.(J) .== 0) / length(J)

    end

    return DataFrame(l1=λ1, energy=e_vec, loss=loss, tot_zf=zf)

end

# Optimize with linear distance specific regularization
function optimize_3d_lin(λ1::Vector{Float64}, data::Vector{Data}, optx::Vector{OptX}, qform::Vector{QForm}, dist::Vector{Float64};
    λ2::Float64=1.0e-6,
    epsconv::Float64=1.0e-10,
    verbose::Bool=true,
    maxiter::Int64=5000)

    # Initialize vectors to store results
    e_vec = zeros(length(λ1))
    loss = zeros(length(λ1))
    zf = zeros(length(λ1))

    # Loop over regularization strengths
    for l in eachindex(λ1)

        l1_reg = SC2Epistasis.threeD_l1(λ1[l], optx, data, dist)
        loss[l] = SC2Epistasis.learn(optx, qform, l1_reg; λ2=λ2, epsconv=epsconv, verbose=verbose, maxiter=maxiter)
        e_vec[l] = SC2Epistasis.energy(optx, qform)
        J = vcat([optx[m].J for m in 1:length(optx)]...)
        zf[l] = sum(abs.(J) .== 0) / length(J)

    end

    return DataFrame(l1=λ1, energy=e_vec, loss=loss, tot_zf=zf)

end

# Optimize with sigmoid distance specific regularization
function optimize_3d_sig(λ1::Vector{Float64}, data::Vector{Data}, optx::Vector{OptX}, qform::Vector{QForm}, dist::Vector{Float64};
    scale::Float64=1.0e-2,
    alpha::Float64=0.8,
    d0::Float64=15.0,
    λ2::Float64=1.0e-6,
    epsconv::Float64=1.0e-10,
    verbose::Bool=true,
    maxiter::Int64=5000)

    # Initialize vectors to store results
    e_vec = zeros(length(λ1))
    loss = zeros(length(λ1))
    zf = zeros(length(λ1))

    # Loop over regularization strengths
    for l in eachindex(λ1)
        l1_reg = SC2Epistasis.sigmoid_l1(scale * λ1[l], λ1[l], optx, data, dist; α=alpha, d0=d0)
        loss[l] = SC2Epistasis.learn(optx, qform, l1_reg; λ2=λ2, epsconv=epsconv, verbose=verbose, maxiter=maxiter)
        e_vec[l] = SC2Epistasis.energy(optx, qform)
        J = vcat([optx[m].J for m in 1:length(optx)]...)
        zf[l] = sum(abs.(J) .== 0) / length(J)
    end

    return DataFrame(l1=λ1, energy=e_vec, loss=loss, tot_zf=zf)

end

##############################--------------------------------

# Script body

# Initialize structs for inference
include("../script/init_inf.jl")

dist = SC2Epistasis.threedist(optx, data, pdbs, af_pdb) # compute 3D distances between residues of putative couplings

# Array of regularization strengths
λ1 = [1.0e-6, 1.0e-5, 5.0e-5, 1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2]

# Fit model for different regularization strategies and strengths
fit_df_unif = optimize_unif(λ1, optx, qform; epsconv=1.0e-8)
fit_df_3d = optimize_3d_lin(λ1, data, optx, qform, dist; epsconv=1.0e-8)
fit_df_3dsig = optimize_3d_sig(λ1, data, optx, qform, dist; scale=1.0e-2, alpha=0.8, d0=15.0, epsconv=1.0e-8)

# Write dataframes to file
CSV.write("results/model_fit/fit_unif.csv", fit_df_unif)
CSV.write("results/model_fit/fit_3d.csv", fit_df_3d)
CSV.write("results/model_fit/fit_3dsig.csv", fit_df_3dsig)