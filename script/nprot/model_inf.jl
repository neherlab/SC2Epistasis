""" Fit model for different regularization strength and functional choices """

# Import packages
using DataFrames, CSV, PdbTool, PyPlot

##############################--------------------------------

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

##############################--------------------------------

# Script body

include("usa_uk_energy.jl") # e0 is defined here

# Initialize structs for inference
include("init_inf.jl")

# Import AlphaFold PDB
af_pdb = PdbTool.parsePdb("data/PDB/Nprot/af_nprot.pdb")

# Map PDB
Nprot = read("data/ref_seq/Nprot.txt", String)
ref_seq = [Nprot]
_ = SC2Epistasis.map_pdb!(af_pdb, ref_seq; mappedTo="data/ref_seq/Nprot.txt")

# Compute distances to define regularization
dist = SC2Epistasis.threedist(optx, data, af_pdb) # compute 3D distances between residues of putative couplings

# Array of regularization strengths
λ1 = [1.0e-6, 1.0e-5, 5.0e-5, 1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2]

# Fit model for different regularization strategies and strengths
fit_df_unif = optimize_unif(λ1, optx, qform; epsconv=1.0e-10)
fit_df_3d = optimize_3d_lin(λ1, data, optx, qform, dist; epsconv=1.0e-10)

# Plot energy trends
labels = ["Uniform", "3d"] # labels for legend
linestyles = ["-", "--"] # line styles for different models
shades = [get_cmap("Blues")(0.4 + 0.2 * i) for i in 0:1] # colors for different models
fit_df_vec = [fit_df_unif, fit_df_3d]

for n in eachindex(fit_df_vec)
    fit_df = fit_df_vec[n]
    plot(fit_df.l1, fit_df.energy, marker=".", markersize=5.0,
        linestyle=linestyles[n], color=shades[n], label=labels[n])
end

hlines(e0, fit_df_vec[end].l1[1], fit_df_vec[end].l1[end], color="black", linestyle="-", label="Estimated noise")
xscale("log")
xlabel("λ₁ (Sparsity Regularization)", fontsize=14)
ylabel("Energy (Mean Squared Error)", fontsize=14)
legend(loc="best", fontsize=12)
tight_layout()
savefig("results/figures/si/fig_s9_B.pdf")
close("all")