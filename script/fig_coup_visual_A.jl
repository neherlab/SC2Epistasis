""" Script to generate panel A of figure 5 """

# Import packages
using DataFrames, CSV, PyPlot, PyCall
include(joinpath(@__DIR__, "plot_coup_map.jl"))

# Load inferred couplings
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)

# Compute Frobenius norm of couplings
Jfrob = SC2Epistasis.frob_norm(Jtab)
# Convert the Frobenius norm dictionary into an array
Jmat = SC2Epistasis.Jmat(Jfrob)

# Plot the coupling map
fig, ax = plot_coup_map(Jmat; idx_max=400) # only plot top 400 couplings
savefig("results/figures/fig_coup_visual_A.pdf")
close(fig)
