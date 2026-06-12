""" Script to generate panel A of figure S2 of the SI """

# Import packages
using DataFrames, CSV, JLD2, PyPlot, PyCall
include("z_plot_raw_zscores_bands.jl")

################################---------------------------------

# Load results
@load "results/fit_zscore_epi.jld2"
z_dict = fit_epi.z_score
s_dict = fit_epi.sites

# Load clade differences
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)
prot = "S" # Spike protein
clade_diff = clade_diff[clade_diff.prot.==prot, :]

# Define clade pairs to be plotted
cp_plt = [("21K", "22B"), ("21K", "23A"), ("21L", "22B"), ("21L", "23A")]

# Make plots
fig, _ = plot_raw_zscores_bands(cp_plt, z_dict, s_dict, clade_diff; x_space=55, rasterized=false)

# Save figure
fig.savefig("results/figures/si/fig_s2_A.pdf", dpi=300)
close(fig)