""" Script to generate panel A of figure S1 of the SI """

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
cp_plt = [("20I", "21J"), ("20I", "21L"), ("21J", "21K"), ("21K", "23I")]

# Make plots
fig, _ = plot_raw_zscores_bands(cp_plt, z_dict, s_dict, clade_diff; x_space=55, rasterized=false)

# Save figure
fig.savefig("results/figures/si/fig_s1_A.pdf", dpi=300)
close(fig)