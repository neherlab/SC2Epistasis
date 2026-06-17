""" Script to generate panel D of figure S9 """

# Import packages
using DataFrames, CSV, PyPlot, PyCall
include(joinpath(@__DIR__, "../plot_coup_map.jl"))

# Load inferred couplings
Jtab = CSV.read("results/jcoup_nprot_l1_1em3_3d.csv", DataFrame)

# Compute Frobenius norm of couplings
Jfrob = SC2Epistasis.frob_norm(Jtab)
# Convert the Frobenius norm dictionary into an array
Jmat = SC2Epistasis.Jmat(Jfrob)

# Plot the coupling map
fig, ax = plot_coup_map(Jmat;
    idx_max=75,  # only plot top 75 couplings
    L=421,
    doms_label=["NTD", "SR", "CTD", "Tetramerization"],
    doms_edges=[48:174, 176:206, 247:364, 400:419],
    doms_col=["cyan", "yellow", "red", "green"],
    tick_step=50
)
fig.savefig("results/figures/si/fig_s9_D.pdf")
close(fig)
