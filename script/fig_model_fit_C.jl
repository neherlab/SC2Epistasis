"""
    Script to generate panel C of figure 4
"""

# Import packages
using DataFrames, CSV, PyPlot

# Load coupling parameters
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)

# Plot and save histogram of coupling parameters
hist(Jtab.J, bins=20, density=false, alpha=0.8, edgecolor="black")
yscale("log")
xlabel(L"J_{ij}(σᵢ,σⱼ)", fontsize=14)
ylabel("Frequency", fontsize=14)
tight_layout()
savefig("results/figures/coup_histogram.pdf", dpi=500)
close("all")