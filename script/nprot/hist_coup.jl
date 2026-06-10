"""
    Script to generate panel C of figure S9
"""

# Import packages
using DataFrames, CSV, PyPlot

# Load coupling parameters
Jtab = CSV.read("results/jcoup_nprot_l1_1em3_3d.csv", DataFrame)

# Plot and save histogram of coupling parameters
hist(Jtab.J, bins=21, density=false, alpha=0.8, edgecolor="black");
yscale("log")
xlabel("Coupling parameters " * L"J_{ij}(σᵢ,σⱼ)", fontsize=14)
ylabel("Frequency", fontsize=14)
tight_layout()
savefig("results/figures/si/fig_s9_C.pdf")
close("all")