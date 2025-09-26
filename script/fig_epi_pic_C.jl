"""
    Script to generate panel C of figure 2
"""

# Import packages
using DataFrames, CSV, PyPlot

# Load and merge the USA-UK mutational fitness dataframes
include("usa_uk_ener.jl")

# Load the fitness discrepancies across clades
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
prot = "S" # Select Spike protein
dfit_prot = delta_fit[delta_fit.prot.==prot, :]

# Compute z-scores for both datasets
rename!(fit_merge, [:delta_fitness => :fit1, :delta_fitness_1 => :fit2, :uncertainty => :std_fit1, :uncertainty_1 => :std_fit2]) # rename columns for compatibility
z_country = SC2Epistasis.z_dfit(fit_merge)
z_clades = SC2Epistasis.z_dfit(dfit_prot)

# Make histogram comparison
hist(z_clades, bins=40, density=true, alpha=0.8, label="Clades", edgecolor="black")
hist(z_country, bins=20, density=true, alpha=0.6, label="USA-UK", edgecolor="black")
legend(fontsize=12)
yscale("log")
xlabel("Normalized fitness discrepancy", fontsize=14);
ylabel("Density", fontsize=14);
tight_layout();
savefig("results/figures/fig_epi_pic_C.pdf", dpi=500);
close("all");
