"""
Script to make histogram comparing z-scores of mutational fitness
discrepancies across clades vs USA-UK regional comparison
"""

# Import packages
using PyPlot

include("usa_uk_energy.jl")

# Load the fitness discrepancies across clades
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
prot = "N" # Select N protein
dfit_prot = delta_fit[delta_fit.prot.==prot, :]

# Compute z-scores for both datasets
rename!(fit_merge, [:delta_fitness => :fit1, :delta_fitness_1 => :fit2, :uncertainty => :std_fit1, :uncertainty_1 => :std_fit2]) # rename columns for compatibility
z_country = SC2Epistasis.z_dfit(fit_merge)
z_clades = SC2Epistasis.z_dfit(dfit_prot)

# Make plot
fig, ax = subplots(figsize=(6.5, 4.5))
ax.hist(z_clades, bins=40, density=true, alpha=0.8, label="Clades", edgecolor="black");
ax.hist(z_country, bins=20, density=true, alpha=0.6, label="USA-UK", edgecolor="black");
ax.legend(fontsize=13);
ax.tick_params(axis="y", labelsize=12);
ax.tick_params(axis="x", labelsize=12);
ax.set_yscale("log");
ax.set_ylabel("Density", fontsize=14);
ax.set_xlabel("Fitness z-score", fontsize=14);
fig.tight_layout();
fig.savefig("results/figures/si/fig_s9_A.pdf");
close(fig)