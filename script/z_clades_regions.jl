"""
    Script to compute the z-scores of mutational fitness discrepancies
    on a either clades or regional (USA-UK) basis
"""

# Import packages
using DataFrames, CSV

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
