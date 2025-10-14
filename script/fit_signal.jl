""" 
    Script to compute the raw epistatic signal between all clade pairs
"""

# Import packages
using DataFrames, CSV, PdbTool, Combinatorics, JLD2, PyPlot

# Import dataframes with fitness discrepancies and clade mismatches
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Protein of interest
prot = "S" # Spike protein

# Protein specific fitness and clade mismatches
delta_fit_prot = delta_fit[delta_fit.prot.==prot, :]
clade_diff_prot = clade_diff[clade_diff.prot.==prot, :]

# Filtering fitness dataframe
idx_nsc = findall(x -> x[end] != '*', delta_fit_prot.aa_mut) #excluding stop codons
delta_fit_prot = delta_fit_prot[idx_nsc, :] #excluding stop codons

# List of Clades
clades = unique(vcat(delta_fit_prot.clade1, delta_fit_prot.clade2))

# Compute z-scores and sites for all clade pairs
z_dict, s_dict = SC2Epistasis.raw_signal(clades, delta_fit_prot, clade_diff_prot)

# Save results
fit_epi = (z_score=z_dict, sites=s_dict)
@save "results/fit_zscore_epi.jld2" fit_epi