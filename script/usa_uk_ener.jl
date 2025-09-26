"""
Script to compute the baseline noise level by comparing
fitness estimates from USA and UK subsets
"""

# Import packages
using DataFrames, CSV, LinearAlgebra

# Load the data
aamut_fit_usa = DataFrame(CSV.File("data/usa_uk_fitness/aamut_fitness_usa.csv.gz"))
aamut_fit_uk = DataFrame(CSV.File("data/usa_uk_fitness/aamut_fitness_uk.csv.gz"))

# Filter data
prot = "S" # Select Spike protein
fit_usa = aamut_fit_usa[aamut_fit_usa.gene.==prot, :]
fit_uk = aamut_fit_uk[aamut_fit_uk.gene.==prot, :]

cnt_thr = 10.0 # Minimum count threshold
fit_usa = fit_usa[fit_usa.predicted_count.>=cnt_thr, :]
fit_uk = fit_uk[fit_uk.predicted_count.>=cnt_thr, :]

# Retain non-synonymous and nonsense mutations
fit_usa = fit_usa[fit_usa.clade_founder_aa.!=fit_usa.mutant_aa, :]
fit_uk = fit_uk[fit_uk.clade_founder_aa.!=fit_uk.mutant_aa, :]

# Merge dataframes
merge_cols = [:cluster, :gene, :clade_founder_aa, :mutant_aa, :aa_site, :aa_mutation]
ret_cols = vcat(merge_cols, [:predicted_count, :actual_count, :tau_squared, :delta_fitness, :uncertainty])
fit_merge = innerjoin(fit_usa[:, ret_cols], fit_uk[:, ret_cols], on=merge_cols, makeunique=true)

# Compute baseline energy
dfit2 = (fit_merge.delta_fitness .- fit_merge.delta_fitness_1) .^ 2
w = 1 ./ (fit_merge.uncertainty .^ 2 + fit_merge.uncertainty_1 .^ 2)
w ./= sum(w)
e0 = dot(dfit2, w)
