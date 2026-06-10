"""
Script to infer the coupling parameters for uniform regularization
"""

# Import packages
using DataFrames, CSV, PdbTool

##############################--------------------------------

# Script body

# Initialize structs for inference
include("../script/init_inf.jl")

# L1 regularization
λ1 = 1.0e-3 # regularization strength
tot_num_j = sum([optx[m].num_j for m in eachindex(muts_prot)]) # total number of coupling parameters
l1_reg = fill(λ1, tot_num_j) # vector of regularizers

# Infer model parameters
loss = SC2Epistasis.learn(optx, qform, l1_reg; λ2=1.0e-6, epsconv=1.0e-10, verbose=true, maxiter=5000) # optimize

# Write model parameters to dataframe
Jtab = SC2Epistasis.epistatic_table(data, optx)
sort!(Jtab, [:i, :σᵢ, :j, :σⱼ])

# Save couplings to CSV file
CSV.write("results/jcoup_l1_1dot5em4_unif.csv", Jtab)