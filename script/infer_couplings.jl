"""
Script to infer the coupling parameters
"""

# Import packages
using DataFrames, CSV, PdbTool

##############################--------------------------------

# Script body

# Initialize structs for inference
include("../script/init_inf.jl")

# Read-in PDB files
pdbs, af_pdb = SC2Epistasis.read_pdbs("data/ref_seq/Spike.txt")

# Compute 3D distances between residues of putative couplings
dist = SC2Epistasis.threedist(optx, data, pdbs, af_pdb)

# L1 regularization
λ1 = 1.0e-3 # regularization strength
l1_reg = SC2Epistasis.threeD_l1(λ1, optx, data, dist) # vector of regularizers

# Infer model parameters
loss = SC2Epistasis.learn(optx, qform, l1_reg; λ2=1.0e-6, epsconv=1.0e-10, verbose=true, maxiter=5000) # optimize

# Write model parameters to dataframe
Jtab = SC2Epistasis.epistatic_table(data, optx)
sort!(Jtab, [:i, :σᵢ, :j, :σⱼ])

# Save couplings to CSV file
CSV.write("results/jcoup_l1_1em3_3d.csv", Jtab)