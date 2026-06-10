"""
Script to infer the coupling parameters
"""

# Import packages
using DataFrames, CSV, PdbTool

##############################--------------------------------

# Script body

# Initialize structs for inference
include("init_inf.jl")

# Import AlphaFold PDB
af_pdb = PdbTool.parsePdb("data/PDB/Nprot/af_nprot.pdb")

# Map PDB
Nprot = read("data/ref_seq/Nprot.txt", String)
ref_seq = [Nprot]
_ = SC2Epistasis.map_pdb!(af_pdb, ref_seq; mappedTo="data/ref_seq/Nprot.txt")

# Compute distances to define regularization
dist = SC2Epistasis.threedist(optx, data, af_pdb) # compute 3D distances between residues of putative couplings

# L1 regularization
λ1 = 1.0e-3 # regularization strength
l1_reg = SC2Epistasis.threeD_l1(λ1, optx, data, dist) # vector of regularizers

# Infer model parameters
loss = SC2Epistasis.learn(optx, qform, l1_reg; λ2=1.0e-6, epsconv=1.0e-10, verbose=true, maxiter=5000) # optimize

# Write model parameters to dataframe
Jtab = SC2Epistasis.epistatic_table(data, optx)
sort!(Jtab, [:i, :σᵢ, :j, :σⱼ])

# Save couplings to CSV file
CSV.write("results/jcoup_nprot_l1_1em3_3d.csv", Jtab)