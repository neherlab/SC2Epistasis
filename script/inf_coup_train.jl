""" Script to infer the couplings in the framework of the train test split for validation
"""

# Import packages
using DataFrames, CSV, JSON

# Read dataframe with mismatches between clade founders
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Read dataframe with fitness discrepancies
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)

# Select Spike protein
prot = "S"
dfit_prot = delta_fit[delta_fit.prot.==prot, :]
cdiff_prot = clade_diff[clade_diff.prot.==prot, :]

# Training clades
clades = unique(vcat(dfit_prot.clade1, dfit_prot.clade2))
train_clades = clades[1:end-1]
dfit_train = filter(r -> r.clade1 in train_clades && r.clade2 in train_clades, dfit_prot)

# Mutations list
muts_prot = unique(dfit_train.aa_mut)
idx_sort = sortperm(map(x -> parse(Int64, x[2:end-1]), muts_prot)) # Sort mutations according to sites
muts_prot = muts_prot[idx_sort]

# Compute structs for inference
data, optx, qform = SC2Epistasis.init_all(muts_prot, dfit_train, cdiff_prot)

# Define distance dependent regularization

# Read-in PDB files
pdbs, af_pdb = SC2Epistasis.read_pdbs("data/ref_seq/Spike.txt")

# Compute 3D distances between residues of putative couplings
dist = SC2Epistasis.threedist(optx, data, pdbs, af_pdb)

# L1 regularization
λ1 = 1.0e-3 # regularization strength
l1_reg = SC2Epistasis.threeD_l1(λ1, optx, data, dist) # vector of regularizers

# Infer coupling parameters
loss = SC2Epistasis.learn(optx, qform, l1_reg; λ2=1.0e-6, epsconv=1.0e-10, verbose=true, maxiter=5000) # optimize

# Write couplings to dataframe and save
Jtab = SC2Epistasis.epistatic_table(muts_prot, optx)
CSV.write("results/jcoup_train.csv", Jtab)