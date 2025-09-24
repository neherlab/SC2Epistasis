""" Script to initialize the structs used for inference """

# Import packages
using DataFrames, CSV

# Import fitness discrepancies and clade differences
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Select protein
prot = "S"
delta_fit_prot = delta_fit[delta_fit.prot.==prot, :]
clade_diff_prot = clade_diff[clade_diff.prot.==prot, :]

# Mutations list
muts_prot = unique(delta_fit_prot[:, :aa_mut])
# Sort mutations according to sites
idx_sort = sortperm(map(x -> parse(Int64, x[2:end-1]), muts_prot))
muts_prot = muts_prot[idx_sort]

# Initialize structs for inference
data, optx, qform = SC2Epistasis.init_all(muts_prot, delta_fit_prot, clade_diff_prot)

# Compute distances to define regularization
pdbs, af_pdb = SC2Epistasis.read_pdbs() # read PDB files