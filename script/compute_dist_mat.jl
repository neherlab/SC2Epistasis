""" Script to compute the distance matrix for all positions in the Spike trimer """

# Import required packages
using PdbTool, DelimitedFiles

# Import PDB files
pdbs, af_pdb = SC2Epistasis.read_pdbs("data/ref_seq/Spike.txt")

# Compute the distance matrix for all positions in the Spike trimer
dist_mat = SC2Epistasis.compute_dist_matrix(pdbs, af_pdb)

# Save the distance matrix to a file
writedlm("results/dist_mat.txt", dist_mat)