""" Script to generate plots related to the RSA analysis. """

# Import necessary packages
using DataFrames, CSV, JLD2, PyPlot, PyCall, PdbTool, DelimitedFiles

# Load RSA data
df_rsa = CSV.read("data/7KRR_per_residue_RSA.csv", DataFrame)

# Assign unique RSA to residues by taking the max RSA across chains
rsa_max = combine(groupby(df_rsa, :residue), :RSA => maximum => :RSA)

# Load Shannon entropy data
site_shannon_ent = CSV.read("data/nextstrain_staging_nextclade_sars-cov-2_diversity.tsv", DataFrame)
rename!(site_shannon_ent, "position" => "residue")

# Populate missin sites with 0 entropy
all_sites = collect(1:1273)
missing_sites = setdiff(all_sites, site_shannon_ent.residue)
df_missing = DataFrame(residue=missing_sites, entropy=zeros(length(missing_sites)))

# Combine the two DataFrames
df_ent = vcat(site_shannon_ent[:, 2:end], df_missing)
sort!(df_ent, :residue)

# Merge RSA and entropy data
df_rsa_ent = innerjoin(rsa_max, df_ent, on=:residue)

# Load coupling parameters
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)
j_res = unique(Jtab.j) # background mismatches

# Find residues which are background mismatches
idx_jmis = findall(x -> x in j_res, df_rsa_ent.residue)

# Compute Frobenius norm of the coupling parameters for each residue pair
Jfrob = SC2Epistasis.frob_norm(Jtab)
Jmat = SC2Epistasis.Jmat(Jfrob)
# Frobenius norm dataframe
Jf_df = DataFrame(i=map(x -> x[1][1], Jmat), j=map(x -> x[1][2], Jmat), J=map(x -> x[2], Jmat))
Jf_n0 = Jf_df[Jf_df.J .> 0, :] #select non-zero couplings
# For each background mismatch j, find the maximum coupling value across all residues i
coup_max = combine(
    groupby(Jf_n0, :j),
    :J => maximum => :J
)
rename!(coup_max, :j => :residue) # rename column j to residue for merging
j_rsa_max = innerjoin(coup_max, rsa_max, on=:residue)

# Make figure
fig, ax = subplots(1, 2, figsize=(10, 5))

# RSA value histogram
ax[1].hist(df_rsa_ent.RSA, bins=20, alpha=0.7, density=true, label="All sites")
ax[1].hist(df_rsa_ent.RSA[idx_jmis], bins=20, alpha=0.6, density=true, label="Background mismatches")
ax[1].hist(df_rsa_ent[df_rsa_ent.entropy .>= 0.01, :RSA], bins=20, alpha=0.5, density=true, label="Variable sites")
ax[1].legend(fontsize=12)
ax[1].set_xlabel("Relative solvent accessibility", fontsize=12)
ax[1].set_ylabel("Density", fontsize=12)

# Coupling Vs RSA scatter plot
ax[2].scatter(j_rsa_max.RSA, j_rsa_max.J, s=12, alpha=0.6)
ax[2].set_xlabel("Relative solvent accessibility", fontsize=12)
ax[2].set_ylabel("Root Mean Square of interactions", fontsize=12)

fig.savefig("results/figures/si/fig_s11.pdf")
close(fig)