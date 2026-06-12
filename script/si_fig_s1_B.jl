""" Script to generate panel B of figure 1 of the SI """

# Import packages
using DataFrames, CSV, JLD2, PyPlot, PdbTool
include("z_plot_sphere_frac.jl")

###############################---------------------------------

# Script body

# Load PDB files
pdbs, af_pdb = SC2Epistasis.read_pdbs("data/ref_seq/Spike.txt")

# Load clade differences
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)
clade_diff = clade_diff[clade_diff.prot.=="S", :]

# Load raw epistatic signal
@load "results/fit_zscore_epi.jld2"
z_dict = fit_epi.z_score
s_dict = fit_epi.sites

# Clade pairs to be plotted
cp_plt = [("20I", "21J"), ("20I", "21L"), ("21J", "21K"), ("21K", "23I")]

# Define set of radii and z-score thresholds
radii = [5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
z_thr = [[1.5, 1.7, 2.0], [2.0, 2.25, 2.5], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]]

# Import dataframe with sites and Shannon entropy
site_shannon_ent = CSV.read("data/nextstrain_staging_nextclade_sars-cov-2_diversity.tsv", DataFrame)
ent_thr = 0.01
var_sites = site_shannon_ent.position[site_shannon_ent.entropy.>=ent_thr]

# Initialize arrays to store fractions
frac_vec = [zeros(Float64, length(radii), length(z_thr[i])) for i in 1:length(cp_plt)]
rnd_frac = [zeros(Float64, length(radii)) for i in 1:length(cp_plt)]

for i in eachindex(cp_plt)
    # Compute fractions in spheres
    frac_vec[i] = SC2Epistasis.frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff, pdbs, af_pdb, radii, z_thr[i])
    # Compute random benchmark
    rnd_frac[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff, pdbs, af_pdb, radii, var_sites; z_thr=z_thr[i][1], nsamp=100)
end

fig, ax = plot_sphere_frac(cp_plt, frac_vec, rnd_frac; radii=radii, z_thr=z_thr)
fig.savefig("results/figures/si/fig_s1_B.pdf")
close(fig)

# --- Fraction vs z-score for two fixed radii ---
radii_sel = [8.0, 15.0]
z_thr_range = collect(0.5:0.1:2.0)

frac_vec_z = [zeros(Float64, length(radii_sel), length(z_thr_range)) for i in eachindex(cp_plt)]
rnd_frac_z = [zeros(Float64, length(radii_sel)) for i in eachindex(cp_plt)]
for i in eachindex(cp_plt)
    frac_vec_z[i] = SC2Epistasis.frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff, pdbs, af_pdb, radii_sel, z_thr_range)
    rnd_frac_z[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff, pdbs, af_pdb, radii_sel, var_sites; z_thr=z_thr_range[1], nsamp=100)
end

fig2, ax2 = plot_sphere_frac_vs_z(cp_plt, frac_vec_z, z_thr_range, rnd_frac_z; radii_sel=radii_sel)
fig2.savefig("results/figures/si/fig_s1_B_vs_z.pdf")
close(fig2)