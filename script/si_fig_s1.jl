""" Composite SI figure S1: panels A (z-score bands) and B (sphere fraction) """

# Import packages
using DataFrames, CSV, JLD2, PyPlot, PyCall, PdbTool
include("z_plot_raw_zscores_bands.jl")
include("z_plot_sphere_frac.jl")

###############################---------------------------------

# Load data for panel A
@load "results/fit_zscore_epi.jld2"
z_dict = fit_epi.z_score
s_dict = fit_epi.sites

clade_diff = CSV.read("results/clade_diff.csv", DataFrame)
clade_diff_S = clade_diff[clade_diff.prot.=="S", :]

# Load data for panel B
pdbs, af_pdb = SC2Epistasis.read_pdbs("data/ref_seq/Spike.txt")
site_shannon_ent = CSV.read("data/nextstrain_staging_nextclade_sars-cov-2_diversity.tsv", DataFrame)
ent_thr = 0.01
var_sites = site_shannon_ent.position[site_shannon_ent.entropy.>=ent_thr]

# Clade pairs
cp_plt = [("20I", "21J"), ("20I", "21L"), ("21J", "21K"), ("21K", "23I")]
n = length(cp_plt)

# Panel B parameters
radii = [5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
z_thr = [[1.5, 1.7, 2.0], [2.0, 2.25, 2.5], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]]

frac_vec = [zeros(Float64, length(radii), length(z_thr[i])) for i in eachindex(cp_plt)]
rnd_frac  = [zeros(Float64, length(radii)) for i in eachindex(cp_plt)]
for i in eachindex(cp_plt)
    frac_vec[i] = SC2Epistasis.frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, pdbs, af_pdb, radii, z_thr[i])
    rnd_frac[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, pdbs, af_pdb, radii, var_sites; z_thr=z_thr[i][1], nsamp=100)
end

###############################---------------------------------

# Build composite figure: n rows × 2 cols, width ratio 12:4.5
fig = figure(figsize=(16.5, 3.5 * n))
gs = fig.add_gridspec(n, 2, width_ratios=[12, 4.5], hspace=0.3, wspace=0.15)

axs_A = Vector{Any}(undef, n)
axs_B = Vector{Any}(undef, n)
axs_A[1] = fig.add_subplot(gs[1, 1])
axs_B[1] = fig.add_subplot(gs[1, 2])
for i in 2:n
    axs_A[i] = fig.add_subplot(gs[i, 1], sharex=axs_A[1])
    axs_B[i] = fig.add_subplot(gs[i, 2], sharex=axs_B[1])
end

# Plot panels
plot_raw_zscores_bands(cp_plt, z_dict, s_dict, clade_diff_S; rasterized=false, fig=fig, axs=axs_A)
plot_sphere_frac(cp_plt, frac_vec, rnd_frac; radii=radii, z_thr=z_thr, fig=fig, axs=axs_B)

# Panel labels
axs_A[1].text(-0.04, 1.02, "A", transform=axs_A[1].transAxes, fontsize=18, fontweight="bold", va="bottom")
axs_B[1].text(-0.08, 1.02, "B", transform=axs_B[1].transAxes, fontsize=18, fontweight="bold", va="bottom")

# Common y-axis labels
fig.text(0.01, 0.5, "Average fitness z-score", va="center", rotation="vertical", fontsize=16)
fig.text(0.99, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < d_{thr}",
    va="center", ha="right", rotation="vertical", fontsize=16)

fig.savefig("results/figures/si/fig_s1.pdf", bbox_inches="tight")
close(fig)
