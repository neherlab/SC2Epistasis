"""
Master script: composite figures (A = z-score bands, B = sphere fraction) for
main figure and SI figures S1 and S2. Outputs both the vs-radius and vs-z-score
versions for each figure set.
"""

# Import packages
using DataFrames, CSV, JLD2, PyPlot, PyCall, PdbTool, DelimitedFiles
include("z_plot_raw_zscores_bands.jl")
include("z_plot_sphere_frac.jl")

###############################---------------------------------

# Load all shared data once

println("Loading z-score data...")
@load "results/fit_zscore_epi.jld2"
z_dict = fit_epi.z_score
s_dict = fit_epi.sites

println("Loading clade differences...")
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)
clade_diff_S = clade_diff[clade_diff.prot .== "S", :]

println("Loading PDB structures...")
pdbs, af_pdb = SC2Epistasis.read_pdbs("data/ref_seq/Spike.txt")

println("Loading distance matrix...")
dist_mat = DelimitedFiles.readdlm("results/dist_mat.txt", Float64)

println("Loading diversity data...")
site_shannon_ent = CSV.read("data/nextstrain_staging_nextclade_sars-cov-2_diversity.tsv", DataFrame)
var_sites = site_shannon_ent.position[site_shannon_ent.entropy .>= 0.01]
L = length(af_pdb.chain["A"].residue) # number of residues in the Spike protein
all_sites = collect(1:L)

radii = [5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
radii_sel = [5.0, 10.0]
z_thr_range = collect(0.5:0.1:2.0)
nsamples = 100

###############################---------------------------------

# Helper: build one composite figure (A left, B right) and return it with axs_B ready to plot
function _make_composite(cp_plt, z_dict, s_dict, clade_diff_S)
    n = length(cp_plt)
    fig = figure(figsize=(16.5, 3.5 * n))
    gs = fig.add_gridspec(n, 2, width_ratios=[12, 4.5], hspace=0.1, wspace=0.05)
    axs_A = Vector{Any}(undef, n)
    axs_B = Vector{Any}(undef, n)
    axs_A[1] = fig.add_subplot(gs[1, 1])
    axs_B[1] = fig.add_subplot(gs[1, 2])
    for i in 2:n
        axs_A[i] = fig.add_subplot(gs[i, 1], sharex=axs_A[1])
        axs_B[i] = fig.add_subplot(gs[i, 2], sharex=axs_B[1])
    end
    plot_raw_zscores_bands(cp_plt, z_dict, s_dict, clade_diff_S; rasterized=false, fig=fig, axs=axs_A)
    axs_A[1].text(-0.04, 1.02, "A", transform=axs_A[1].transAxes, fontsize=18, fontweight="bold", va="bottom")
    axs_B[1].text(-0.08, 1.02, "B", transform=axs_B[1].transAxes, fontsize=18, fontweight="bold", va="bottom")
    fig.text(0.08, 0.5, "Average fitness z-score", va="center", rotation="vertical", fontsize=16)
    fig.text(0.95, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < r",
        va="center", ha="right", rotation="vertical", fontsize=16)
    return fig, axs_B
end

# Helper: compute B-panel data for a given set of clade pairs and z-score thresholds
function _compute_frac(cp_plt, z_thr, z_dict, s_dict, clade_diff_S, all_sites, var_sites,
    dist_mat, radii, radii_sel, z_thr_range)
    frac_vec = [zeros(Float64, length(radii), length(z_thr[i])) for i in eachindex(cp_plt)]
    rnd_frac = [zeros(Float64, length(radii), 1) for i in eachindex(cp_plt)]
    frac_vec_z = [zeros(Float64, length(radii_sel), length(z_thr_range)) for i in eachindex(cp_plt)]
    rnd_frac_z = [zeros(Float64, length(radii_sel), length(z_thr_range)) for i in eachindex(cp_plt)]
    rnd_frac_var = [zeros(Float64, length(radii), 1) for i in eachindex(cp_plt)]
    rnd_frac_z_var = [zeros(Float64, length(radii_sel), length(z_thr_range)) for i in eachindex(cp_plt)]
    for i in eachindex(cp_plt)
        println("  $(cp_plt[i][1])-$(cp_plt[i][2]): computing sphere fractions...")
        frac_vec[i] = SC2Epistasis.frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, dist_mat, radii, z_thr[i])
        rnd_frac[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, dist_mat, radii, [z_thr[i][1]]; var_sites=all_sites, nsamp=nsamples)
        frac_vec_z[i] = SC2Epistasis.frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, dist_mat, radii_sel, z_thr_range)
        rnd_frac_z[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, dist_mat, radii_sel, z_thr_range; var_sites=all_sites, nsamp=nsamples)
        rnd_frac_var[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, dist_mat, radii, [z_thr[i][1]]; var_sites=var_sites, nsamp=nsamples)
        rnd_frac_z_var[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff_S, dist_mat, radii_sel, z_thr_range; var_sites=var_sites, nsamp=nsamples)
    end
    return frac_vec, rnd_frac, frac_vec_z, rnd_frac_z, rnd_frac_var, rnd_frac_z_var
end

###############################---------------------------------

# Figure configurations: (clade pairs, z-score thresholds, output path prefix)
configs = [
    (
        [("21J", "21L"), ("21K", "21L"), ("21L", "22E"), ("21L", "23I")],
        [[1.75, 2.0, 2.25], [1.25, 1.5, 1.75], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]],
        "results/figures/fig_raw_signal"
    ),
    (
        [("20I", "21J"), ("20I", "21L"), ("21J", "21K"), ("21K", "23I")],
        [[1.5, 1.7, 2.0], [2.0, 2.25, 2.5], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]],
        "results/figures/si/fig_s1"
    ),
    (
        [("21K", "22B"), ("21K", "23A"), ("21L", "22B"), ("21L", "23A")],
        [[1.5, 1.75, 2.0], [2.0, 2.25, 2.5], [1.25, 1.5, 1.75], [2.0, 2.25, 2.5]],
        "results/figures/si/fig_s2"
    ),
]

# Make figures: vs radius and vs z-score for each configuration. Including both random tests
for (k, (cp_plt, z_thr, outpath)) in enumerate(configs)

    println("\n[$k/$(length(configs))] Computing data for $(basename(outpath))...")
    frac_vec, rnd_frac, frac_vec_z, rnd_frac_z, rnd_frac_var, rnd_frac_z_var = _compute_frac(cp_plt, z_thr, z_dict, s_dict, clade_diff_S, all_sites, var_sites, dist_mat, radii, radii_sel, z_thr_range)

    # Version 1: B = fraction vs radius
    println("  Rendering $(basename(outpath)).pdf...")
    fig, axs_B = _make_composite(cp_plt, z_dict, s_dict, clade_diff_S)
    plot_sphere_frac(cp_plt, frac_vec, rnd_frac, rnd_frac_var; radii=radii, z_thr=z_thr, fig=fig, axs=axs_B)
    fig.tight_layout(rect=[0.06, 0, 0.96, 1])
    fig.savefig(outpath * ".pdf", bbox_inches="tight")
    close(fig)
    println("  Saved $(outpath).pdf")

    # Version 2: B = fraction vs z-score threshold
    println("  Rendering $(basename(outpath))_vs_z.pdf...")
    fig, axs_B = _make_composite(cp_plt, z_dict, s_dict, clade_diff_S)
    plot_sphere_frac_vs_z(cp_plt, frac_vec_z, z_thr_range, rnd_frac_z, rnd_frac_z_var; radii_sel=radii_sel, fig=fig, axs=axs_B)
    fig.tight_layout(rect=[0.06, 0, 0.96, 1])
    fig.savefig(outpath * "_vs_z.pdf", bbox_inches="tight")
    close(fig)
    println("  Saved $(outpath)_vs_z.pdf")

end

println("\nDone.")
