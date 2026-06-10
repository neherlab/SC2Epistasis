""" Script to generate panel B of figure 2 of the SI """

# Import packages
using DataFrames, CSV, JLD2, PyPlot, PdbTool

# Plotting function
function plot_sphere_frac(cp_plt::Vector{Tuple{S1,S1}}, frac::Vector{Matrix{Float64}}, rnd_frac::Vector{Vector{Float64}};
    radii::Vector{Float64}=[5.0, 8.0, 10.0, 12.0, 15.0, 20.0],
    z_thr::Vector{Vector{Float64}}=[[1.75, 2.0, 2.25], [1.25, 1.5, 1.75], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]]) where {S1<:AbstractString}

    n = length(cp_plt)  # number of clade pairs to be plotted

    fig, axs = subplots(n, 1, figsize=(4.5, 5.4 * n), sharex=true, layout="constrained")
    if n == 1
        axs = [axs]
    end

    for (i, cpair) in enumerate(cp_plt)
        ax = axs[i]

        # Plot fraction in sphere and random benchmark
        for k in eachindex(z_thr[i])
            ax.plot(radii, frac[i][:, k], ".-")
        end
        ax.plot(radii, rnd_frac[i], ".-")

        # Subplot clade pair as y-axis label
        ax.set_ylabel(cpair[1] * "-" * cpair[2], fontsize=14, labelpad=18, rotation=90, ha="left", va="center")

        # Move y-tick labels to the right and increase fontsize
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", labelsize=11)

        # set common y-axis limit
        ax.set_ylim(0.0, 1.05)

        # Subplot legend
        ax.legend(vcat("z ≥ " .* string.(z_thr[i]), "Random"), fontsize=13)

        # only clear tick labels for non-bottom subplots
        if i != n
            ax.tick_params(labelbottom=false)
        end

    end

    # Then adjust margins to accommodate the common y-axis label and right-side labels
    fig.tight_layout(rect=(0.02, 0.025, 0.87, 0.98))

    # common y-axis label on the right side
    fig.text(0.92, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < d_{thr}", va="center", rotation="vertical", fontsize=16)

    # bottom x-axis
    axs[end].set_xlabel("Sphere radius (Å)", fontsize=16)
    axs[end].tick_params(axis="x", labelsize=11)

    return fig, axs

end

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
cp_plt = [("21K", "22B"), ("21K", "23A"), ("21L", "22B"), ("21L", "23A")]

# Define set of radii and z-score thresholds
radii = [5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
z_thr = [[1.5, 1.75, 2.0], [2.0, 2.25, 2.5], [1.25, 1.5, 1.75], [2.0, 2.25, 2.5]]

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
fig.savefig("results/figures/si/fig_s2_B.pdf")
close(fig)