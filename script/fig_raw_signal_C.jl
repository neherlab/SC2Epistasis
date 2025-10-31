""" Script to generate panel C of figure 3 """

# Import packages
using DataFrames, CSV, JLD2, PyPlot, PdbTool

# Plotting function
function plot_sphere_frac(cp_plt::Vector{Tuple{S1,S1}}, z_dict::Dict{Tuple{S2,S2},Vector{Float64}}, s_dict::Dict{Tuple{S2,S2},Vector{Int}},
    clade_diff_prot::DataFrame, pdbs::Vector{PdbTool.Pdb}, af_pdb::PdbTool.Pdb;
    radii::Vector{Float64}=[5.0, 8.0, 10.0, 12.0, 15.0, 20.0],
    z_thr::Vector{Vector{Float64}}=[[1.75, 2.0, 2.25], [1.25, 1.5, 1.75], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]],
    nsamp::Int=100) where {S1<:AbstractString,S2<:AbstractString}

    n = length(cp_plt)  # number of clade pairs to be plotted

    fig, axs = subplots(n, 1, figsize=(7, 5.25 * n), sharex=true)
    if n == 1
        axs = [axs]
    end

    for (i, cpair) in enumerate(cp_plt)
        ax = axs[i]

        # Compute and plot fraction in sphere and random benchmark
        frac = SC2Epistasis.frac_in_sphere(cpair, z_dict, s_dict, clade_diff_prot, pdbs, af_pdb, radii, z_thr[i])
        rnd_frac, _ = SC2Epistasis.rand_frac_in_sphere(cpair, z_dict, s_dict, clade_diff_prot, pdbs, af_pdb, radii, z_thr=z_thr[i][1], nsamp=nsamp)
        for k in eachindex(z_thr[i])
            ax.plot(radii, frac[:, k], ".-")
        end
        ax.plot(radii, rnd_frac, ".-")

        # Subplot clade pair as y-axis label
        ax.set_ylabel(cpair[1] * "-" * cpair[2], fontsize=14, labelpad=10, rotation=90, ha="left", va="center")
        ax.yaxis.set_label_position("right")  # moves labels to the right side of the subplots

        # set common y-axis limit
        ax.set_ylim(0.0, 1.05)

        # Subplot legend
        ax.legend(vcat("z ≥ " .* string.(z_thr[i]), "Random"), fontsize=14)

        # only clear tick labels for non-bottom subplots
        if i != n
            ax.tick_params(labelbottom=false)
        end

    end

    # Then adjust margins to accommodate the common y-axis label and right-side labels
    fig.tight_layout(rect=(0.06, 0.025, 0.98, 0.98))

    # common y-axis label
    fig.text(0.01, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < d_{thr}", va="center", rotation="vertical", fontsize=14)

    # bottom x-axis
    axs[end].set_xlabel("Sphere radius (Å)", fontsize=14)

    return fig, axs

end

function plot_sphere_frac(cp_plt::Vector{Tuple{S1,S1}}, frac::Vector{Matrix{Float64}}, rnd_frac::Vector{Vector{Float64}};
    radii::Vector{Float64}=[5.0, 8.0, 10.0, 12.0, 15.0, 20.0],
    z_thr::Vector{Vector{Float64}}=[[1.75, 2.0, 2.25], [1.25, 1.5, 1.75], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]]) where {S1<:AbstractString}

    n = length(cp_plt)  # number of clade pairs to be plotted

    fig, axs = subplots(n, 1, figsize=(5, 6 * n), sharex=true)
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
        ax.set_ylabel(cpair[1] * "-" * cpair[2], fontsize=14, labelpad=10, rotation=90, ha="left", va="center")
        ax.yaxis.set_label_position("right")  # moves labels to the right side of the subplots

        # set common y-axis limit
        ax.set_ylim(0.0, 1.05)

        # Subplot legend
        ax.legend(vcat("z ≥ " .* string.(z_thr[i]), "Random"), fontsize=14)

        # only clear tick labels for non-bottom subplots
        if i != n
            ax.tick_params(labelbottom=false)
        end

    end

    # Then adjust margins to accommodate the common y-axis label and right-side labels
    fig.tight_layout(rect=(0.06, 0.025, 0.98, 0.98))

    # common y-axis label
    fig.text(0.01, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < d_{thr}", va="center", rotation="vertical", fontsize=14)

    # bottom x-axis
    axs[end].set_xlabel("Sphere radius (Å)", fontsize=14)

    return fig, axs

end

###############################---------------------------------

# Script body

# Load PDB files
pdbs, af_pdb = SC2Epistasis.read_pdbs()

# Load clade differences
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)
clade_diff = clade_diff[clade_diff.prot.=="S", :]  # only keep spike protein differences

# Load raw epistatic signal
@load "results/fit_zscore_epi.jld2"
z_dict = fit_epi.z_score
s_dict = fit_epi.sites

# Clade pairs to be plotted
cp_plt = [("21J", "21L"), ("21K", "21L"), ("21L", "22E"), ("21L", "23I")]

# Define set of radii and z-score thresholds
radii = [5.0, 8.0, 10.0, 12.0, 15.0, 20.0]
z_thr = [[1.75, 2.0, 2.25], [1.25, 1.5, 1.75], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]]

# Initialize arrays to store fractions
frac_vec = [zeros(Float64, length(radii), length(z_thr[i])) for i in 1:length(cp_plt)]
rnd_frac = [zeros(Float64, length(radii)) for i in 1:length(cp_plt)]

for i in eachindex(cp_plt)
    # Compute fractions in spheres
    frac_vec[i] = SC2Epistasis.frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff, pdbs, af_pdb, radii, z_thr[i])
    # Compute random benchmark
    rnd_frac[i], _ = SC2Epistasis.rand_frac_in_sphere(cp_plt[i], z_dict, s_dict, clade_diff, pdbs, af_pdb, radii, z_thr=z_thr[i][1], nsamp=100)
end

fig, ax = plot_sphere_frac(cp_plt, z_dict, s_dict, clade_diff, pdbs, af_pdb)