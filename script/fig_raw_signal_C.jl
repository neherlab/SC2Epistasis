""" Script to generate panel C of figure 3 """

# Import packages
using DataFrames, CSV, JLD2, PyPlot

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
        frac = SC2Epistasis.frac_in_sphere(cpair, z_dict, s_dict, clade_diff_prot, pdbs, af_pdb, radii, z_thr)
        rnd_frac, _ = SC2Epistasis.rand_frac_in_sphere(cpair, z_dict, s_dict, clade_diff_prot, pdbs, af_pdb, radii, z_thr=z_thr[1], nsamp=nsamp)
        for k in eachindex(z_thr[i])
            ax.plot(radii, frac[:, k], ".-")
        end
        ax.plot(radii, rnd_frac, ".-")

        # Subplot title
        ax.set_title(cpair[1] * "-" * cpair[2], fontsize=14)

        # set common y-axis limit
        ax.set_ylim(0.0, 1.25)

        # Subplot legend
        ax.legend(vcat("z ≥ " .* string.(z_thr[i]), "Random"), fontsize=14)

        # only clear tick labels for non-bottom subplots
        if i != n
            ax.tick_params(labelbottom=false)
        end

    end

    # common y-axis label
    fig.text(0.04, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < d_{thr}", va="center", rotation="vertical", fontsize=16)

    # bottom x-axis
    axs[end].set_xlabel("Sphere radius (Å)", fontsize=16)

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

fig, ax = plot_sphere_frac(cp_plt, z_dict, s_dict, clade_diff, pdbs, af_pdb)