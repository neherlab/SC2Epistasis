""" Script to generate panel B of figure 7 """

# Import packages
using DataFrames, CSV, PyPlot

# Plotting function
function plot_model_pred(av_z::Vector{Float64}, av_dz::Vector{Float64}; nrank::Int=20)

    # Isolating well predicted mutations
    idx_pred = findall(av_dz .< av_z)
    av_z_pred = av_z[idx_pred]
    av_dz_pred = av_dz[idx_pred]

    idx_sort = sortperm((av_z_pred .^ 2 .- av_dz_pred .^ 2), rev=true)
    label = ["ϵ²", "z²"]

    fig, ax = subplots(figsize=(5, 6))
    for n in 1:nrank
        ax.plot([0, 1], vcat(av_dz_pred[idx_sort[n]]^2, av_z_pred[idx_sort[n]]^2), "o-", alpha=0.7)
        if n == 1
            ax.set_xticks([0, 1])
            ax.set_xticklabels(label, fontsize=15)
        end
    end
    ax.set_xlim(-0.1, 1.1)
    fig.tight_layout()

    return fig, ax

end

# Read raw and model curated z-scores for test mutations
z_df = CSV.read("results/model_test_zscore.csv", DataFrame)
av_z = z_df.av_z
av_dz = z_df.av_dz

fig, ax = plot_model_pred(av_z, av_dz)
fig.savefig("results/figures/fig_coup_valid_B.pdf")
close(fig)

