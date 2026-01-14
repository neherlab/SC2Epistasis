""" Script to generate panel B of figure 4 """

# Import packages
using DataFrames, CSV, PyPlot

# Plotting function

function plot_model_fit(fit_df_vec::Vector{DataFrame}, e0::Float64;
    labels::Vector{String}=["Uniform", "3d", "3d-sigmoid"],
    linestyles::Vector{String}=["-", "--", "-."],
    shades1::Vector{NTuple{4,Float64}}=[get_cmap("Blues")(0.4 + 0.2 * i) for i in 0:2],
    shades2::Vector{NTuple{4,Float64}}=[get_cmap("Reds")(0.4 + 0.2 * i) for i in 0:2])

    @assert length(fit_df_vec) == length(labels) == length(linestyles) == length(shades1) == length(shades2)

    fig, ax = subplots(2, 1, figsize=(5, 10))

    for n in eachindex(fit_df_vec)
        fit_df = fit_df_vec[n]
        ax[1].plot(fit_df.l1, fit_df.energy, marker=".", markersize=5.0,
            linestyle=linestyles[n], color=shades1[n], label=labels[n])
    end

    ax[1].hlines(e0, fit_df_vec[end].l1[1], fit_df_vec[end].l1[end], color="black", linestyle="-", label="Estimated noise")
    ax[1].set_xscale("log")
    ax[1].set_ylabel("Energy", fontsize=14)
    ax[1].legend(loc="best", fontsize=12)
    ax[1].tick_params(labelbottom=false)
    for n in eachindex(fit_df_vec)
        fit_df = fit_df_vec[n]
        ax[2].plot(fit_df.l1, 1 .- fit_df.tot_zf, marker=".", markersize=5.0,
            linestyle=linestyles[n], color=shades2[n], label=labels[n])
    end
    ax[2].set_xscale("log")
    ax[2].set_yscale("log")
    ax[2].set_xlabel("λ₁", fontsize=16)
    ax[2].set_ylabel("Non-zero fraction", fontsize=14)
    ax[2].legend(loc="best", fontsize=12)
    fig.tight_layout()

    return fig, ax

end

##############################--------------------------------

# Script body

# Read dataframes with model fits
path_model_fit = "results/model_fit/" # path to folder with model fit results
files = readdir(path_model_fit, join=true) # list of files in the folder
fit_df_vec = [DataFrame(CSV.File(f)) for f in files]; # vector of dataframes

# Compute baseline noise level
include("usa_uk_ener.jl") # e0 is defined here

# Plot model fits
fig, ax = plot_model_fit(fit_df_vec[[3, 1, 2]], e0); # re-order to match legend
fig.savefig("results/figures/fig_model_fit_B.pdf", dpi=500)
close(fig)
