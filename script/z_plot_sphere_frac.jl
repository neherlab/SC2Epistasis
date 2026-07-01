""" Shared plotting function: fraction of positions in sphere vs radius """

function plot_sphere_frac(cp_plt::Vector{Tuple{S1,S1}}, frac::Vector{Matrix{Float64}}, rnd_frac::Vector{Matrix{Float64}}, rnd_frac_var::Vector{Matrix{Float64}};
    radii::Vector{Float64}=[5.0, 8.0, 10.0, 12.0, 15.0, 20.0],
    z_thr::Vector{Vector{Float64}}=[[1.75, 2.0, 2.25], [1.25, 1.5, 1.75], [1.5, 1.75, 2.0], [1.5, 1.75, 2.0]],
    fig=nothing, axs=nothing) where {S1<:AbstractString}

    standalone = isnothing(fig)
    n = length(cp_plt)  # number of clade pairs to be plotted

    if standalone
        fig, axs = subplots(n, 1, figsize=(4.5, 3.5 * n), sharex=true)
        if n == 1
            axs = [axs]
        end
    end

    for (i, cpair) in enumerate(cp_plt)
        ax = axs[i]

        # Plot fraction in sphere and random benchmark
        _markers = ["o", "s", "^", "D"]
        _lstyles = ["-", "--", ":", "-."]
        for k in eachindex(z_thr[i])
            ax.plot(radii, frac[i][:, k], marker=_markers[k], linestyle=_lstyles[k], label="z ≥ $(z_thr[i][k])")
        end
        ax.plot(radii, rnd_frac[i], marker="x", linestyle="-.", label="Random")
        ax.plot(radii, rnd_frac_var[i], marker="P", linestyle=(0, (3, 1, 1, 1)), label="Random var")

        # Move y-tick labels to the right and increase fontsize
        ax.yaxis.tick_right()
        ax.tick_params(axis="y", labelsize=11)

        # set common y-axis limit
        ax.set_ylim(0.0, 1.05)

        # Subplot legend
        ax.legend(fontsize=13)

        # only clear tick labels for non-bottom subplots
        if i != n
            ax.tick_params(labelbottom=false)
        end

    end

    if standalone
        fig.tight_layout(rect=(0.02, 0.025, 0.87, 0.98))
        fig.text(0.92, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < r", va="center", rotation="vertical", fontsize=16)
    end

    # bottom x-axis
    axs[end].set_xlabel("Sphere radius (Å)", fontsize=16)
    axs[end].tick_params(axis="x", labelsize=11)

    return fig, axs

end

function plot_sphere_frac_vs_z(cp_plt::Vector{Tuple{S1,S1}}, frac_z::Vector{Matrix{Float64}},
    z_thr_range::Vector{Float64}, rnd_frac::Vector{Matrix{Float64}}, rnd_frac_var::Vector{Matrix{Float64}};
    radii_sel::Vector{Float64}=[8.0, 15.0],
    fig=nothing, axs=nothing) where {S1<:AbstractString}

    standalone = isnothing(fig)
    n = length(cp_plt)

    if standalone
        fig, axs = subplots(n, 1, figsize=(4.5, 3.5 * n), sharex=true)
        if n == 1
            axs = [axs]
        end
    end

    for (i, cpair) in enumerate(cp_plt)
        ax = axs[i]

        _markers = ["o", "s"]
        _lstyles = ["-", "--"]
        for (k, r) in enumerate(radii_sel)
            frac_label = i == 1 ? "r = $(Int(r)) Å" : "_nolegend_"
            ax.plot(z_thr_range, frac_z[i][k, :], marker=_markers[k], linestyle=_lstyles[k], label=frac_label, color="C$(k-1)")
            ax.plot(z_thr_range, rnd_frac[i][k, :], linestyle=":", color="C$(k-1)", alpha=0.6)
            ax.plot(z_thr_range, rnd_frac_var[i][k, :], linestyle="-.", color="C$(k-1)", alpha=0.6)
        end

        ax.yaxis.tick_right()
        ax.tick_params(axis="y", labelsize=11)
        ax.set_ylim(0.0, 1.05)
        if i == 1
            ax.legend(fontsize=13)
        elseif i == 2
            ax.plot(Float64[], Float64[], linestyle=":", color="k", label="Random")
            ax.plot(Float64[], Float64[], linestyle="-.", color="k", label="Random var.")
            ax.legend(fontsize=13)
        end

        if i != n
            ax.tick_params(labelbottom=false)
        end
    end

    if standalone
        fig.tight_layout(rect=(0.02, 0.025, 0.87, 0.98))
        fig.text(0.92, 0.5, "Fraction of " * L"i\;" * "s.t. " * L"\exists j:\:d(i,j) < r", va="center", rotation="vertical", fontsize=16)
    end
    axs[end].set_xlabel("z-score threshold", fontsize=16)
    axs[end].tick_params(axis="x", labelsize=11)

    return fig, axs

end
