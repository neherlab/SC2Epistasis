""" Shared plotting function: z-score bands with domain annotations """

using PyCall
@isdefined(mpatches) || (const mpatches = pyimport("matplotlib.patches"))

function plot_raw_zscores_bands(cp_list::Vector{Tuple{S1,S1}},
    z_dict::Dict{Tuple{S2,S2},Vector{Float64}},
    s_dict::Dict{Tuple{S2,S2},Vector{Int}},
    cdiff::DataFrame;
    x_space::Int=100,
    doms_label::Vector{S1}=["NTD", "RBD", "CTD1", "CTD2", "FP", "HR1", "CH", "CD", "HR2"],
    doms_edges::Vector{UnitRange{Int}}=[14:305, 319:541, 542:590, 591:690, 788:806, 910:984, 985:1034, 1035:1067, 1163:1212],
    doms_col::Vector{S1}=["cyan", "blue", "orange", "yellow", "red", "green", "lightgreen", "hotpink", "violet"],
    box_height::Float64=0.12,
    rasterized=true,
    fig=nothing, axs=nothing) where {S1<:AbstractString,S2<:AbstractString}

    standalone = isnothing(fig)
    n = length(cp_list)
    L = maximum([maximum(s_dict[cp]) for cp in cp_list])  # protein length

    domains = [(doms_edges[n][1], doms_edges[n][end], doms_label[n], doms_col[n]) for n in eachindex(doms_label)]

    # --- compute global max across all clade pairs ---
    global_max = maximum([maximum(z_dict[cp]) for cp in cp_list])
    y_max = global_max * 1.05   # add 5% margin

    if standalone
        fig, axs = subplots(n, 1, figsize=(12, 3.5 * n), sharex=true)
        if n == 1
            axs = [axs]
        end
    end

    for (i, cp_plt) in enumerate(cp_list)
        ax = axs[i]

        # Add domain bands to each subplot
        for (start, stop, label, color) in domains
            ax.axvspan(start, stop, alpha=0.20, color=color, zorder=0)
        end

        site_diff = cdiff[(cdiff.clade1.==cp_plt[1]).&(cdiff.clade2.==cp_plt[2]), :site]
        sd_str = string.(site_diff)
        sd_str[1:end-1] .*= ", "
        if length(sd_str) >= 30
            insert!(sd_str, div(length(sd_str) + 5, 2), "\n    ")
        end

        ax.plot(s_dict[cp_plt], z_dict[cp_plt], ".-", zorder=2, rasterized=rasterized)
        ax.vlines(site_diff, 0.0, global_max, color="orange", zorder=3)
        # ax.set_ylabel(cp_plt[1] * "-" * cp_plt[2], fontsize=14, labelpad=10, rotation=90, ha="left", va="center")
        ax.yaxis.set_label_position("right")  # moves labels to the right side of the subplots
        ax.text(0.98, 0.75, cp_plt[1] * "-" * cp_plt[2], transform=ax.transAxes,
            fontsize=16, fontweight="bold", ha="right", va="top")

        # set common y-axis limit
        ax.set_ylim(-0.2, y_max)

        # only clear tick labels for non-bottom subplots
        if i != n
            ax.tick_params(labelbottom=false)
        end

        legend_str = "Positions: " * join(sd_str)
        ax.legend([legend_str], loc="best", markerscale=0, handlelength=0, fontsize=10)
    end

    # common y-axis label (standalone only; composite caller places it)
    if standalone
        fig.text(0.08, 0.5, "Average fitness z-score", va="center", rotation="vertical", fontsize=16)
    end

    # bottom x-axis
    axs[end].set_xlabel("Residue", fontsize=16, labelpad=60)  # push further down
    axs[end].set_xticks(collect(0:x_space:L))
    axs[end].xaxis.set_tick_params(labelsize=12)

    # --- domain boxes positioning ---
    y_bottom = -(box_height + 0.0025)  # position so top edge touches the axis

    # --- grey background box ---
    bg = mpatches.Rectangle(
        (1, y_bottom),
        L - 1,
        box_height,
        facecolor="lightgrey",
        edgecolor="none",
        transform=axs[end].get_xaxis_transform(),
        clip_on=false,
        zorder=0
    )
    axs[end].add_patch(bg)

    # --- colored domain boxes ---
    for (start, stop, label, color) in domains
        rect = mpatches.FancyBboxPatch(
            (start, y_bottom),
            stop - start,
            box_height,
            boxstyle="round4,pad=0.0",
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            transform=axs[end].get_xaxis_transform(),
            clip_on=false,
            zorder=1
        )
        axs[end].add_patch(rect)

        # domain label INSIDE the box
        axs[end].text(
            (start + stop) / 2,
            y_bottom + box_height / 2,
            label,
            ha="center", va="center",
            transform=axs[end].get_xaxis_transform(),
            fontsize=10, fontweight="bold",
            color="white",
            path_effects=[pyimport("matplotlib.patheffects").withStroke(linewidth=2, foreground="black")]
        )
    end

    # --- create secondary invisible x-axis for shifted tick labels ---
    ax_ticks = axs[end].twiny()

    # make it share the same x-limits
    ax_ticks.set_xlim(axs[end].get_xlim())

    # explicit ticks at 0, 100, 200, ... L
    xticks = collect(0:x_space:L)
    ax_ticks.set_xticks(xticks)
    ax_ticks.set_xticklabels([string(t) for t in xticks])

    # position this new axis at the bottom edge of the boxes
    ax_ticks.spines["top"].set_position(("axes", y_bottom - 0.015))

    # hide all spines including the top one
    for spine in ["top", "bottom", "left", "right"]
        ax_ticks.spines[spine].set_visible(false)
    end

    # show ticks on top (which is positioned at bottom of boxes) and add padding for labels
    ax_ticks.tick_params(axis="x", which="both", top=true, bottom=false, pad=-21, labelsize=14)

    # --- hide ticks and labels from the main x-axis (keep the spine line) ---
    axs[end].tick_params(axis="x", which="both", bottom=false, labelbottom=false)

    # --- main axis label further below everything ---
    axs[end].set_xlabel("Residue", fontsize=16, labelpad=55)

    if standalone
        fig.subplots_adjust(bottom=0.35, left=0.12, hspace=0.3)
        fig.tight_layout(rect=[0.06, 0, 1, 1])
    end

    return fig, axs

end
