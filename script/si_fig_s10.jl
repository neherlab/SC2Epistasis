""" Script to generate figure S10 of the SI"""

# Import packages
using DataFrames, CSV, PyPlot, PyCall
mpatches = pyimport("matplotlib.patches") # import from matplotlib to draw boxes

# Plotting function
function plot_coup_map(Jmat::Vector{Vector{Union{Float64,Vector{Int64}}}};
    L::Int=1273,
    idx_max::Int=400,
    doms_label::Vector{S}=["NTD", "RBD", "CTD1", "CTD2", "FP", "HR1", "CH", "CD", "HR2"],
    doms_edges::Vector{UnitRange{Int}}=[14:305, 319:541, 542:590, 591:690, 788:806, 910:984, 985:1034, 1035:1067, 1163:1212],
    doms_col::Vector{S}=["cyan", "blue", "orange", "yellow", "red", "green", "lightgreen", "hotpink", "violet"],
    box_height::Float64=0.06,
    box_width::Float64=0.06) where {S<:AbstractString}

    idx_i = [Jmat[k][1][1] for k in 1:idx_max]
    idx_j = [Jmat[k][1][2] for k in 1:idx_max]
    Jval = [Jmat[k][2] for k in 1:idx_max]

    domains = [(doms_edges[n][1], doms_edges[n][end], doms_label[n], doms_col[n]) for n in eachindex(doms_label)]

    fig, ax = subplots(figsize=(10, 8))

    # --- scatter plot ---
    pc = ax.scatter(idx_i, idx_j, c=log.(Jval), s=5, alpha=0.8)
    ax.set_xlim([1, L + 12])
    ax.set_ylim([1, L + 12])
    ax.set_xlabel("Mutated residue i", fontsize=16, labelpad=60)
    ax.set_ylabel("Background residue j", fontsize=16, labelpad=60)
    ax.tick_params(axis="both", labelsize=12)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label(L"\log J_{ij}", fontsize=16)
    cbar.ax.tick_params(labelsize=12)

    # --- X-axis domain boxes positioning ---
    x_bottom = -(box_height + 0.002)  # position so top edge touches the axis

    # --- grey background box for x-axis ---
    bgx = mpatches.Rectangle(
        (1, x_bottom),
        L - 1,
        box_height,
        facecolor="lightgrey",
        edgecolor="none",
        transform=ax.get_xaxis_transform(),
        clip_on=false,
        zorder=0
    )
    ax.add_patch(bgx)

    # --- colored domain boxes for x-axis ---
    for (start, stop, label, color) in domains
        rect = mpatches.FancyBboxPatch(
            (start, x_bottom),
            stop - start,
            box_height,
            boxstyle="round4,pad=0.0",
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            transform=ax.get_xaxis_transform(),
            clip_on=false,
            zorder=1
        )
        ax.add_patch(rect)

        # domain label INSIDE the x-axis box
        ax.text(
            (start + stop) / 2,
            x_bottom + box_height / 2,
            label,
            ha="center", va="center",
            transform=ax.get_xaxis_transform(),
            fontsize=10, fontweight="bold",
            color="white",
            path_effects=[pyimport("matplotlib.patheffects").withStroke(linewidth=2, foreground="black")]
        )
    end

    # --- Y-axis domain boxes positioning ---
    y_left = -(box_width + 0.002)  # position so right edge touches the axis

    # --- grey background box for y-axis ---
    bgy = mpatches.Rectangle(
        (y_left, 1),
        box_width,
        L - 1,
        facecolor="lightgrey",
        edgecolor="none",
        transform=ax.get_yaxis_transform(),
        clip_on=false,
        zorder=0
    )
    ax.add_patch(bgy)

    # --- colored domain boxes for y-axis ---
    for (start, stop, label, color) in domains
        rect = mpatches.FancyBboxPatch(
            (y_left, start),
            box_width,
            stop - start,
            boxstyle="round4,pad=0.0",
            facecolor=color,
            edgecolor="black",
            linewidth=0.5,
            transform=ax.get_yaxis_transform(),
            clip_on=false,
            zorder=1
        )
        ax.add_patch(rect)

        # domain label INSIDE the y-axis box
        ax.text(
            y_left + box_width / 2,
            (start + stop) / 2,
            label,
            ha="center", va="center", rotation=90,
            transform=ax.get_yaxis_transform(),
            fontsize=10, fontweight="bold",
            color="white",
            path_effects=[pyimport("matplotlib.patheffects").withStroke(linewidth=2, foreground="black")]
        )
    end

    # --- create secondary x-axis for ticks ---
    ax_x_ticks = ax.twiny()
    ax_x_ticks.set_xlim(ax.get_xlim())
    ax_x_ticks.set_xticks(collect(0:100:L))
    ax_x_ticks.set_xticklabels(collect(0:100:L), fontsize=12)
    ax_x_ticks.spines["top"].set_position(("axes", x_bottom - 0.008))

    for spine in ["top", "bottom", "left", "right"]
        ax_x_ticks.spines[spine].set_visible(false)
    end
    ax_x_ticks.tick_params(axis="x", which="both", top=true, bottom=false, pad=-15)

    # --- create secondary y-axis for ticks ---
    ax_y_ticks = ax.twinx()
    ax_y_ticks.set_ylim(ax.get_ylim())
    ax_y_ticks.set_yticks(collect(0:100:L))
    ax_y_ticks.set_yticklabels(collect(0:100:L), fontsize=12)
    ax_y_ticks.spines["right"].set_position(("axes", y_left - 0.008))

    for spine in ["top", "bottom", "left", "right"]
        ax_y_ticks.spines[spine].set_visible(false)
    end
    ax_y_ticks.tick_params(axis="y", which="both", right=true, left=false, pad=-15, labelrotation=90)

    # --- hide original ticks and labels ---
    ax.tick_params(axis="x", which="both", bottom=false, labelbottom=false)
    ax.tick_params(axis="y", which="both", left=false, labelleft=false)

    fig.subplots_adjust(bottom=0.12, left=0.12)
    fig.tight_layout()

    return fig, ax

end

# Load inferred couplings
Jtab = CSV.read("results/jcoup_l1_1dot5em4_unif.csv", DataFrame)

# Compute Frobenius norm of couplings
Jfrob = SC2Epistasis.frob_norm(Jtab)
# Convert the Frobenius norm dictionary into an array
Jmat = SC2Epistasis.Jmat(Jfrob)

# Plot the coupling map
fig, ax = plot_coup_map(Jmat; idx_max=400) # only plot top 400 couplings
savefig("results/figures/si/fig_s10.pdf")
close(fig)