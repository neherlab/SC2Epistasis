""" Script to generate composite figure S9 (panels A–D) """

# Import packages
using DataFrames, CSV, PdbTool, PyPlot, PyCall, LinearAlgebra

const mpatches = pyimport("matplotlib.patches")
const mpatheffects = pyimport("matplotlib.patheffects")

# ---------------------------------------------------------------------------
# Shared data dependencies
# ---------------------------------------------------------------------------
include(joinpath(@__DIR__, "usa_uk_energy.jl"))  # defines fit_merge, e0
include(joinpath(@__DIR__, "init_inf.jl"))        # defines delta_fit, data, optx, qform

# ---------------------------------------------------------------------------
# Optimization functions (panel B)
# ---------------------------------------------------------------------------

function optimize_unif(λ1::Vector{Float64}, optx::Vector{OptX}, qform::Vector{QForm};
    λ2::Float64=1.0e-6,
    epsconv::Float64=1.0e-10,
    verbose::Bool=true,
    maxiter::Int64=5000)

    tot_num_j = sum([optx[m].num_j for m in eachindex(optx)])
    e_vec = zeros(length(λ1))
    loss = zeros(length(λ1))
    zf = zeros(length(λ1))

    for l in eachindex(λ1)
        l1_reg = fill(λ1[l], tot_num_j)
        loss[l] = SC2Epistasis.learn(optx, qform, l1_reg; λ2=λ2, epsconv=epsconv, verbose=verbose, maxiter=maxiter)
        e_vec[l] = SC2Epistasis.energy(optx, qform)
        J = vcat([optx[m].J for m in 1:length(optx)]...)
        zf[l] = sum(abs.(J) .== 0) / length(J)
    end

    return DataFrame(l1=λ1, energy=e_vec, loss=loss, tot_zf=zf)
end

function optimize_3d_lin(λ1::Vector{Float64}, data::Vector{Data}, optx::Vector{OptX}, qform::Vector{QForm}, dist::Vector{Float64};
    λ2::Float64=1.0e-6,
    epsconv::Float64=1.0e-10,
    verbose::Bool=true,
    maxiter::Int64=5000)

    e_vec = zeros(length(λ1))
    loss = zeros(length(λ1))
    zf = zeros(length(λ1))

    for l in eachindex(λ1)
        l1_reg = SC2Epistasis.threeD_l1(λ1[l], optx, data, dist)
        loss[l] = SC2Epistasis.learn(optx, qform, l1_reg; λ2=λ2, epsconv=epsconv, verbose=verbose, maxiter=maxiter)
        e_vec[l] = SC2Epistasis.energy(optx, qform)
        J = vcat([optx[m].J for m in 1:length(optx)]...)
        zf[l] = sum(abs.(J) .== 0) / length(J)
    end

    return DataFrame(l1=λ1, energy=e_vec, loss=loss, tot_zf=zf)
end

# ---------------------------------------------------------------------------
# Coupling-map helper: draws onto an existing axis (adapted from plot_coup_map.jl)
# ---------------------------------------------------------------------------

function plot_coup_map_on_ax(fig, ax, Jmat::Vector{Vector{Union{Float64,Vector{Int64}}}};
    L::Int=1273,
    idx_max::Int=400,
    doms_label::Vector{S}=["NTD", "RBD", "CTD1", "CTD2", "FP", "HR1", "CH", "CD", "HR2"],
    doms_edges::Vector{UnitRange{Int}}=[14:305, 319:541, 542:590, 591:690, 788:806, 910:984, 985:1034, 1035:1067, 1163:1212],
    doms_col::Vector{S}=["cyan", "blue", "orange", "yellow", "red", "green", "lightgreen", "hotpink", "violet"],
    box_height::Float64=0.06,
    box_width::Float64=0.06,
    tick_step::Int=100) where {S<:AbstractString}

    idx_i = [Jmat[k][1][1] for k in idx_max:-1:1]
    idx_j = [Jmat[k][1][2] for k in idx_max:-1:1]
    Jval  = [Jmat[k][2]    for k in idx_max:-1:1]

    domains = [(doms_edges[n][1], doms_edges[n][end], doms_label[n], doms_col[n]) for n in eachindex(doms_label)]

    mcolors = pyimport("matplotlib.colors")
    mcm = pyimport("matplotlib.cm")
    np = pyimport("numpy")
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "trunc_YlOrRd", mcm.YlOrRd(np.linspace(0.3, 1.0, 256)))
    pc = ax.scatter(idx_i, idx_j, c=log.(Jval), s=70*max.(0.1, log.(Jval)) + 5.0*ones(length(Jval)), alpha=0.6, cmap=cmap)
    ax.set_xlim([1, L + 12])
    ax.set_ylim([1, L + 12])
    ax.set_xlabel("Mutated residue i", fontsize=12, labelpad=50)
    ax.set_ylabel("Background residue j", fontsize=12, labelpad=50)
    ax.tick_params(axis="both", labelsize=10)
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label(L"\log J_{ij}", fontsize=12)
    cbar.ax.tick_params(labelsize=10)

    x_bottom = -(box_height + 0.002)

    bgx = mpatches.Rectangle(
        (1, x_bottom), L - 1, box_height,
        facecolor="lightgrey", edgecolor="none",
        transform=ax.get_xaxis_transform(), clip_on=false, zorder=0)
    ax.add_patch(bgx)

    for (start, stop, label, color) in domains
        rect = mpatches.FancyBboxPatch(
            (start, x_bottom), stop - start, box_height,
            boxstyle="round4,pad=0.0", facecolor=color, edgecolor="black",
            linewidth=0.5, transform=ax.get_xaxis_transform(), clip_on=false, zorder=1)
        ax.add_patch(rect)
        ax.text((start + stop) / 2, x_bottom + box_height / 2, label,
            ha="center", va="center", transform=ax.get_xaxis_transform(),
            fontsize=10, fontweight="bold", color="white",
            path_effects=[mpatheffects.withStroke(linewidth=2, foreground="black")])
    end

    y_left = -(box_width + 0.002)

    bgy = mpatches.Rectangle(
        (y_left, 1), box_width, L - 1,
        facecolor="lightgrey", edgecolor="none",
        transform=ax.get_yaxis_transform(), clip_on=false, zorder=0)
    ax.add_patch(bgy)

    for (start, stop, label, color) in domains
        rect = mpatches.FancyBboxPatch(
            (y_left, start), box_width, stop - start,
            boxstyle="round4,pad=0.0", facecolor=color, edgecolor="black",
            linewidth=0.5, transform=ax.get_yaxis_transform(), clip_on=false, zorder=1)
        ax.add_patch(rect)
        ax.text(y_left + box_width / 2, (start + stop) / 2, label,
            ha="center", va="center", rotation=90, transform=ax.get_yaxis_transform(),
            fontsize=8, fontweight="bold", color="white",
            path_effects=[mpatheffects.withStroke(linewidth=2, foreground="black")])
    end

    ax_x_ticks = ax.twiny()
    ax_x_ticks.set_xlim(ax.get_xlim())
    ax_x_ticks.set_xticks(collect(0:tick_step:L))
    ax_x_ticks.set_xticklabels(collect(0:tick_step:L), fontsize=10)
    ax_x_ticks.spines["top"].set_position(("axes", x_bottom - 0.008))
    for spine in ["top", "bottom", "left", "right"]
        ax_x_ticks.spines[spine].set_visible(false)
    end
    ax_x_ticks.tick_params(axis="x", which="both", top=true, bottom=false, pad=-15)

    ax_y_ticks = ax.twinx()
    ax_y_ticks.set_ylim(ax.get_ylim())
    ax_y_ticks.set_yticks(collect(0:tick_step:L))
    ax_y_ticks.set_yticklabels(collect(0:tick_step:L), fontsize=10)
    ax_y_ticks.spines["right"].set_position(("axes", y_left - 0.008))
    for spine in ["top", "bottom", "left", "right"]
        ax_y_ticks.spines[spine].set_visible(false)
    end
    ax_y_ticks.tick_params(axis="y", which="both", right=true, left=false, pad=-15, labelrotation=90)

    ax.tick_params(axis="x", which="both", bottom=false, labelbottom=false)
    ax.tick_params(axis="y", which="both", left=false, labelleft=false)
end

# ---------------------------------------------------------------------------
# Compute panel A data (z-score histograms)
# ---------------------------------------------------------------------------

rename!(fit_merge, [:delta_fitness => :fit1, :delta_fitness_1 => :fit2,
                    :uncertainty => :std_fit1, :uncertainty_1 => :std_fit2])
z_country = SC2Epistasis.z_dfit(fit_merge)
z_clades  = SC2Epistasis.z_dfit(delta_fit_prot)

# ---------------------------------------------------------------------------
# Compute panel B data (model inference energy vs λ)
# ---------------------------------------------------------------------------

af_pdb = PdbTool.parsePdb("data/PDB/Nprot/af_nprot.pdb")
Nprot_seq = read("data/ref_seq/Nprot.txt", String)
_ = SC2Epistasis.map_pdb!(af_pdb, [Nprot_seq]; mappedTo="data/ref_seq/Nprot.txt")
dist = SC2Epistasis.threedist(optx, data, af_pdb)

λ1 = [1.0e-6, 1.0e-5, 5.0e-5, 1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2]
fit_df_unif = optimize_unif(λ1, optx, qform; epsconv=1.0e-10)
fit_df_3d   = optimize_3d_lin(λ1, data, optx, qform, dist; epsconv=1.0e-10)

# ---------------------------------------------------------------------------
# Load coupling parameters for panels C and D
# ---------------------------------------------------------------------------

Jtab = CSV.read("results/jcoup_nprot_l1_1em3_3d.csv", DataFrame)
Jfrob = SC2Epistasis.frob_norm(Jtab)
Jmat_data = SC2Epistasis.Jmat(Jfrob)

# ---------------------------------------------------------------------------
# Build composite figure
# ---------------------------------------------------------------------------

fig, axes = subplots(2, 2, figsize=(14, 12))
ax_A, ax_B = axes[1, 1], axes[1, 2]
ax_C, ax_D = axes[2, 1], axes[2, 2]

# Panel A — z-score histogram
ax_A.hist(z_clades, bins=40, density=true, alpha=0.8, label="Clades", edgecolor="black")
ax_A.hist(z_country, bins=20, density=true, alpha=0.6, label="USA-UK", edgecolor="black")
ax_A.legend(fontsize=11)
ax_A.set_yscale("log")
ax_A.set_ylabel("Density", fontsize=12)
ax_A.set_xlabel("Fitness z-score", fontsize=12)
ax_A.tick_params(axis="both", labelsize=10)

# Panel B — energy vs regularization
labels_B = ["Uniform", "3d"]
linestyles_B = ["-", "--"]
shades_B = [get_cmap("Blues")(0.4 + 0.2 * i) for i in 0:1]
for (n, fit_df) in enumerate([fit_df_unif, fit_df_3d])
    ax_B.plot(fit_df.l1, fit_df.energy, marker=".", markersize=5.0,
        linestyle=linestyles_B[n], color=shades_B[n], label=labels_B[n])
end
ax_B.hlines(e0, fit_df_3d.l1[1], fit_df_3d.l1[end], color="black", linestyle="-", label="Estimated noise")
ax_B.set_xscale("log")
ax_B.set_xlabel("λ₁ (Sparsity Regularization)", fontsize=12)
ax_B.set_ylabel("Energy (Mean Squared Error)", fontsize=12)
ax_B.legend(loc="best", fontsize=10)
ax_B.tick_params(axis="both", labelsize=10)

# Panel C — histogram of coupling parameters
ax_C.hist(Jtab.J, bins=21, density=false, alpha=0.8, edgecolor="black")
ax_C.set_yscale("log")
ax_C.set_xlabel("Coupling parameters " * L"J_{ij}(σᵢ,σⱼ)", fontsize=12)
ax_C.set_ylabel("Frequency", fontsize=12)
ax_C.tick_params(axis="both", labelsize=10)

# Panel D — coupling map
plot_coup_map_on_ax(fig, ax_D, Jmat_data;
    idx_max=75,
    L=421,
    doms_label=["NTD", "SR", "CTD", "Tetramerization"],
    doms_edges=[48:174, 176:206, 247:364, 400:419],
    doms_col=["cyan", "yellow", "red", "green"],
    tick_step=50
)

# Panel labels
for (ax, label) in zip([ax_A, ax_B, ax_C, ax_D], ["A", "B", "C", "D"])
    ax.text(-0.12, 1.05, label, transform=ax.transAxes,
        fontsize=18, fontweight="bold", va="top", ha="right")
end

fig.tight_layout()
fig.savefig("results/figures/si/fig_s9.pdf")
close(fig)
