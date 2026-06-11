""" Combined figure: panels A, B, C of figure 2 """

# Import packages
using DataFrames, CSV, PyPlot, Printf, Statistics
using PyCall

# Preliminary data for panel C
include("z_struct_prots.jl")
include("z_clades_regions.jl")

# Helper functions for panel B
function format_legend(label_vec::Vector{S}, group_size::Int=10) where {S<:AbstractString}
    groups = [join(label_vec[i:min(i + group_size - 1, end)], ", ")
              for i in 1:group_size:length(label_vec)]
    return join(groups, "\n")
end

function plot_jdfit!(ax, jdfit; markersize=60, alpha=0.7, label_fs=12)
    n_c = length(jdfit.clades)
    aa = unique(jdfit.sj)
    n_aa = length(aa)
    idx = [zeros(Bool, n_c) for a in 1:n_aa]
    x = fill(0.0, n_c)
    x[jdfit.sj.==aa[1]] .= 1.0
    for k in 1:n_aa
        idx[k] .= (jdfit.sj .== aa[k])
        ax.scatter(x[idx[k]], jdfit.fit[idx[k]], s=markersize, alpha=alpha)
        av = mean(jdfit.fit[idx[k]])
        unc = sqrt(sum(jdfit.s_fit[idx[k]] .^ 2) / sum(idx[k]))
        ax.errorbar(x[idx[k]][1], av, yerr=unc, fmt="*", markersize=11, elinewidth=2.5, capsize=5, alpha=0.6)
    end
    ax.tick_params(axis="y", labelsize=12)
    ax.set_xticks([x[end], x[1]])
    ax.set_xticklabels(["σⱼ=S", "σⱼ=R"], fontsize=label_fs)
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylabel("Mutation effect Δf($(jdfit.si_wt) → $(jdfit.si))", fontsize=label_fs)
    ax.set_xlabel("Background", fontsize=label_fs)
end

###########################################------------------------------------

# Load data for panels A and B
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
prot = "S"
dfit_prot = delta_fit[delta_fit.prot.==prot, :]

# Panel A data
cnt_thr = 10.0
dfit_prot_A = dfit_prot[(dfit_prot.exp_count1.>=cnt_thr).&(dfit_prot.exp_count2.>=cnt_thr), :]
clade_pair = ["21J", "21L"]
dfit_c1c2 = dfit_prot_A[(dfit_prot_A.clade1.==clade_pair[1]).&(dfit_prot_A.clade2.==clade_pair[2]), :]
mut = "P384L"
idx_mut = findall(dfit_c1c2.aa_mut .== mut)
fmax = maximum(vcat(dfit_c1c2.fit1, dfit_c1c2.fit2))
fmin = minimum(vcat(dfit_c1c2.fit1, dfit_c1c2.fit2))
r = cor(dfit_c1c2.fit1, dfit_c1c2.fit2)
n = length(dfit_c1c2.fit1)

# Panel B data
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)
cdiff_prot = clade_diff[clade_diff.prot.==prot, :]
res = 408
jdfit = SC2Epistasis.coup_dfit(Jtab, dfit_prot, cdiff_prot, [mut], [res])[1]

###########################################------------------------------------

label_fs = 14  # axis label fontsize for all panels

# Build figure: 2 rows x 2 cols
# Col 0: A (top), B (bottom)
# Col 1: C_top and C_bot tightly stacked, sharing x-axis
fig = figure(figsize=(12, 10))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.4)
gs_right = py"$gs[:, 1]".subgridspec(2, 1, hspace=0.05)  # tightly stacked C panels, full right column

ax_A  = fig.add_subplot(gs[1, 1])
ax_B  = fig.add_subplot(gs[2, 1])
ax_C1 = fig.add_subplot(gs_right[1, 1])
ax_C2 = fig.add_subplot(gs_right[2, 1], sharex=ax_C1)

# Panel A
shift = 0.2
diag = collect(fmin:0.2:fmax)
ax_A.scatter(dfit_c1c2.fit1, dfit_c1c2.fit2, s=12, alpha=0.4, rasterized=true)
ax_A.scatter(dfit_c1c2.fit1[idx_mut], dfit_c1c2.fit2[idx_mut], s=12, color="red")
ax_A.scatter(dfit_c1c2.fit1[idx_mut], dfit_c1c2.fit2[idx_mut], s=120,
    facecolors="none", edgecolors="red", linewidths=1.5)
ax_A.annotate(mut, (dfit_c1c2.fit1[idx_mut][1], dfit_c1c2.fit2[idx_mut][1]),
    xytext=(8, 8), textcoords="offset points", fontsize=11, color="red")
ax_A.set_xlim(fmin - shift, fmax + shift)
ax_A.set_ylim(fmin - shift, fmax + shift)
ax_A.plot(diag, diag, color="orange")
ax_A.set_xlabel(clade_pair[1] * " (Delta) fitness effect", fontsize=label_fs)
ax_A.set_ylabel(clade_pair[2] * " (BA.2) fitness effect", fontsize=label_fs)
ax_A.text(0.2, 0.85, "r = " * @sprintf("%.2f", r) * "\nn = " * string(n),
    fontsize=12, transform=ax_A.transAxes)

# Panel B
plot_jdfit!(ax_B, jdfit; markersize=60, alpha=0.7, label_fs=label_fs)

# Panel C — top histogram (no x-tick labels, shared with bottom)
ax_C1.hist(z_clades, bins=40, density=true, alpha=0.8, label="Clades", edgecolor="black")
ax_C1.hist(z_country, bins=20, density=true, alpha=0.6, label="USA-UK", edgecolor="black")
ax_C1.legend(fontsize=13)
ax_C1.tick_params(axis="y", labelsize=14)
ax_C1.tick_params(axis="x", labelbottom=false)
ax_C1.set_yscale("log")

# Panel C — bottom histogram
ax_C2.hist(z_nsyn, bins=50, alpha=0.6, density=true, label="Non-synonymous")
ax_C2.hist(z_syn, bins=50, alpha=0.5, density=true, label="Synonymous")
ax_C2.hist(z_stop, bins=50, alpha=0.4, density=true, label="Stop codon")
ax_C2.set_xlabel("Fitness z-score", fontsize=label_fs)
ax_C2.legend(fontsize=11)
ax_C2.tick_params(axis="both", labelsize=14)
ax_C2.set_yscale("log")

# Common y-label for C panels (positioned at mid-height of the right column)
fig.text(0.5, 0.5, "Density", va="center", rotation="vertical", fontsize=label_fs)

# Panel labels
for (ax, lbl) in zip([ax_A, ax_B, ax_C1], ["A", "B", "C"])
    ax.text(-0.12, 1.02, lbl, transform=ax.transAxes, fontsize=16, fontweight="bold", va="bottom")
end

fig.savefig("results/figures/fig_epi_pic.pdf")
close(fig)
