""" Script for generating panel A of figure 2 """

# Import packages
using DataFrames, CSV, PyPlot, Printf, Statistics, PyCall
@pyimport adjustText

# Plotting functions
function format_legend(label_vec::Vector{S}, group_size::Int=10) where {S<:AbstractString}

    groups = [join(label_vec[i:min(i + group_size - 1, end)], ", ")
              for i in 1:group_size:length(label_vec)]
    return join(groups, "\n")

end

function plot_jdfit(jdfit)

    fig, ax = subplots(figsize=(6, 8))
    n_c = length(jdfit.clades)
    aa = unique(jdfit.sj)
    n_aa = length(aa)
    idx = [zeros(Bool, n_c) for a in 1:n_aa]
    x = fill(0.0, n_c)
    x[jdfit.sj.==aa[1]] .= 1.0
    leg = String[]
    for k in 1:n_aa
        idx[k] .= (jdfit.sj .== aa[k])
        ax.scatter(x[idx[k]], jdfit.fit[idx[k]])
        av = mean(jdfit.fit[idx[k]])
        unc = sqrt(sum(jdfit.s_fit[idx[k]] .^ 2) / sum(idx[k]))
        ax.errorbar(x[idx[k]][1], av, yerr=unc, fmt="*", markersize=11, elinewidth=2.5, capsize=5, alpha=0.6)
        label_vec = jdfit[1][idx[k]]
        label = format_legend(label_vec, 10) * ": " * aa[k]
        push!(leg, label)
    end
    ax.set_xticks([x[1], x[end]])
    ax.set_xticklabels(["σⱼ=S", "σⱼ=R"], fontsize=16)
    ax.set_xlim(-0.25, 1.25)
    ax.set_ylabel("Δf($(jdfit.si_wt) → $(jdfit.si))", fontsize=16)
    ax.legend(leg, fontsize=13, loc="upper left")

    fig.tight_layout()

    return fig, ax

end

# Script body
################################---------------------------------

# Load inferred couplings
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)
# Load dataframe of fitness discrepancies
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
# Load dataframe of clade differences
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Select protein
prot = "S"
dfit_prot = delta_fit[delta_fit.prot.==prot, :]
cdiff_prot = clade_diff[clade_diff.prot.==prot, :]

# Mutation
mut = "P384L"
# Interacting site
res = 408

# Compute struct for plotting
jdfit = SC2Epistasis.coup_dfit(Jtab, dfit_prot, cdiff_prot, [mut], [res])[1]

# Make the plot
fig, ax = plot_jdfit(jdfit)
fig.savefig("results/figures/fig_epi_pic_B.pdf")
close("all")