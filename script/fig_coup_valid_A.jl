""" Script to generate figure panel A of figure 6 """

# Import packages
using DataFrames, CSV, PyPlot, PyCall, Statistics
@pyimport adjustText

# Plotting functions
function format_legend(label_vec::Vector{S}, group_size::Int=10) where {S<:AbstractString}

    groups = [join(label_vec[i:min(i + group_size - 1, end)], ", ")
              for i in 1:group_size:length(label_vec)]
    return join(groups, "\n")

end

function plot_jdfit(jdfit_muts::Vector{Any};
    nrows::Int=2,
    ncols::Int=2,
    xsize::Float64=6.0,
    ysize::Float64=8.0,
    wspace::Float64=-0.1,
    hspace::Float64=0.1)

    @assert nrows * ncols == length(jdfit_muts)
    fig, ax = subplots(nrows, ncols, figsize=(xsize * ncols, ysize * nrows))

    for n in eachindex(jdfit_muts)

        jdfit = jdfit_muts[n]
        n_c = length(jdfit.clades)
        aa = unique(jdfit.sj)
        n_aa = length(aa)
        idx = [zeros(Bool, n_c) for a in 1:n_aa]
        leg = String[]
        texts = PyObject[]
        for k in 1:n_aa
            idx[k] .= (jdfit.sj .== aa[k])
            ax[n].scatter(jdfit.J[idx[k]], jdfit.fit[idx[k]])
            av = mean(jdfit.fit[idx[k]])
            unc = sqrt(sum(jdfit.s_fit[idx[k]] .^ 2) / sum(idx[k]))
            ax[n].errorbar(jdfit.J[idx[k]][1], av, yerr=unc, fmt="*", markersize=11, elinewidth=2.5, capsize=5, alpha=0.6)
            for c in eachindex(jdfit.clades[idx[k]])
                push!(texts, ax[n].text(jdfit.J[idx[k]][c] + 0.1, jdfit.fit[idx[k]][c], jdfit.clades[idx[k]][c], fontsize=10))
            end
            label_vec = jdfit[1][idx[k]]
            label = format_legend(label_vec, 10) * ": " * aa[k]
            push!(leg, label)
        end
        # Call adjust_text separately for each subplot
        adjustText.adjust_text(texts, ax=ax[n], only_move="xy")
        ax[n].set_box_aspect(4 / 3)
        ax[n].set_xlabel("i=$(jdfit.i) " * ", j=$(jdfit.j)", fontsize=14)
        ax[n].set_ylabel("Δf($(jdfit.si_wt) → $(jdfit.si))", fontsize=14)
        ax[n].tick_params(labelsize=12)
        ax[n].legend(leg, fontsize=10)

    end

    fig.supxlabel("Coupling J", fontsize=16, ha="center")
    fig.supylabel("Mutation Δf", fontsize=16, va="center")

    # Minimize white space between subplots
    fig.subplots_adjust(wspace=wspace, hspace=hspace)

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

# List of mutations to be plotted
muts_list = ["R21G", "P384L", "P499L", "R683Q"]
# List of interacting background mismatching residues
int_res = [19, 408, 445, 681]

# Compute structs for plotting
jdfit_muts = SC2Epistasis.coup_dfit(Jtab, dfit_prot, cdiff_prot, muts_list, int_res)

# Make the plots
fig, ax = plot_jdfit(jdfit_muts; nrows=2, ncols=2, wspace=-0.1, hspace=0.1)
fig.tight_layout()
fig.savefig("results/figures/fig_coup_valid_A.pdf")
close(fig)