""" Script to generate figure 6 of the SI """

# Import packages
using DataFrames, CSV, PyPlot, PyCall, Statistics
@pyimport adjustText

function plot_jdfit(jdfit_muts::Vector{Any};
    nrows::Int=2,
    ncols::Int=2,
    xsize::Float64=6.0,
    ysize::Float64=8.0,
    wspace::Float64=-0.3,
    hspace::Float64=0.15,
    left::Float64=0.02,
    right::Float64=0.98,
    top::Float64=0.97,
    bottom::Float64=0.08,
    x_supy::Float64=0.01,
    y_supx::Float64=0.07)

    @assert nrows * ncols == length(jdfit_muts)
    fig, ax = subplots(nrows, ncols, figsize=(xsize * ncols, ysize * nrows))

    for n in eachindex(jdfit_muts)

        jdfit = jdfit_muts[n]
        n_c = length(jdfit.clades)
        aa = unique(jdfit.sj)
        n_aa = length(aa)
        idx = [zeros(Bool, n_c) for a in 1:n_aa]
        texts = PyObject[]
        for k in 1:n_aa
            idx[k] .= (jdfit.sj .== aa[k])
            ax[n].scatter(jdfit.J[idx[k]], jdfit.fit[idx[k]])
            av = mean(jdfit.fit[idx[k]])
            unc = sqrt(sum(jdfit.s_fit[idx[k]] .^ 2) / sum(idx[k]))
            ax[n].errorbar(jdfit.J[idx[k]][1], av, yerr=unc, fmt="*", markersize=11, elinewidth=2.5, capsize=5, alpha=0.6)

            # Add initial horizontal offset based on position
            x_range = maximum(jdfit.J) - minimum(jdfit.J)
            for c in eachindex(jdfit.clades[idx[k]])
                # Determine offset direction: right if on left half, left if on right half
                x_offset = (jdfit.J[idx[k]][c] < mean(jdfit.J)) ? 0.08 * x_range : -0.08 * x_range
                push!(texts, ax[n].text(jdfit.J[idx[k]][c] + x_offset, jdfit.fit[idx[k]][c],
                    jdfit.clades[idx[k]][c], fontsize=10))
            end
        end
        # Call adjust_text separately for each subplot with stronger repulsion
        adjustText.adjust_text(texts, ax=ax[n], only_move=Dict("text" => "xy", "points" => ""),
            force_text=(0.8, 0.8), force_points=(1.5, 1.5),
            expand_text=(1.5, 1.5), expand_points=(2.0, 2.0))
        ax[n].set_box_aspect(4 / 3)
        ax[n].set_xlabel("i=$(jdfit.i) " * ", j=$(jdfit.j)", fontsize=14)
        ax[n].set_ylabel("Δf($(jdfit.si_wt) → $(jdfit.si))", fontsize=14)
        ax[n].tick_params(labelsize=12)

    end

    # Minimize white space between subplots
    fig.subplots_adjust(wspace=wspace, hspace=hspace, left=left, right=right, top=top, bottom=bottom)

    fig.supxlabel("Coupling J", fontsize=16, ha="center", y=x_supy)
    fig.supylabel("Mutation Δf", fontsize=16, va="center", x=y_supx)

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
muts_list = ["I68V", "P139S", "A419S", "S1003I"]
# List of interacting background mismatching residues
int_res = [69, 83, 417, 764]

# Compute structs for plotting
jdfit_muts = SC2Epistasis.coup_dfit(Jtab, dfit_prot, cdiff_prot, muts_list, int_res)

# Make the plots
fig, ax = plot_jdfit(jdfit_muts; nrows=2, ncols=2, wspace=-0.2, hspace=0.15, left=0.02, right=0.98, y_supx=0.02)
fig.savefig("results/figures/si/fig_s6.pdf")
close(fig)