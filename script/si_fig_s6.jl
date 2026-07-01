""" Script to generate figure 6 of the SI """

# Import packages
using DataFrames, CSV, PyPlot, PyCall, Statistics
@pyimport adjustText
@isdefined(mpatches) || (const mpatches = PyPlot.PyCall.pyimport("matplotlib.patches"))

function plot_jdfit!(axs_flat, jdfit_muts, background_states)

    for n in eachindex(jdfit_muts)
        ax = axs_flat[n]
        jdfit = jdfit_muts[n]
        n_c = length(jdfit.clades)
        aa = unique(jdfit.sj)
        n_aa = length(aa)
        idx = [zeros(Bool, n_c) for a in 1:n_aa]
        texts = PyObject[]
        for k in 1:n_aa
            idx[k] .= (jdfit.sj .== aa[k])
            ax.scatter(jdfit.J[idx[k]], jdfit.fit[idx[k]])
            av = mean(jdfit.fit[idx[k]])
            unc = sqrt(sum(jdfit.s_fit[idx[k]] .^ 2) / sum(idx[k]))
            ax.errorbar(jdfit.J[idx[k]][1], av, yerr=unc, fmt="*", markersize=11, elinewidth=2.5, capsize=5, alpha=0.6)
            x_range = maximum(jdfit.J) - minimum(jdfit.J)
            for (ci, c) in enumerate(eachindex(jdfit.clades[idx[k]]))
                if jdfit.J[idx[k]][c] < -0.5
                    x_offset = 0.03 * x_range
                elseif jdfit.J[idx[k]][c] > 0.5
                    x_offset = -0.12 * x_range
                else
                    x_offset = ci%2==0 ? 0.03 * x_range : -0.12 * x_range
                end
                push!(texts, ax.text(jdfit.J[idx[k]][c] + x_offset, jdfit.fit[idx[k]][c],
                    jdfit.clades[idx[k]][c], fontsize=10))
            end
        end
        adjustText.adjust_text(texts, ax=ax, only_move=Dict("text" => "xy", "points" => ""),
            force_text=(0.2, 1.2), force_points=(0.5, 0.8),
            expand_text=(1.2, 1.2), expand_points=(2.0, 2.0))
        ax.set_box_aspect(3 / 3)
        states_str = join(background_states[n], ",")
        ax.set_xlabel("i=$(jdfit.i), j=$(jdfit.j),  \$\\sigma_j = $(states_str)\$", fontsize=14)
        ax.set_ylabel("Δf($(jdfit.si_wt) → $(jdfit.si))", fontsize=14)
        ax.tick_params(labelsize=12)
    end

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
dfit_prot = delta_fit[delta_fit.prot .== prot, :]
cdiff_prot = clade_diff[clade_diff.prot .== prot, :]

# List of mutations to be plotted
muts_list = ["I68V", "P139S", "A419S", "S1003I"]
# List of interacting background mismatching residues
int_res = [69, 83, 417, 764]

# Compute structs for plotting
jdfit_muts = SC2Epistasis.coup_dfit(Jtab, dfit_prot, cdiff_prot, muts_list, int_res)

# List of mismatches between group of clades
background_states = [["-", "H"], ["V", "A"], ["N", "K"], ["N", "K"]]

# Make the plots

# Initialize figure and axes
fig, ax = subplots(2, 2, figsize=(10, 10))
ax_flat = vec(ax)

# Plot the data on the axes
plot_jdfit!(ax_flat, jdfit_muts, background_states)

# Super-labels for figure
fig.text(0.5, 0.01, L"Interaction parameters $J_{ij}$", ha="center", fontsize=16)
fig.text(0.02, 0.5, L"Mutation $\Delta f_i$", va="center", rotation="vertical", fontsize=16)

# Tight layout
fig.tight_layout(rect=[0.04, 0.03, 0.98, 0.98])

# Save and close the figure
fig.savefig("results/figures/si/fig_s6.pdf")
close(fig)