""" Script to generate figure 6 """

# Import packages
using DataFrames, CSV, PyPlot, PyCall, Statistics
@pyimport adjustText

# Functions
function coup_dfit(Jtab::DataFrame, dfit::DataFrame, cdiff::DataFrame, muts::Vector{S}, res::Vector{Int}) where {S<:AbstractString}

    j_dfit = []
    for m in eachindex(muts)
        site = parse(Int, muts[m][2:end-1])
        # Extract couplings involving the mutation site
        Jmut = Jtab[(Jtab.σᵢ_wt.==string(muts[m][1])).&(Jtab.i.==site).&(Jtab.σᵢ.==string(muts[m][end])).&(Jtab.j.==res[m]), :]
        # Extract fitness discrepancies for the mutation
        dfit_mut = dfit[dfit.aa_mut.==muts[m], :]
        push!(j_dfit, coup_dfit(Jmut, dfit_mut, cdiff, res[m]))
    end

    return j_dfit

end

function coup_dfit(Jmut::DataFrame, dfit_mut::DataFrame, cdiff::DataFrame, res::Int)

    cp_mut = SC2Epistasis.clade_pair_mut(dfit_mut, cdiff)
    cp_mut = cp_mut[cp_mut.site.==res, :]

    clades = unique(vcat(dfit_mut.clade1, dfit_mut.clade2))
    fit = fill(0.0, length(clades))
    std_fit = fill(0.0, length(clades))
    sj = fill("", length(clades))
    J = fill(0.0, length(clades))

    for (i, c) in enumerate(clades)
        if c in dfit_mut.clade1
            idx_c = findall(dfit_mut.clade1 .== c)
            fit[i] = dfit_mut.fit1[idx_c][1]
            std_fit[i] = dfit_mut.std_fit1[idx_c][1]
        else
            idx_c = findall(dfit_mut.clade2 .== c)
            fit[i] = dfit_mut.fit2[idx_c][1]
            std_fit[i] = dfit_mut.std_fit2[idx_c][1]
        end
        if sum(cp_mut.clade1 .== c) != 0
            aa = cp_mut[cp_mut.clade1.==c, :aa_c1][1]
        else
            aa = cp_mut[cp_mut.clade2.==c, :aa_c2][1]
        end
        sj[i] = aa
        J[i] = Jmut[string.(Jmut.σⱼ).==aa, :J][1]
    end

    return (clades=clades, fit=fit, s_fit=std_fit, sj=sj, J=J, si_wt=Jmut.σᵢ_wt[1], si=Jmut.σᵢ[1], i=Jmut.i[1], j=res)

end

# Plotting functions
function format_legend(label_vec::Vector{S}, group_size::Int=10) where {S<:AbstractString}

    groups = [join(label_vec[i:min(i + group_size - 1, end)], ", ")
              for i in 1:group_size:length(label_vec)]
    return join(groups, "\n")

end

function plot_jdfit(jdfit_muts::Vector{Any};
    xsize::Float64=6.0,
    ysize::Float64=8.0,
    wspace::Float64=-0.4,
    hspace::Float64=0.1)

    #@assert nrows * ncol == length(jdfit_muts)
    fig, ax = subplots(2, 4, figsize=(xsize * 4, ysize * 2))

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
        ax[n].legend(leg, fontsize=10)

    end

    fig.supxlabel("Coupling J", fontsize=16, ha="center")
    fig.supylabel("Mutation Δf", fontsize=16, va="center")

    fig.tight_layout()

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
muts_list = ["R21G", "I68V", "P139S", "P384L", "A419S", "P499L", "R683Q", "S1003I"]
# List of interacting background mismatching residues
int_res = [19, 69, 83, 408, 417, 445, 681, 764]

# Compute structs for Plotting
jdfit_muts = coup_dfit(Jtab, dfit_prot, cdiff_prot, muts_list, int_res)

# Make the plots
fig, ax = plot_jdfit(jdfit_muts; xsize=6.0, ysize=8.0, wspace=-0.4, hspace=0.1)
fig.tight_layout()
fig.supylabel("Mutation Δf", fontsize=16, va="center", x=0.065)
fig.savefig("figures/fig_coup_fit.pdf")
