""" Composite figure: panels A, B, C of figure 6 """

# Import packages
using DataFrames, CSV, PyPlot, PyCall, Statistics
@pyimport adjustText
@isdefined(mpatches) || (const mpatches = PyPlot.PyCall.pyimport("matplotlib.patches"))

###############################---------------------------------

# Panel A plotting function
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

# Panel B plotting function
function plot_hist!(ax, av_z, av_dz)
    ax.hist(av_z, bins=30, density=true, alpha=0.7, edgecolor="black")
    ax.hist(av_dz, bins=30, density=true, alpha=0.5, edgecolor="black")
    ax.set_xlabel("Fitness discrepancy z-score", fontsize=14)
    ax.set_ylabel("Density", fontsize=14)
    ax.legend(["Raw discrepancies z", "Model predicted ϵ"], fontsize=13)
    ax.tick_params(labelsize=12)
end

# Panel C helper: merge DMS and fitness discrepancy dataframes
function merge_dms_dfit(dms_shift::DataFrame, dfit::DataFrame, clade_pair::Tuple{S,S};
    cnt_thr1::Float64=10.0, cnt_thr2::Float64=10.0) where {S<:AbstractString}

    dfit_c1c2 = filter(row -> row.clade1 == clade_pair[1] && row.clade2 == clade_pair[2], dfit)
    if isempty(dfit_c1c2)
        dfit_c1c2 = filter(row -> row.clade1 == clade_pair[2] && row.clade2 == clade_pair[1], dfit)
    end
    rename!(dfit_c1c2, :aa_mut => :mutation)
    df_merge = innerjoin(dms_shift, dfit_c1c2, on=:mutation)
    df_merge = df_merge[(df_merge.exp_count1 .>= cnt_thr1) .& (df_merge.exp_count2 .>= cnt_thr2), :]
    return df_merge

end

# Panel C plotting function
function plot_shift!(ax, dms_shift::DataFrame, Jtab::DataFrame,
    dfit_prot::DataFrame, cdiff_prot::DataFrame, clade_pair::Tuple{S,S};
    cnt_thr1::Float64=40.0, cnt_thr2::Float64=20.0,
    doms_label::Vector{S1}=["NTD", "RBD"],
    doms_edges::Vector{UnitRange{Int}}=[14:305, 319:541],
    doms_col::Vector{S1}=["orange", "blue"],
    doms_marker::Vector{S1}=["s", "^"]) where {S<:AbstractString,S1<:AbstractString}

    shift_dfit = merge_dms_dfit(dms_shift, dfit_prot, clade_pair; cnt_thr1=cnt_thr1, cnt_thr2=cnt_thr2)

    muts = shift_dfit.mutation
    res_index_vec = parse.(Int, map(x -> x[2:(end-1)], muts))
    j_ddf = SC2Epistasis.mutcp_ddf(Jtab, muts, cdiff_prot, clade_pair[1], clade_pair[2])
    shift = shift_dfit[!, "score_"*clade_pair[1]] .- shift_dfit[!, "score_"*clade_pair[2]]

    domains = [(doms_edges[n][1], doms_edges[n][end], doms_label[n], doms_col[n]) for n in eachindex(doms_label)]
    domain_masks = [(res_index_vec .>= start) .& (res_index_vec .<= stop) for (start, stop, _, _) in domains]
    other_domain = .!foldl((acc, mask) -> acc .| mask, domain_masks)

    if any(other_domain)
        ρ = cor(shift[other_domain], j_ddf[other_domain])
        ax.scatter(j_ddf[other_domain], shift[other_domain], alpha=0.7, s=15, color="lightgreen", marker="o", label="Other" * " (ρ=" * string(round(ρ, digits=2)) * ")")
    end
    for (n, (start, stop, label, color)) in enumerate(domains)
        in_domain = (res_index_vec .>= start) .& (res_index_vec .<= stop)
        ρ = cor(shift[in_domain], j_ddf[in_domain])
        if any(in_domain)
            ax.scatter(j_ddf[in_domain], shift[in_domain], alpha=0.7, s=15, color=color, marker=doms_marker[n], label=label * " (ρ=" * string(round(ρ, digits=2)) * ")")
        end
    end

    ax.legend(loc="upper left", fontsize=11, frameon=true)

    ax.set_xlabel("Model prediction ΔΔϕ", fontsize=14)
    ax.set_ylabel("Experimental fitness shift", fontsize=14)
    ax.tick_params(labelsize=12)

end

###############################---------------------------------

# Load data

# Shared: couplings, fitness, clade differences
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

prot = "S"
dfit_prot = delta_fit[delta_fit.prot .== prot, :]
cdiff_prot = clade_diff[clade_diff.prot .== prot, :]

# Panel A
muts_list = ["R21G", "P384L", "P499L", "R683Q"]
int_res = [19, 408, 445, 681]
background_states = [["I", "T"], ["P", "R"], ["V", "H", "P"], ["P", "H", "R"]]
jdfit_muts = SC2Epistasis.coup_dfit(Jtab, dfit_prot, cdiff_prot, muts_list, int_res)

# Panel B
z_df = CSV.read("results/model_test_zscore.csv", DataFrame)

# Panel C
dms_shift = CSV.read("data/dms_shift/dms_shift_ba1_ba2_21j.csv", DataFrame)
cpair = ("21J", "21K")

###############################---------------------------------

# Build composite figure: 2 rows × 3 cols
# Cols 0-1: A's 2×2 subplots; col 2: B (top) and C (bottom)
fig = figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, width_ratios=[5, 5, 6], hspace=0.24, wspace=0.22)

axs_A = [fig.add_subplot(gs[i, j]) for i in 1:2, j in 1:2]
axs_A_flat = vec(axs_A)  # row-major: [A[1,1], A[1,2], A[2,1], A[2,2]]
ax_B = fig.add_subplot(gs[1, 3])
ax_C = fig.add_subplot(gs[2, 3])

# Plot each panel
plot_jdfit!(axs_A_flat, jdfit_muts, background_states)
plot_hist!(ax_B, z_df.av_z, z_df.av_dz)
plot_shift!(ax_C, dms_shift, Jtab, dfit_prot, cdiff_prot, cpair)

# Super-labels for A
fig.text(0.38, 0.01, L"Interaction parameters $J_{ij}$", ha="center", fontsize=16)
fig.text(0.05, 0.5, L"Mutation $\Delta f_i$", va="center", rotation="vertical", fontsize=16)

# Panel labels
axs_A_flat[1].text(-0.10, 1.02, "A", transform=axs_A_flat[1].transAxes, fontsize=18, fontweight="bold", va="bottom")
ax_B.text(-0.12, 1.02, "B", transform=ax_B.transAxes, fontsize=18, fontweight="bold", va="bottom")
ax_C.text(-0.12, 1.02, "C", transform=ax_C.transAxes, fontsize=18, fontweight="bold", va="bottom")

fig.tight_layout(rect=[0.04, 0.03, 1, 1])
fig.savefig("results/figures/fig_coup_valid.pdf", bbox_inches="tight")
close(fig)
