""" Script to generate figure 8 of the SI """

# Import packages
using DataFrames, CSV, PyPlot, Statistics

# Functions

# Merge clade-pair specific and epistatic shifts dataframes
function merge_dms_dfit(dms_shift::DataFrame, dfit::DataFrame, clade_pair::Tuple{S,S};
    cnt_thr1::Float64=10.0,
    cnt_thr2::Float64=10.0) where S<:AbstractString

    dfit_c1c2 = filter(row -> row.clade1 == clade_pair[1] && row.clade2 == clade_pair[2], dfit)
    if isempty(dfit_c1c2)
        dfit_c1c2 = filter(row -> row.clade1 == clade_pair[2] && row.clade2 == clade_pair[1], dfit)
    end

    # Merge dataframes
    df_merge = innerjoin(dms_shift, dfit_c1c2, on=:mutation)
    df_merge = df_merge[(df_merge.exp_count1 .>= cnt_thr1) .& (df_merge.exp_count2 .>= cnt_thr2), :]

    return df_merge

end

# Plotting function
function plot_shift(dms_shift_ba1_ba2_21j::DataFrame, dms_shift_ba2_xbb::DataFrame, Jtab::DataFrame,
    dfit_prot::DataFrame, cdiff_prot::DataFrame, clade_pairs::Vector{Tuple{S,S}};
    ncols::Int=2,
    cnt_thr1::Vector{Float64}=[20.0, 20.0],
    cnt_thr2::Vector{Float64}=[20.0, 15.0],
    doms_label::Vector{S}=["NTD", "RBD"],
    doms_edges::Vector{UnitRange{Int}}=[14:305, 319:541],
    doms_col::Vector{S}=["orange", "blue"],
    doms_marker::Vector{S}=["s", "^"]) where S<:AbstractString

    @assert length(clade_pairs) == ncols

    # Initialize figure
    fig, ax = subplots(1, ncols, figsize=(ncols * 4 + 1, 4))

    # Cycle over clade pairs
    for (i, cpair) in enumerate(clade_pairs)

        # Merge DMS shifts and fitness discrepancies
        shift_dfit = merge_dms_dfit(
            cpair[2] == "21K" ? dms_shift_ba1_ba2_21j : dms_shift_ba2_xbb,
            dfit_prot,
            cpair;
            cnt_thr1=cnt_thr1[i],
            cnt_thr2=cnt_thr2[i]
        )

        muts = shift_dfit.mutation # overlap mutations
        res_index_vec = parse.(Int, map(x -> x[2:(end-1)], muts)) # residue indices
        j_ddf = SC2Epistasis.mutcp_ddf(Jtab, muts, cdiff_prot, cpair[1], cpair[2]) # model predicted fitness discrepancies
        cld_shift = cpair[1]
        cld_comp = cpair[2]
        key1 = "score_" * cld_shift
        key2 = "score_" * cld_comp
        shift = shift_dfit[!, key1] .- shift_dfit[!, key2] # experimental shifts
        ρ = cor(shift, j_ddf) # correlation coefficient

        if length(clade_pairs) == 1
            ax = [ax] # ensure ax is an array for single subplot case
        end

        domains = [(doms_edges[n][1], doms_edges[n][end], doms_label[n], doms_col[n]) for n in eachindex(doms_label)]
        domain_masks = [(res_index_vec .>= start) .& (res_index_vec .<= stop) for (start, stop, _, _) in domains]
        other_domain = .!foldl((acc, mask) -> acc .| mask, domain_masks)

        if any(other_domain)
            ρ = cor(shift[other_domain], j_ddf[other_domain])
            ax[i].scatter(j_ddf[other_domain], shift[other_domain], alpha=0.7, s=15, color="lightgreen", marker="o", label="Other" * " (ρ=" * string(round(ρ, digits=2)) * ")")
        end
        for (n, (start, stop, label, color)) in enumerate(domains)
            in_domain = (res_index_vec .>= start) .& (res_index_vec .<= stop)
            ρ = cor(shift[in_domain], j_ddf[in_domain])
            if any(in_domain)
                ax[i].scatter(j_ddf[in_domain], shift[in_domain], alpha=0.7, s=15, color=color, marker=doms_marker[n], label=label * " (ρ=" * string(round(ρ, digits=2)) * ")")
            end
        end

        ax[i].legend(loc="upper left", fontsize=11, frameon=true)
        ax[i].set_title("$(cpair[1])-$(cpair[2])", fontsize=14)

    end

    fig.supylabel("Experimental fitness shift", fontsize=14, va="center")
    fig.supxlabel("Model prediction ΔΔϕ", fontsize=14, ha="center")
    fig.tight_layout()

    return fig, ax

end

###########################################------------------------------------

# Script body

# Load dataframes with estimated shifts from DMS experiments
dms_shift_ba1_ba2_21j = CSV.read("data/dms_shift/dms_shift_ba1_ba2_21j.csv", DataFrame)
dms_shift_ba2_xbb = CSV.read("data/dms_shift/dms_shift_ba2_xbb.csv", DataFrame)

rename!(dms_shift_ba2_xbb, [:shift_XBB => :shift_23A, :predicted_func_score_BA2 => :score_21L, :predicted_func_score_XBB => :score_23A]) # switch to nextstrain convention

# Load inferred coupling parameters
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)

# Load dataframe with fitness discrepancies
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
rename!(delta_fit, :aa_mut => :mutation) # rename column for consistency

# Load dataframe with clade founder mismatches
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Select Spike protein
prot = "S"
dfit_prot = delta_fit[delta_fit.prot .== prot, :]
cdiff_prot = clade_diff[clade_diff.prot .== prot, :]

# Clade pairs to analyze
cpairs = [("21L", "21K"), ("23A", "21L")]

# Predicted counts vectors
cnt_thr1 = [20.0, 20.0]
cnt_thr2 = [20.0, 15.0]

# Make plot
fig, ax = plot_shift(dms_shift_ba1_ba2_21j, dms_shift_ba2_xbb, Jtab, dfit_prot, cdiff_prot, cpairs;
    cnt_thr1=cnt_thr1, cnt_thr2=cnt_thr2)

# Save figure
savefig("results/figures/si/fig_s8.pdf")
close(fig)