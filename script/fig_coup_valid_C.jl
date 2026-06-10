""" Script to generate panel C of figure 6 """

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
    df_merge = df_merge[(df_merge.exp_count1.>=cnt_thr1).&(df_merge.exp_count2.>=cnt_thr2), :]

    return df_merge

end

# Plotting function

function plot_shift(dms_shift::DataFrame, Jtab::DataFrame,
    dfit_prot::DataFrame, cdiff_prot::DataFrame, clade_pair::Tuple{S,S};
    cnt_thr1::Float64=40.0,
    cnt_thr2::Float64=20.0,
    doms_label::Vector{S1}=["NTD", "RBD", "CTD1", "CTD2", "FP", "HR1", "CH", "CD", "HR2"],
    doms_edges::Vector{UnitRange{Int}}=[14:305, 319:541, 542:590, 591:690, 788:806, 910:984, 985:1034, 1035:1067, 1163:1212],
    doms_col::Vector{S1}=["cyan", "blue", "orange", "yellow", "red", "green", "lightgreen", "hotpink", "violet"]) where {S<:AbstractString,S1<:AbstractString}

    # Merge DMS shifts and fitness discrepancies
    shift_dfit = merge_dms_dfit(
        dms_shift,
        dfit_prot,
        clade_pair;
        cnt_thr1=cnt_thr1,
        cnt_thr2=cnt_thr2
    )

    muts = shift_dfit.mutation # overlap mutations
    res_index_vec = parse.(Int, map(x -> x[2:end-1], muts))
    j_ddf = SC2Epistasis.mutcp_ddf(Jtab, muts, cdiff_prot, clade_pair[1], clade_pair[2]) # model predicted fitness discrepancies
    cld_shift = clade_pair[1]
    cld_comp = clade_pair[2]
    key1 = "score_" * cld_shift
    key2 = "score_" * cld_comp
    shift = shift_dfit[!, key1] .- shift_dfit[!, key2] # experimental shifts
    ρ = cor(shift, j_ddf) # correlation coefficient

    domains = [(doms_edges[n][1], doms_edges[n][end], doms_label[n], doms_col[n]) for n in eachindex(doms_label)]

    # Color each point by its domain assignment
    domain_masks = [(res_index_vec .>= start) .& (res_index_vec .<= stop) for (start, stop, _, _) in domains]
    for (start, stop, label, color) in domains
        in_domain = (res_index_vec .>= start) .& (res_index_vec .<= stop)
        if any(in_domain)
            scatter(j_ddf[in_domain], shift[in_domain], alpha=0.7, s=10, color=color, label=label)
        end
    end

    other_domain = .!foldl((acc, mask) -> acc .| mask, domain_masks)
    if any(other_domain)
        scatter(j_ddf[other_domain], shift[other_domain], alpha=0.7, s=10, color="lightgrey", label="Other")
    end
    #title("$(clade_pair[1])-$(clade_pair[2])", fontsize=14)

    # Get current axis limits
    xlim_vals = xlim()
    ylim_vals = ylim()
    x_range = xlim_vals[2] - xlim_vals[1]
    y_range = ylim_vals[2] - ylim_vals[1]

    # Try different positions: top-left, top-right, bottom-left, bottom-right
    positions = [
        (0.05, 0.95, "top", "left"),
        (0.95, 0.95, "top", "right"),
        (0.05, 0.05, "bottom", "left"),
        (0.95, 0.05, "bottom", "right")
    ]

    # Find position with least overlap
    best_pos = positions[1]
    min_overlap = Inf

    for pos in positions
        x_text = xlim_vals[1] + pos[1] * x_range
        y_text = ylim_vals[1] + pos[2] * y_range
        # Count nearby points (within 10% of range)
        overlap = sum((abs.(j_ddf .- x_text) .< 0.1 * x_range) .& (abs.(shift .- y_text) .< 0.1 * y_range))
        if overlap < min_overlap
            min_overlap = overlap
            best_pos = pos
        end
    end

    # Get current axes
    ax = gca()
    text(best_pos[1], best_pos[2], "ρ = " * string(round(ρ, digits=2)),
        fontsize=12, transform=ax.transAxes,
        verticalalignment=best_pos[3], horizontalalignment=best_pos[4])

    # Legend for domain-color mapping
    legend_handles = [mpatches.Patch(color=color, label=label) for (_, _, label, color) in domains]
    legend_entries = [label for (_, _, label, _) in domains]
    if any(other_domain)
        push!(legend_handles, mpatches.Patch(color="lightgrey", label="Other"))
        push!(legend_entries, "Other")
    end
    legend(legend_handles, legend_entries, title="Domain", fontsize=10, title_fontsize=11,
        loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0, frameon=true)

    xlabel("ΔΔϕ", fontsize=14)
    ylabel("Experimental fitness shift", fontsize=14)
    tight_layout()

end

###########################################------------------------------------

# Script body

# Load dataframes with estimated shifts from DMS experiments
dms_shift_ba1_ba2_21j = CSV.read("data/dms_shift/dms_shift_ba1_ba2_21j.csv", DataFrame)

# Load inferred coupling parameters
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)

# Load dataframe with fitness discrepancies
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
rename!(delta_fit, :aa_mut => :mutation) # rename column for consistency

# Load dataframe with clade founder mismatches
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Select Spike protein
prot = "S"
dfit_prot = delta_fit[delta_fit.prot.==prot, :]
cdiff_prot = clade_diff[clade_diff.prot.==prot, :]

# Clade pairs to analyze
cpair = ("21J", "21K")

# Predicted counts vectors
cnt_thr1 = 40.0
cnt_thr2 = 20.0

# Initialize figure
plot_shift(dms_shift_ba1_ba2_21j, Jtab, dfit_prot, cdiff_prot, cpair; cnt_thr1=cnt_thr1, cnt_thr2=cnt_thr2)
savefig("results/figures/fig_coup_valid_C.pdf")
close("all")
