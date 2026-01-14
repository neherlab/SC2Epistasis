""" Script to generate figure 8 of the manuscript """

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

function plot_shift(dms_shift_ba1_ba2_21j::DataFrame, dms_shift_ba2_xbb::DataFrame, Jtab::DataFrame,
    dfit_prot::DataFrame, cdiff_prot::DataFrame, clade_pairs::Vector{Tuple{S,S}};
    cnt_thr1::Vector{Float64}=[40.0, 20.0, 40.0],
    cnt_thr2::Vector{Float64}=[20.0, 20.0, 10.0]) where S<:AbstractString

    # Initialize figure
    fig, ax = subplots(1, 3, figsize=(13, 4))

    # Cycle over clade pairs
    for (i, cpair) in enumerate(clade_pairs)

        # Merge DMS shifts and fitness discrepancies
        shift_dfit = merge_dms_dfit(
            i == 1 || i == 2 ? dms_shift_ba1_ba2_21j : dms_shift_ba2_xbb,
            dfit_prot,
            cpair;
            cnt_thr1=cnt_thr1[i],
            cnt_thr2=cnt_thr2[i]
        )

        muts = shift_dfit.mutation # overlap mutations
        j_ddf = SC2Epistasis.mutcp_ddf(Jtab, muts, cdiff_prot, cpair[1], cpair[2]) # model predicted fitness discrepancies
        cld_shift = cpair[1]
        cld_comp = cpair[2]
        key1 = "score_" * cld_shift
        key2 = "score_" * cld_comp
        shift = shift_dfit[!, key1] .- shift_dfit[!, key2] # experimental shifts
        ρ = cor(shift, j_ddf) # correlation coefficient

        ax[i].scatter(j_ddf, shift, alpha=0.7, s=10)
        ax[i].set_title("$(cpair[1])-$(cpair[2])", fontsize=14)

        # Position text to avoid overlap with points
        xlim = ax[i].get_xlim()
        ylim = ax[i].get_ylim()
        x_range = xlim[2] - xlim[1]
        y_range = ylim[2] - ylim[1]

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
            x_text = xlim[1] + pos[1] * x_range
            y_text = ylim[1] + pos[2] * y_range
            # Count nearby points (within 10% of range)
            overlap = sum((abs.(j_ddf .- x_text) .< 0.1 * x_range) .& (abs.(shift .- y_text) .< 0.1 * y_range))
            if overlap < min_overlap
                min_overlap = overlap
                best_pos = pos
            end
        end

        ax[i].text(best_pos[1], best_pos[2], "ρ = " * string(round(ρ, digits=2)),
            fontsize=12, transform=ax[i].transAxes,
            verticalalignment=best_pos[3], horizontalalignment=best_pos[4])

    end

    fig.supxlabel("ΔΔϕ", fontsize=14, ha="center")
    fig.supylabel("Experimental fitness shift", fontsize=14, va="center")
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
dfit_prot = delta_fit[delta_fit.prot.==prot, :]
cdiff_prot = clade_diff[clade_diff.prot.==prot, :]

# Clade pairs to analyze
clade_pairs = [("21J", "21K"), ("21L", "21K"), ("23A", "21L")]

# Predicted counts vectors
cnt_thr1 = [40.0, 20.0, 20.0]
cnt_thr2 = [20.0, 20.0, 15.0]

# Initialize figure
fig, ax = plot_shift(dms_shift_ba1_ba2_21j, dms_shift_ba2_xbb, Jtab, dfit_prot, cdiff_prot, clade_pairs; cnt_thr1=cnt_thr1, cnt_thr2=cnt_thr2)
fig.savefig("results/figures/fig_dms_comp.pdf")
close(fig)
