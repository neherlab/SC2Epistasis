""" Compute correlation between DMS shifts and model predictions on a domain-wise basis """

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

# Merge DMS shifts and fitness discrepancies
shift_dfit = merge_dms_dfit(
    dms_shift_ba1_ba2_21j,
    dfit_prot,
    cpair;
    cnt_thr1=cnt_thr1,
    cnt_thr2=cnt_thr2
)

# Boundaries of Spike domains
dom_edges = [14:305, 319:541, 542:690] # NTD, RBD, CTD1+CTD2
regions = ["NTD", "RBD", "CTD1+CTD2", "Other"]

shift_dfit[!, :site] = map(x -> parse(Int, x[2:end-1]), shift_dfit.mutation) # extract site number from mutation string
idx_reg = [map(x -> x in dom_edges[i], shift_dfit.site) for i in eachindex(dom_edges)] # indices of mutations in each domain

idx_all_reg = zeros(Bool, size(shift_dfit, 1)) # indices of mutations in all large domains
for i in eachindex(idx_reg)
    idx_all_reg .+= idx_reg[i]
end

idx_rest = .!idx_all_reg # indices of mutations in other regions

# Combine indices for all regions (large domains + rest)
idx_cor = vcat(idx_reg, [idx_rest])

rho = zeros(Float64, length(idx_cor)) # vector to store correlation coefficients for each region
for i in eachindex(idx_cor)
    shift = shift_dfit[idx_cor[i], :score_21J] .- shift_dfit[idx_cor[i], :score_21K]
    j_ddf = SC2Epistasis.mutcp_ddf(Jtab, shift_dfit[idx_cor[i], :mutation], cdiff_prot, cpair[1], cpair[2])
    rho[i] = cor(shift, j_ddf)
end

print("Correlation coefficients for each region: " * "NTD-" * string(rho[1]) * ", RBD-" * string(rho[2]) * ", CTD1+CTD2-" * string(rho[3]) * ", Other-" * string(rho[4]) * "\n")