"""Script to generate panel C of figure 5"""

# Import packages
using DataFrames, CSV, PyPlot

symbol_size = 20

# Plotting function
function plot_int_num(dist::Vector{Int}, int_num::Vector{Vector{Float64}}, coeffs::Vector{Float64})

    @assert length(int_num) == length(coeffs)

    scatter(dist, int_num[1], s=symbol_size, alpha=0.8, label=L"\Delta\Delta\phi_{\mathrm{thr}}=1,\,\propto 3")
    plot(dist, coeffs[1] .* dist, color="C0", alpha=0.7, label="_nolegend_")

    scatter(dist, int_num[2], s=symbol_size, alpha=0.8, label=L"\Delta\Delta\phi_{\mathrm{thr}}=1.5,\,\propto 1.5")
    plot(dist, coeffs[2] .* dist, color="C1", alpha=0.7, label="_nolegend_")

    scatter(dist, int_num[3], s=symbol_size, alpha=0.8, label=L"\Delta\Delta\phi_{\mathrm{thr}}=2,\,\propto 0.7")
    plot(dist, coeffs[3] .* dist, color="C2", alpha=0.7, label="_nolegend_")

    xlabel("Number of mismatches", fontsize=12)
    ylabel("Effective number of interactions", fontsize=12)

    legend(fontsize=13)

    tight_layout()

end

# Script body

# Load coupling parameters
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)
# Load fitness discrepancies
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
# Load clade pairs differences
clade_diff = CSV.read("results/clade_diff.csv", DataFrame)

# Filter for the protein of interest
prot = "S"
dfit_prot = delta_fit[(delta_fit.prot.==prot), :]
cdiff_prot = clade_diff[(clade_diff.prot.==prot), :]

# Retain major clades for the analysis
clades = ["20I", "21J", "21K", "21L", "22B", "22E", "23A", "23I"]

# Corresponding clade pairs
cpairs = NTuple{2,String}[]
for k1 in 1:length(clades)-1
    for k2 in k1+1:length(clades)
        push!(cpairs, (clades[k1], clades[k2]))
    end
end

# Subset of shared mutations
muts = dfit_prot[(dfit_prot.clade1.==cpairs[1][1]).&(dfit_prot.clade2.==cpairs[1][2]), :aa_mut]
for c in eachindex(cpairs)[2:end]
    m = dfit_prot[(dfit_prot.clade1.==cpairs[c][1]).&(dfit_prot.clade2.==cpairs[c][2]), :aa_mut]
    global muts = intersect(muts, m)
end

# Populate ΔJ vectors for each clade pair
dj_vec = [zeros(length(muts)) for c in eachindex(cpairs)]
# Populate mismatches vectors for each clade pair
dist = zeros(Int64, length(cpairs))

# Loop over clade pairs
for c in eachindex(cpairs)
    cdiff_c1c2 = cdiff_prot[(cdiff_prot.clade1.==cpairs[c][1]).&(cdiff_prot.clade2.==cpairs[c][2]), :]
    n_mis = size(cdiff_c1c2, 1)
    dist[c] = n_mis
    dj_vec[c] = SC2Epistasis.mutcp_ddf(Jtab, muts, cdiff_prot, cpairs[c][1], cpairs[c][2])
end

# Threshold for significant |ΔJ|
thr_vec = [1.0, 1.5, 2.0]

larger_thresh = [zeros(length(dist)) for n in eachindex(thr_vec)]
for n in eachindex(thr_vec)
    larger_thresh[n] .= map(x -> sum(abs.(x) .> thr_vec[n]), dj_vec)
end

# Best fit coefficients without intercept
coeffs = [sum(larger_thresh[n] .* dist) / sum(dist .^ 2) for n in eachindex(thr_vec)]

# Make plot
plot_int_num(dist, larger_thresh, coeffs)
savefig("results/figures/fig_coup_visual_C.pdf")
close("all")