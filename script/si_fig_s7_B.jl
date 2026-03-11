""" Script to generate panel B of figure 7  of the SI"""

# Import packages
using DataFrames, CSV, PyPlot, Statistics, JLD2

# Read clade-pair specific z-scores
@load "results/zmut_vec.jld2" z_vec
z_mut_vec = z_vec.z
zj_mut_vec = z_vec.zj

# Compute raw and model curated average z-scores over clades
av_z = sqrt.(mean.([(z_mut_vec[m]) .^ 2 for m in eachindex(z_mut_vec)]));
av_dz = sqrt.(mean.([(z_mut_vec[m] .- zj_mut_vec[m]) .^ 2 for m in eachindex(z_mut_vec)]));

# Indexes of predicted and not predicted mutations
idx_pred = findall(av_dz .< av_z) # Predicted mutations
idx_npred = findall(av_dz .>= av_z) # Not predicted mutations

# Compute correlations between raw and predicted z-scores
z_pred = z_mut_vec[idx_pred]
z_npred = z_mut_vec[idx_npred]
zj_pred = zj_mut_vec[idx_pred]
zj_npred = zj_mut_vec[idx_npred]
rho_vec_pred = [cor(z_pred[length.(z_pred).>2][n], zj_pred[length.(z_pred).>2][n]) for n in eachindex(idx_pred[length.(z_pred).>2])]
rho_vec_npred = [cor(z_npred[length.(z_npred).>2][n], zj_npred[length.(z_npred).>2][n]) for n in eachindex(idx_npred[length.(z_npred).>2])]

# Plot histograms of Pearson correlation coefficients
hist(rho_vec_pred, bins=20, density=true, alpha=0.7, edgecolor="black");
hist(rho_vec_npred, bins=20, density=true, alpha=0.5, edgecolor="black");
xlabel("Pearson correlation", fontsize=14)
ylabel("Density", fontsize=14)
legend(["ϵ < z", "ϵ ≥ z"], fontsize=13)
tight_layout()
savefig("results/figures/si/fig_s7_B.pdf")
close("all")