""" Script to generate panel A of figure 7 """

# Import packages
using DataFrames, CSV, PyPlot

# Read raw and model curated z-scores for test mutations
z_df = CSV.read("results/model_test_zscore.csv", DataFrame)
av_z = z_df.av_z
av_dz = z_df.av_dz

# Make histograms
hist(av_z, bins=30, density=true, alpha=0.7, edgecolor="black")
hist(av_dz, bins=30, density=true, alpha=0.5, edgecolor="black")
xlabel("Fitness discrepancy z-score", fontsize=14)
ylabel("Density", fontsize=14)
legend(["Raw discrepancies z", "Model predicted ϵ"], fontsize=13)
tight_layout()
savefig("results/figures/fig_coup_valid_B.pdf")
close("all")