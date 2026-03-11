""" Script to generate panel A of figure 7 of the SI"""

# Import packages
using DataFrames, CSV, JSON, JLD2, PyPlot, Statistics

# Function to compute the empirical cumulative density function
ecdf(x) = collect(1:1:length(x)) ./ length(x)

# Load raw and model curated average z-scores over clades
zscore_df = CSV.read("results/model_test_zscore.csv", DataFrame)
av_z = zscore_df.av_z
av_dz = zscore_df.av_dz

# Plot cumulative density of raw and model curated average z-scores
plot(sort(av_z), ecdf(av_z), label="Raw z")
plot(sort(av_dz), ecdf(av_dz), label="Predicted ϵ")
xscale(:log)
xlabel("Fitness discrepancy", fontsize=14)
ylabel("CDF", fontsize=14)
legend(fontsize=12)
tight_layout()
savefig("results/figures/si/fig_s7_A.pdf")
close("all")