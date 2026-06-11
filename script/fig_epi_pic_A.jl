""" Script to generate panel A of figure 2 """

# Import packages
using DataFrames, CSV, PyPlot, Printf, Statistics

# Import fitness data
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)
prot = "S" # Protein of interest
dfit_prot = delta_fit[delta_fit.prot.==prot, :]
cnt_thr = 10.0 # Minimum count threshold for inclusion
dfit_prot = dfit_prot[(dfit_prot.exp_count1.>=cnt_thr).&(dfit_prot.exp_count2.>=cnt_thr), :]

# Specific clade pair to plot
clade_pair = ["21J", "21L"]
dfit_c1c2 = dfit_prot[(dfit_prot.clade1.==clade_pair[1]).&(dfit_prot.clade2.==clade_pair[2]), :]
mut = "P384L"
idx_mut = findall(dfit_c1c2.aa_mut .== mut)

# Range of scatter plot
fmax = maximum(vcat(dfit_c1c2.fit1, dfit_c1c2.fit2))
fmin = minimum(vcat(dfit_c1c2.fit1, dfit_c1c2.fit2))

# Correlation and number of points
r = cor(dfit_c1c2.fit1, dfit_c1c2.fit2)
n = length(dfit_c1c2.fit1)

# Make scatter plot
fig, ax = subplots(figsize=(5, 5))
ax.scatter(dfit_c1c2.fit1, dfit_c1c2.fit2, s=12, alpha=0.4, rasterized=true)
ax.scatter(dfit_c1c2.fit1[idx_mut], dfit_c1c2.fit2[idx_mut], s=12, color="red")
ax.scatter(dfit_c1c2.fit1[idx_mut], dfit_c1c2.fit2[idx_mut], s=120, facecolors="none", edgecolors="red", linewidths=1.5)
ax.annotate(mut, (dfit_c1c2.fit1[idx_mut][1], dfit_c1c2.fit2[idx_mut][1]),
    xytext=(8, 8), textcoords="offset points", fontsize=11, color="red")
diag = collect(fmin:0.2:fmax)
shift = 0.2
ax.set_xlim(fmin - shift, fmax + shift)
ax.set_ylim(fmin - shift, fmax + shift)
ax.plot(diag, diag, color="orange")
ax.set_xlabel(clade_pair[1] * " (Delta) fitness effect", fontweight="bold", fontsize=12)
ax.set_ylabel(clade_pair[2] * " (Omicron BA.2) fitness effect", fontweight="bold", fontsize=12)
fig.text(0.2, 0.85, "r = " * @sprintf("%.2f", r) * "\nn = " * string(n), fontsize=12)
fig.tight_layout()
fig.savefig("results/figures/fig_epi_pic_A.pdf")
close(fig)