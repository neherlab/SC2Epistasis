"""
    Script to generate panel C of figure 2
"""

# Run preliminary scripts
include("z_struct_prots.jl")
include("z_clades_regions.jl")


# Make plot
fig, ax = subplots(2, 1, figsize=(6, 9), sharex=true)

# Top histogram
ax[1].hist(z_clades, bins=40, density=true, alpha=0.8, label="Clades", edgecolor="black");
ax[1].hist(z_country, bins=20, density=true, alpha=0.6, label="USA-UK", edgecolor="black");
ax[1].legend(fontsize=13);
ax[1].tick_params(axis="y", labelsize=14);
ax[1].set_yscale("log");

# Bottom histogram
ax[2].hist(z_nsyn, bins=50, alpha=0.6, density=true, label="Non-synonymous");
ax[2].hist(z_syn, bins=50, alpha=0.5, density=true, label="Synonymous");
ax[2].hist(z_stop, bins=50, alpha=0.4, density=true, label="Stop codon");
ax[2].set_xlabel("Fitness z-score", fontsize=16);
ax[2].legend(fontsize=13);
ax[2].tick_params(axis="both", labelsize=14);
ax[2].set_yscale("log");

fig.tight_layout(rect=[0.05, 0, 1, 1]);

# Add common y-axis label
fig.text(0.02, 0.5, "Density", va="center", rotation="vertical", fontsize=16);

# Save figure
fig.savefig("results/figures/fig_epi_pic_C.pdf");
close(fig)
