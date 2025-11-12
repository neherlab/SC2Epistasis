"""
    Script to generate panel A of figure 3
"""

# Import packages
using DataFrames, CSV, PyPlot, JSON

#################################---------------------------------

# Functions

function synonymous_dfit(; prots::Vector{String}=["E", "M", "N", "S"], cnt_thr::Float64=10.0)

    # Import by clade fitness estimates of amino acid mutations
    path_aamut_fit = "data/aamut_fit_by_clade/" # path to folder with fitness estimates
    files = readdir(path_aamut_fit, join=true) # list of files in the folder
    aamut_fit = vcat([DataFrame(CSV.File(f)) for f in files]...) # concatenate all dataframes
    rename!(aamut_fit, :cluster => :clade)

    # Read the json file containing clade founders mutations
    clade_muts = JSON.Parser.parsefile("data/clade_founders_neher.json")
    cmuts_roemer = JSON.Parser.parsefile("data/clade_founders_roemer.json")

    # Late clades not included in Neher's file
    pango_clades = ["XBB.1.5", "XBB.1.16", "CH.1.1", "XBB.1.9", "XBB.2.3", "EG.5.1", "XBB.1.5.70", "HK.3", "BA.2.86", "JN.1"]

    # Adding late clades to Neher's dictionary
    for c in pango_clades
        cld = cmuts_roemer[c]["nextstrainClade"]
        aamuts = cmuts_roemer[c]["aaSubstitutions"]
        nuc_muts = cmuts_roemer[c]["nucSubstitutions"]
        push!(clade_muts, cld => Dict{String,Any}("aa" => aamuts, "nuc" => nuc_muts))
    end

    # Deletions as background mismacthes between clades 
    clade_del = DataFrame(CSV.File("data/clade_del.csv"))

    # Exclude mutations of deleted residues
    all_del_clades = DataFrame(CSV.File("data/all_del_clades.csv"))
    aamut_fit = antijoin(aamut_fit, all_del_clades, on=[:clade, :gene, :aa_site])

    # Select synonymous mutations
    aamut_fit_syn = aamut_fit[aamut_fit.clade_founder_aa.==aamut_fit.mutant_aa, :]

    # Select clades to be included in the analysis
    clades = unique(aamut_fit_syn.clade)

    # Create dataframe for mutational fitness estimates in mismatching clades
    delta_fit_syn = SC2Epistasis.delta_fit(prots, clades, aamut_fit_syn, clade_muts, clade_del; cnt_thr=cnt_thr, exclude_stop=true, exclude_syn=false)

    return delta_fit_syn

end

#################################---------------------------------

# Script body

# Load dataframe of fitness discrepancies
delta_fit = CSV.read("results/delta_fit.csv", DataFrame)

cnt_thr = 20.0 # minimum count threshold for mutations to be included

delta_fit = delta_fit[(delta_fit.exp_count1.>=cnt_thr).&(delta_fit.exp_count2.>=cnt_thr), :]

# Stratify mutations as: nonsynonymous, synonymous, stop codon

# Stop-codon mutations
idx_sc = findall(x -> x[end] == '*', delta_fit.aa_mut)
dfit_stop = delta_fit[idx_sc, :]

# Non-synonymous mutations
idx_nsc = findall(x -> (x[end] != '*') && (x[1] != x[end]), delta_fit.aa_mut)
dfit_nsyn = delta_fit[idx_nsc, :]

# Fitness discrepancies for synonymous mutations
dfit_syn = synonymous_dfit(cnt_thr=cnt_thr)

# Compute z-scores
z_stop = SC2Epistasis.z_dfit(dfit_stop)
z_nsyn = SC2Epistasis.z_dfit(dfit_nsyn)
z_syn = SC2Epistasis.z_dfit(dfit_syn)

# Plot histogram
figure(figsize=(8, 4))
hist(z_nsyn, bins=50, alpha=0.6, density=true, label="Non-synonymous");
hist(z_syn, bins=50, alpha=0.5, density=true, label="Synonymous");
hist(z_stop, bins=50, alpha=0.4, density=true, label="Stop codon");
yscale("log")
legend(fontsize=12)
xlabel("Fitness z-score", fontsize=14);
ylabel("Density", fontsize=14);
tight_layout()
savefig("results/figures/fig_epi_pic_D.pdf", dpi=500);
close("all");