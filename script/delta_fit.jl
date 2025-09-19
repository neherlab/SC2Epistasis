""" 
    Script to create a dataframe with mutational fitness estimates
    for all pairs of clades that differ for at least one background
    mismatch in the structural proteins
"""

# Import packages
using DataFrames, CSV, JSON, HypothesisTests

# Import by clade fitness estimates of amino acid mutations
path_aamut_fit = "data/aamut_fit_by_clade/" # path to folder with fitness estimates
files = readdir(path_aamut_fit, join=true) # list of files in the folder
aamut_fit = vcat([DataFrame(CSV.File(f)) for f in files]...); # concatenate all dataframes
rename!(aamut_fit, :cluster => :clade)

# Read the json file containing clade founders mutations
clade_muts = JSON.Parser.parsefile("data/clade_founders_neher.json");
cmuts_roemer = JSON.Parser.parsefile("data/clade_founders_roemer.json");

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

# Select structural proteins
prots = ["E", "M", "N", "S"]

# Dictionary from Nexst-strain to pango clades
next_pango_dic = JSON.Parser.parsefile("data/next_pango.json");

# List of deletions in clade founders
all_del_clades = DataFrame(CSV.File("data/all_del_clades.csv"))

# Exclude mutations of deleted residues
aamut_fit = antijoin(aamut_fit, all_del_clades, on=[:clade, :gene, :aa_site])

# Select clades to be included in the analysis
clades = unique(aamut_fit.clade)

# Create dataframe for mutational fitness estimates in mismatching clades
cnt_thr = 10.0 # minimum count threshold for mutations to be included
delta_fit = SC2Epistasis.delta_fit(prots, clades, aamut_fit, clade_muts, clade_del; cnt_thr=cnt_thr, exclude_stop=false, exclude_syn=true)

# Create dataframe with mismatches between clade founders in structural proteins
clade_diff = SC2Epistasis.clade_diff_list(prots, clades, clade_muts; clade_del=clade_del)

# Write dataframes to CSV file
CSV.write("results/delta_fit.csv", delta_fit)
CSV.write("results/clade_diff.csv", clade_diff)