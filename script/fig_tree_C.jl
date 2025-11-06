""" Script to generate panel C of figure 1 """

# Import packages
using CSV, DataFrames, PyPlot, JSON

# Script body
#######################################---------------------------------------

# Compute number of mutations within clades

# Path to folder with mutation annotated sequence data by clade
dir_path = "data/muts_by_clade/"
# List of files
file_list = readdir(dir_path)

prot = "S"  # Protein of interest
mut_counts = Int64[]
clades = String[]

for f in file_list
    # Add clade to the list
    clade = split(f, "_")[end][1:3]
    push!(clades, clade)
    # Read the CSV file into a DataFrame
    df = CSV.read(joinpath(dir_path, f), DataFrame)

    df_muts = df[df.leaves_sharing_mutations.==1, :]
    muts = map(x -> split(x, ";"), df_muts.aa_mutations)

    mut_counts_clade = zeros(Int, length(muts))
    for (k, m) in enumerate(muts)
        for i in eachindex(m)
            gene, mut = split(m[i], ":")
            if gene == prot && mut[1] != mut[end]
                mut_counts_clade[k] += 1
            end
        end
    end
    mut_counts = vcat(mut_counts, mut_counts_clade)
end

# Compute mutations across clade founders

# Read the json file containing clade founders mutations
clade_muts = Dict(JSON.parsefile("data/clade_founders_neher.json"));
cmuts_roemer = Dict(JSON.parsefile("data/clade_founders_roemer.json"));

pango_clades = ["XBB.1.5", "XBB.1.16", "CH.1.1", "XBB.1.9", "XBB.2.3", "EG.5.1", "XBB.1.5.70", "HK.3", "BA.2.86", "JN.1"] # pango clades not in Neher's file

# Adding late clades to Neher's dictionary
for c in pango_clades
    cld = cmuts_roemer[c]["nextstrainClade"]
    aamuts = cmuts_roemer[c]["aaSubstitutions"]
    nuc_muts = cmuts_roemer[c]["nucSubstitutions"]
    push!(clade_muts, cld => Dict{String,Any}("aa" => aamuts, "nuc" => nuc_muts))
end

# Deletions as background mismacthes between clades 
clade_del = DataFrame(CSV.File("data/clade_del.csv"))

# Create dataframe with mismatches between clade founders on Spike protein
clade_diff = SC2Epistasis.clade_diff_list([prot], clades[2:end], clade_muts; clade_del=clade_del) # Exclude root for pairwise comparisons

# Dictionary to count number of mismatches between clade founders
mis_cf_dic = Dict{Vector{String},Int}()

for r in eachrow(clade_diff)
    cp = [r.clade1, r.clade2]
    if haskey(mis_cf_dic, cp)
        mis_cf_dic[cp] += 1
    elseif !haskey(mis_cf_dic, cp)
        push!(mis_cf_dic, cp => 1)
    end
end

# Plot histogram
hist(collect(values(mis_cf_dic)), bins=10, density=true, alpha=0.8, edgecolor="black")
hist(mut_counts, density=true, alpha=0.6, edgecolor="black")
yscale("log")
ylabel("Density", fontsize=14)
xlabel("Number of mutations", fontsize=14)
legend(["Between clade founders", "Within clades"], fontsize=12)
tight_layout()
savefig("results/figures/fig_tree_C.pdf")