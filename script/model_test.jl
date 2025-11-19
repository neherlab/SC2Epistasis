""" Compute and save model prediction on test clades for fitness z-scores """

# Import packages
using DataFrames, CSV, JSON, JLD2

###########################################---------------------------------------

# Create data frame with mutational fitness discrepancies including JN.1

# Import by clade fitness estimates of amino acid mutations
path_aamut_fit = "data/aamut_fit_by_clade/" # path to folder with fitness estimates
files = readdir(path_aamut_fit, join=true) # list of files in the folder
push!(files, "data/jn1_fitness/24A_aamut_fitness.csv.gz") # add JN.1 fitness estimates
aamut_fit = vcat([DataFrame(CSV.File(f)) for f in files]...); # concatenate all dataframes
rename!(aamut_fit, :cluster => :clade)
prot = "S"
aamut_fit_prot = aamut_fit[aamut_fit.gene.==prot, :]

# Import mismacthes between clade founders

# Read the json file containing clade founders mutations
clade_muts = Dict(JSON.parsefile("data/clade_founders_neher.json"));
cmuts_roemer = Dict(JSON.parsefile("data/clade_founders_roemer.json"));

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

# Dictionary from Nextstrain to pango clades
next_pango_dic = Dict(JSON.parsefile("data/next_pango.json"));

# Clades list including JN.1
clades = unique(aamut_fit.clade)

# List of deletions in clade founders
all_del_clades = DataFrame(clade=String[], gene=String[], aa_site=Int[])
for c in clades

    png_cld = next_pango_dic[c]
    aa_del = cmuts_roemer[png_cld]["aaDeletions"]
    if aa_del[1] == ""
        continue
    end
    for aa in aa_del
        p, mut = split(aa, ":")
        if p == prot
            s = parse(Int, mut[2:end-1])
            push!(all_del_clades, [c, p, s])
        else
            continue
        end
    end

end

# Exclude mutations of deleted residues
aamut_fit = antijoin(aamut_fit_prot, all_del_clades, on=[:clade, :gene, :aa_site])

# Counts threshold
cnt_thr1 = 10.0 # minimum count threshold for mutations to be included
cnt_thr2 = 5.0 # minimum count threshold for mutations to be included in JN.1

aamut_fit = aamut_fit[((aamut_fit.clade.!="24A").&(aamut_fit.predicted_count.>=cnt_thr1)).|((aamut_fit.clade.=="24A").&(aamut_fit.predicted_count.>=cnt_thr2)), :]
idx_nsc = findall(aamut_fit.mutant_aa .!= "*") # exclude stop codons
aamut_fit = aamut_fit[idx_nsc, :]
idx_no_syn = findall(x -> x[1] != x[end], aamut_fit.aa_mutation) # exclude synonymous mutations
aamut_fit = aamut_fit[idx_no_syn, :]

dfit = SC2Epistasis.delta_fit(prot, clades, aamut_fit, clade_muts, clade_del)

# Create dataframe with mismatches between clade founders in structural proteins including JN.1
clade_diff = SC2Epistasis.clade_diff_list([prot], clades, clade_muts; clade_del=clade_del)

###########################################---------------------------------------

# Load couplings inferred from train set
Jtab = CSV.read("results/jcoup_train.csv", DataFrame)

# Train-test split
train_clades = clades[1:26]
test_clades = clades[27:end]
comp_clust = copy(train_clades)

# Train dataset
dfit_train = filter(r -> r.clade1 in train_clades && r.clade2 in train_clades, dfit)
# Test dataset
dfit_test = filter(r -> (r.clade1 in comp_clust && r.clade2 in test_clades) || (r.clade1 in test_clades && r.clade2 in comp_clust) || (r.clade1 in test_clades && r.clade2 in test_clades), dfit)

# Mutations in train-test sets
muts_train = unique(dfit_train.aa_mut)
idx_sort_mut = sortperm(map(x -> parse(Int64, x[2:end-1]), muts_train)); # Sort mutations
muts_train = muts_train[idx_sort_mut];

# Select mutations which are both in train and test set
dfit_test = filter(r -> r.aa_mut in muts_train, dfit_test)
idx_ovrlp = findall(x -> x in dfit_test.aa_mut, muts_train)
muts_test = muts_train[idx_ovrlp]

# Compute vectors of raw and model curated z-scores for mutations in the test set
z_mut_vec, zj_mut_vec = SC2Epistasis.zscore_mut(muts_test, Jtab, dfit_train, dfit_test, clade_diff)

# Excluding empty z-scores
idx_no_empty = findall(x -> !isempty(x), zj_mut_vec)
muts_nemp = muts_test[idx_no_empty]
z_mut_vec = z_mut_vec[idx_no_empty]
zj_mut_vec = zj_mut_vec[idx_no_empty]

# Write clade-pair specific z-scores to JLD file
z_vec = (z=z_mut_vec, zj=zj_mut_vec)
@save "results/zmut_vec.jld2" z_vec

# Compute raw and model curated average z-scores over clades
av_z = sqrt.(mean.([(z_mut_vec[m]) .^ 2 for m in eachindex(z_mut_vec)]));
av_dz = sqrt.(mean.([(z_mut_vec[m] .- zj_mut_vec[m]) .^ 2 for m in eachindex(z_mut_vec)]));

# Create dataframe with: mutation, average raw z-score, average model curated z-score
zscore_df = DataFrame(aa_mut=muts_nemp, av_z=av_z, av_dz=av_dz)
CSV.write("results/model_test_zscore.csv", zscore_df)