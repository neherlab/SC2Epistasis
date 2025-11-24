""" Functions to compute predictions from inferred models """

#####################################---------------------------------------

# Determining mismatching residues and amino acids observed in the train and test sets

# Collect mismatching sites and amino acids for a set of mutations and associated clades
function mismatch(dfit::DataFrame, clade_diff::DataFrame)

    muts = unique(dfit.aa_mut)
    mis_dic = [Dict{Int,Vector{String}}() for n in eachindex(muts)]
    for n in eachindex(muts)
        mut = muts[n]
        mis_dic[n] = mismatch(mut, dfit, clade_diff)
    end

    return mis_dic

end

# Collect mismatching sites and amino acids for a given mutation and associated clades
function mismatch(mut::S, dfit::DataFrame, clade_diff::DataFrame) where {S<:AbstractString}

    mis_dic = Dict{Int,Vector{String}}()
    dfit_mut = dfit[dfit.aa_mut.==mut, :]

    clade_pairs = SC2Epistasis.unique_clade_pairs(dfit_mut)

    for cp in clade_pairs

        c1 = cp[1]
        c2 = cp[2]

        cdiff_c1c2 = clade_diff[(clade_diff.clade1.==c1).&(clade_diff.clade2.==c2), :]
        mismatch(mis_dic, cdiff_c1c2)

    end

    return mis_dic

end

# Populate dictionary `mis_dic` with mismatching sites and amino acids from
# reference dataframe `clade_diff` containing mismatches between clade founders
function mismatch(mis_dic::Dict{Int,Vector{String}}, clade_diff::DataFrame)

    for r in eachrow(clade_diff)
        site = r.site
        aa1 = r.aa_c1
        aa2 = r.aa_c2

        if !haskey(mis_dic, site)
            push!(mis_dic, site => [aa1, aa2])
        elseif haskey(mis_dic, site)
            if !(aa1 in mis_dic[site])
                push!(mis_dic[site], aa1)
            end
            if !(aa2 in mis_dic[site])
                push!(mis_dic[site], aa2)
            end
        end

    end

end

# Given a dataframe containing mismacthes between clade founders, returns a dictionary with:
# - Keys: sites with mismatches
# - Values: list of amino acids observed at the site across clade founders
function mismatch(clade_diff::DataFrame)

    mis_dic = Dict{Int,Vector{String}}()

    for r in eachrow(clade_diff)
        site = r.site
        aa1 = r.aa_c1
        aa2 = r.aa_c2

        if !haskey(mis_dic, site)
            push!(mis_dic, site => [aa1, aa2])
        elseif haskey(mis_dic, site)
            if !(aa1 in mis_dic[site])
                push!(mis_dic[site], aa1)
            end
            if !(aa2 in mis_dic[site])
                push!(mis_dic[site], aa2)
            end
        end

    end

    return mis_dic

end

# Taking as input two dictionaries of mismatches between clade founders 
# `mis_dic_c1c2` refers to the mismatches between a specific pair of clades
# that belong to the test set, while `mis_dic` refers to the mismatches between
# all clade founders in the training set.
# The function returns the list of sites where the amino acids corresponding to
# mismatches in the test set are a subset of those in the training set.
function retain_sites(mis_dic_c1c2::Dict{Int,Vector{S}}, mis_dic::Dict{Int,Vector{S}}) where {S<:AbstractString}

    sites = sort(collect(keys(mis_dic_c1c2)))
    idx_sites = zeros(Bool, length(sites))

    for k in eachindex(sites)
        s = sites[k]
        if haskey(mis_dic, s)
            idx_sites[k] = issubset(Set(mis_dic_c1c2[s]), Set(mis_dic[s]))
        else
            continue
        end
    end

    return sites[idx_sites]

end

####################################---------------------------------------

# Compute predicted fitness discrepancies and z-scores

# Compute the raw and model curated z-scores for the list of mutations in the test set
function zscore_mut(muts::Vector{S1}, Jtab::DataFrame, dfit_train::DataFrame, dfit_test::DataFrame, cdiff::DataFrame) where {S1<:AbstractString}

    # Initialize vector of z-scores
    z_mut_vec = [Float64[] for n in eachindex(muts)]
    zj_mut_vec = [Float64[] for n in eachindex(muts)]

    for n in eachindex(muts)
        mut = muts[n]
        mis_mut_train = mismatch(mut, dfit_train, cdiff)
        z_mut_vec[n], zj_mut_vec[n] = zscore_mut(mut, mis_mut_train, Jtab, dfit_test, cdiff)
    end

    return z_mut_vec, zj_mut_vec

end

# Compute raw and model curated z-scores of fitness discrepancies for a given mutation
# Input:
# - mut: mutation
# - mis_mut_train: dictionary of mismatches observed in the training set
# - Jtab: dataframe of inferred couplings
# - dfit: dataframe of fitness discrepancies in the test set
# - cdiff: dataframe of mismatches between clade founders
# Output:
# - z_mut_vec: vector of raw z-scores of fitness discrepancies averaged over clade pairs
# - zj_mut_vec: vector of model curated z-scores of fitness discrepancies averaged over clade pairs
function zscore_mut(mut::S1, mis_mut_train::Dict{Int,Vector{S2}}, Jtab::DataFrame, dfit::DataFrame, cdiff::DataFrame) where {S1<:AbstractString,S2<:AbstractString}

    Jmut = Jtab[(Jtab.σᵢ_wt.==string(mut[1])).&(Jtab.i.==parse(Int, mut[2:end-1])).&(Jtab.σᵢ.==string(mut[end])), :]
    dfit_mut = dfit[dfit.aa_mut.==mut, :]
    cpairs = SC2Epistasis.unique_clade_pairs(dfit_mut)
    idx_ret = zeros(Bool, length(cpairs))

    zj_mut_vec = Float64[]

    for (n, cp) in enumerate(cpairs)
        c1 = cp[1]
        c2 = cp[2]
        cdiff_c1c2 = cdiff[(cdiff.clade1.==c1).&(cdiff.clade2.==c2), :]
        mis_mut = mismatch(cdiff_c1c2)
        sites = retain_sites(mis_mut, mis_mut_train)
        idx_ret[n] = !isempty(sites) ? true : false
        if idx_ret[n]
            Jmut_sites = filter(r -> r.j in sites, Jmut)
            if all(Jmut_sites.J .== 0.0)
                idx_ret[n] = false
                continue
            else
                push!(zj_mut_vec, SC2Epistasis.mutcp_ddf(Jmut_sites, mut, cdiff, c1, c2))
            end
        end
    end

    dfit_mut = dfit_mut[idx_ret, :]
    z_mut_vec = SC2Epistasis.z_dfit(dfit_mut)
    s = sqrt.(dfit_mut.std_fit1 .^ 2 .+ dfit_mut.std_fit2 .^ 2)
    zj_mut_vec = zj_mut_vec ./ s
    @assert length(z_mut_vec) == length(zj_mut_vec)

    return z_mut_vec, zj_mut_vec

end

# Function to compute the model predicted z-score for a specific mutation and clade pair
function mutcp_z(J::DataFrame, mut::S, dfit::DataFrame, cdiff::DataFrame) where {S<:AbstractString}

    N = size(dfit, 1)
    ddf = zeros(Float64, N)
    s = zeros(Float64, N)

    clade1 = string.(dfit.clade1)
    clade2 = string.(dfit.clade2)

    s .= sqrt.(dfit.std_fit1 .^ 2 .+ dfit.std_fit2 .^ 2)
    ddf .= mutcp_ddf(J, mut, cdiff, clade1, clade2)

    z = ddf ./ s

    return z

end

# Function to compute the model predicted fitness discrepancies for a specific mutation and a list of clade pairs
function mutcp_ddf(J::DataFrame, mut::S1, cdiff::DataFrame, clade1::Vector{S2}, clade2::Vector{S2}) where {S1<:AbstractString,S2<:AbstractString}

    @assert length(clade1) == length(clade2)
    ddf = zeros(Float64, length(clade1))

    for c in eachindex(clade1)
        ddf[c] = mutcp_ddf(J, mut, cdiff, clade1[c], clade2[c])
    end

    return ddf

end

# Function to compute the model predicted fitness discrepancies for a list of mutations and a pair of clades
function mutcp_ddf(J::DataFrame, mut::Vector{S1}, cdiff::DataFrame, clade1::S2, clade2::S2) where {S1<:AbstractString,S2<:AbstractString}

    ddf = zeros(Float64, length(mut))

    Threads.@threads for n in eachindex(mut)
        ddf[n] = mutcp_ddf(J, mut[n], cdiff, clade1, clade2)
    end

    return ddf

end

# Functions to compute the model predicted fitness discrepancies for a specific mutation and a pair of clades
function mutcp_ddf(J::DataFrame, mut::S1, cdiff::DataFrame, c1::S2, c2::S2) where {S1<:AbstractString,S2<:AbstractString}

    if sum((cdiff.clade1 .== c1) .& (cdiff.clade2 .== c2)) != 0
        cdiff_c1c2 = cdiff[(cdiff.clade1.==c1).&(cdiff.clade2.==c2), :]
    else
        cdiff_c1c2 = cdiff[(cdiff.clade1.==c2).&(cdiff.clade2.==c1), :]
    end

    site = parse(Int, mut[2:end-1])
    s_wt = string(mut[1])
    s_mut = string(mut[end])
    Jmut = J[(J.i.==site).&(J.σᵢ_wt.==s_wt).&(J.σᵢ.==s_mut), :]
    Jmut_c1c2 = Jmut[findall(x -> x in cdiff_c1c2.site, Jmut.j), :]

    ddf = mutcp_ddf(Jmut_c1c2, cdiff_c1c2, c1, c2)

    return ddf

end

function mutcp_ddf(J::DataFrame, cdiff::DataFrame, c1::S, c2::S) where {S<:AbstractString}

    @assert (cdiff.clade1[1] == c1 && cdiff.clade2[1] == c2) || (cdiff.clade1[1] == c2 && cdiff.clade2[1] == c1)

    ddf = 0.0
    fact = cdiff.clade1[1] == c1 ? 1 : -1

    for r in eachrow(J)
        ddf += fact * r.J * (r.σⱼ == cdiff[(cdiff.site.==r.j), :aa_c1][1])
        ddf += -fact * r.J * (r.σⱼ == cdiff[(cdiff.site.==r.j), :aa_c2][1])
    end

    return ddf

end

# Function for computing the energy of an inferred model
function energy(J::DataFrame, dfit::DataFrame, cdiff::DataFrame)

    energy = zeros(Float64, size(dfit, 1))
    w = SC2Epi.weight(dfit.std_fit1, dfit.std_fit2)

    @assert length(energy) == length(w)

    for n in 1:size(dfit, 1)

        mut = dfit.aa_mut[n]
        site = parse(Int, mut[2:end-1])
        s_wt = string(mut[1])
        s_mut = string(mut[end])
        Jmut = J[(J.i.==site).&(J.σᵢ_wt.==s_wt).&(J.σᵢ.==s_mut), :]
        c1 = string(dfit.clade1[n])
        c2 = string(dfit.clade2[n])
        j_ddf = mutcp_ddf(Jmut, mut, cdiff, c1, c2)
        ddf = dfit.fit1[n] - dfit.fit2[n]
        energy[n] = (ddf - j_ddf)^2

    end

    return dot(w, energy)

end