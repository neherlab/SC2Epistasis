"""
    Functions for finding mistmatches between clade founders and
    generating a dataframe with mutational fitness for mismatching
    clades
"""

# Function for finding mutations for a specific protein and clade
function clade_prot_mut(prot::S1, c::S2, clade_muts::Dict) where {S1<:AbstractString,S2<:AbstractString}

    aamut_c = clade_muts[c]["aa"]

    aamut_c_p = []
    for m in aamut_c
        p, mut = split(m, ":")
        if p == prot
            push!(aamut_c_p, (string(mut[1]), parse(Int, mut[2:end-1]), string(mut[end])))
        end
    end

    return aamut_c_p

end

# Function for finding mutations and deletions for a specific protein and clade
function clade_prot_mut(prot::S1, c::S2, clade_muts::Dict, clade_del::DataFrame) where {S1<:AbstractString,S2<:AbstractString}

    aamut_c = clade_muts[c]["aa"]
    del_pc = filter(row -> row.protein == prot && row.clade == c, clade_del)

    aamut_c_p = Tuple{String,Int,String}[]
    for m in aamut_c
        p, mut = split(m, ":")
        if p == prot
            push!(aamut_c_p, (string(mut[1]), parse(Int, mut[2:end-1]), string(mut[end])))
        end
    end

    for r in eachrow(del_pc)
        push!(aamut_c_p, (r["aa1"], r["position"], r["aa2"]))
    end

    return aamut_c_p

end

# Finding mismatching residues between two clades for a specific protein
function clade_pair_diff(prot::S1, clade1::S2, clade2::S2, clade_muts::Dict, clade_del::DataFrame) where {S1<:AbstractString,S2<:AbstractString}

    aamut_c1 = clade_prot_mut(prot, clade1, clade_muts, clade_del)
    aamut_c2 = clade_prot_mut(prot, clade2, clade_muts, clade_del)

    site_mut = clade_pair_diff(aamut_c1, aamut_c2)

    return site_mut

end

# Finding index of mismatch between two clades
function clade_pair_diff(aamut_c1, aamut_c2)

    nmut_1 = length(aamut_c1)
    nmut_2 = length(aamut_c2)

    daamut_c1 = Dict([aamut_c1[k][2] => aamut_c1[k][3] for k in 1:nmut_1])
    daamut_c2 = Dict([aamut_c2[k][2] => aamut_c2[k][3] for k in 1:nmut_2])

    site_mut = sort(unique(vcat([aamut_c1[k][2] for k in 1:nmut_1], [aamut_c2[k][2] for k in 1:nmut_2])))

    idx_diff = [true for k in 1:length(site_mut)]

    for s in eachindex(site_mut)
        if haskey(daamut_c1, site_mut[s]) && haskey(daamut_c2, site_mut[s]) && daamut_c1[site_mut[s]] == daamut_c2[site_mut[s]]
            idx_diff[s] = false
        end
    end

    return site_mut[idx_diff]

end

# Given a list of proteins and clades (with mutations/deletions) provides a dataframe for mistmacthes
function clade_diff_list(prot_list::Vector{S1}, clades::Vector{S2}, clade_muts::Dict; clade_del::DataFrame=DataFrame()) where {S1<:AbstractString,S2<:AbstractString}

    Nc = length(clades)

    # Defining dataframe with mismatches between clades 
    clades_diff = DataFrame(prot=String[], site=Int64[], clade1=String[], aa_c1=String[], clade2=String[], aa_c2=String[])

    for p in prot_list

        for c1 in 1:Nc-1
            for c2 in c1+1:Nc

                if isempty(clade_del)
                    mut_c1 = clade_prot_mut(p, clades[c1], clade_muts)
                    mut_c2 = clade_prot_mut(p, clades[c2], clade_muts)
                else
                    mut_c1 = clade_prot_mut(p, clades[c1], clade_muts, clade_del)
                    mut_c2 = clade_prot_mut(p, clades[c2], clade_muts, clade_del)
                end


                idx_mis_c1c2 = clade_pair_diff(mut_c1, mut_c2)

                if !isempty(idx_mis_c1c2)
                    dmut_c1 = Dict([mut_c1[k][2] => [mut_c1[k][1], mut_c1[k][3]] for k in eachindex(mut_c1)])
                    dmut_c2 = Dict([mut_c2[k][2] => [mut_c2[k][1], mut_c2[k][3]] for k in eachindex(mut_c2)])
                    for i in idx_mis_c1c2
                        if haskey(dmut_c1, i) && haskey(dmut_c2, i)
                            aa_c1 = dmut_c1[i][2]
                            aa_c2 = dmut_c2[i][2]
                            @assert aa_c1 != aa_c2
                        elseif !haskey(dmut_c1, i)
                            aa_c1 = dmut_c2[i][1]
                            aa_c2 = dmut_c2[i][2]
                        else
                            aa_c1 = dmut_c1[i][2]
                            aa_c2 = dmut_c1[i][1]
                        end
                        push!(clades_diff, [p, i, clades[c1], aa_c1, clades[c2], aa_c2])
                    end
                end

            end
        end

    end

    return clades_diff

end

function delta_fit(prots::Vector{S1}, clades::Vector{S2}, aamut_fit::DataFrame, clade_muts::Dict{String,Any}, clade_del::DataFrame;
    cnt_thr::Float64=5.0,
    exclude_stop::Bool=true,
    exclude_syn::Bool=true) where {S1<:AbstractString,S2<:AbstractString}

    # Retaining only mutations with a number of expected counts larger than threshold
    aamut_fit = aamut_fit[aamut_fit.predicted_count.>=cnt_thr, :]

    # Selecting only non-synonymous mutations and no stop-codon
    if exclude_syn
        aamut_fit = aamut_fit[findall(x -> (x[1] != x[end]), aamut_fit.aa_mutation), :]
    end
    if exclude_stop
        aamut_fit = aamut_fit[findall(x -> (x[end] != '*'), aamut_fit.aa_mutation), :]
    end

    # Sort clades
    sort_clades!(clades)

    # Initialize empty dataframe
    delta_fit = DataFrame(prot=String[], clade1=String[], clade2=String[], aa_mut=String[], fit1=Float64[], fit2=Float64[], std_fit1=Float64[], std_fit2=Float64[],
        actual_count1=Int64[], actual_count2=Int64[], exp_count1=Float64[], exp_count2=Float64[], p_value=Float64[])

    # Iterate over proteins
    for p in prots
        aamut_fit_p = aamut_fit[aamut_fit.gene.==p, :]
        dfit_mut = SC2Epistasis.delta_fit(p, clades, aamut_fit_p, clade_muts, clade_del)
        delta_fit = vcat(delta_fit, dfit_mut)
    end

    return delta_fit

end

function delta_fit(protein::S1, clades::Vector{S2}, aamut_fit::DataFrame, clade_muts::Dict{String,Any}, clade_del::DataFrame) where {S1<:AbstractString,S2<:AbstractString}

    N_c = length(clades)

    # Initialize empty arrays for dataframe columns
    prot = String[]
    clade1 = String[]
    clade2 = String[]
    aa_mut = String[]
    fit1 = Float64[]
    fit2 = Float64[]
    std_fit1 = Float64[]
    std_fit2 = Float64[]
    actual_count1 = Int64[]
    actual_count2 = Int64[]
    exp_count1 = Float64[]
    exp_count2 = Float64[]
    p_value = Float64[]

    @assert protein == (length(unique(aamut_fit.gene)) == 1 ? unique(aamut_fit.gene)[1] : error("Multiple proteins for mutation dataframe"))

    # Populate columns
    for c1 in 1:N_c-1
        for c2 in c1+1:N_c

            # Getting mistmatches between clade founders
            mut_c1 = clade_prot_mut(protein, clades[c1], clade_muts, clade_del)
            mut_c2 = clade_prot_mut(protein, clades[c2], clade_muts, clade_del)

            # Mistmatching residues between clades
            idx_mis_c1c2 = clade_pair_diff(mut_c1, mut_c2)

            if isempty(idx_mis_c1c2)
                continue
            end

            # Isolate fitness data for clades
            fit_c1 = aamut_fit[aamut_fit.clade.==clades[c1], [:clade_founder_aa, :mutant_aa, :aa_site, :aa_mutation, :predicted_count, :actual_count, :delta_fitness, :uncertainty]]
            fit_c2 = aamut_fit[aamut_fit.clade.==clades[c2], [:clade_founder_aa, :mutant_aa, :aa_site, :aa_mutation, :predicted_count, :actual_count, :delta_fitness, :uncertainty]]

            # Exclude mutations from mistmatching sites
            fit_c1 = fit_c1[findall(x -> !(x in idx_mis_c1c2), fit_c1.aa_site), :]
            fit_c2 = fit_c2[findall(x -> !(x in idx_mis_c1c2), fit_c2.aa_site), :]

            # Common mutations between clades
            idx_c1 = findall(x -> x in fit_c2.aa_mutation, fit_c1.aa_mutation)
            idx_c2 = Int[]
            for i in idx_c1
                push!(idx_c2, findall(x -> x == fit_c1.aa_mutation[i], fit_c2.aa_mutation)[1])
            end

            @assert length(idx_c1) == length(idx_c2)

            # Populate dataframe entries
            for i in eachindex(idx_c1)

                @assert fit_c1[idx_c1[i], :aa_mutation] == fit_c2[idx_c2[i], :aa_mutation]
                aamut = fit_c1[idx_c1[i], :aa_mutation]
                f1 = fit_c1[idx_c1[i], :delta_fitness]
                f2 = fit_c2[idx_c2[i], :delta_fitness]
                std1 = fit_c1[idx_c1[i], :uncertainty]
                std2 = fit_c2[idx_c2[i], :uncertainty]
                cnt_c1 = fit_c1[idx_c1[i], :actual_count]
                cnt_c2 = fit_c2[idx_c2[i], :actual_count]
                expcnt_c1 = fit_c1[idx_c1[i], :predicted_count]
                expcnt_c2 = fit_c2[idx_c2[i], :predicted_count]
                pval = pvalue(FisherExactTest(cnt_c1, cnt_c2, Int(round(expcnt_c1)), Int(round(expcnt_c2))))
                push!(prot, protein)
                push!(clade1, clades[c1])
                push!(clade2, clades[c2])
                push!(aa_mut, aamut)
                push!(fit1, f1)
                push!(fit2, f2)
                push!(std_fit1, std1)
                push!(std_fit2, std2)
                push!(actual_count1, cnt_c1)
                push!(actual_count2, cnt_c2)
                push!(exp_count1, expcnt_c1)
                push!(exp_count2, expcnt_c2)
                push!(p_value, pval)

            end

        end
    end

    # Convert into dataframe
    delta_fit = DataFrame(prot=prot, clade1=clade1, clade2=clade2, aa_mut=aa_mut, fit1=fit1, fit2=fit2, std_fit1=std_fit1, std_fit2=std_fit2,
        actual_count1=actual_count1, actual_count2=actual_count2, exp_count1=exp_count1, exp_count2=exp_count2, p_value=p_value)

    return delta_fit

end

function sort_clades!(clades::Vector{S}) where {S<:AbstractString}
    # Custom comparator function
    function clade_comparator(a::S, b::S) where {S<:AbstractString}
        # Extract numeric and alphabetic parts
        num_a, alpha_a = match(r"(\d+)(\D+)", a).captures
        num_b, alpha_b = match(r"(\d+)(\D+)", b).captures

        # Compare numeric parts first, then alphabetic parts
        if parse(Int, num_a) == parse(Int, num_b)
            return alpha_a < alpha_b
        else
            return parse(Int, num_a) < parse(Int, num_b)
        end
    end

    # Sort using the custom comparator
    sort!(clades, lt=clade_comparator)
end