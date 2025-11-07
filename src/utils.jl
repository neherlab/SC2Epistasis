""" Utility functions for performing inference and testing """

# Creates a table of epistatic interactions:
# i --> mutated site
# σʷᵗᵢ --> w.t. amino acid
# σᵢ --> arrival amino acid
# j --> interacting site
# σⱼ --> interacting amino acid
# J --> coupling parameter
function epistatic_table(data::Vector{Data}, optx::Vector{OptX})

    mutations = [data[m].dfit_mut.aa_mut[1] for m in eachindex(data)]

    return epistatic_table(mutations, optx)

end

function epistatic_table(mutations::Vector{S}, optx::Vector{OptX}) where {S<:AbstractString}

    epi = [optx[m].epi for m in eachindex(optx)]

    epistatic_table(mutations, epi)

end

function epistatic_table(mutations::Vector{S}, couplings::Vector{Dict{Vector{Union{Int,String}},Float64}}) where {S<:AbstractString}

    n_row = 0
    for m in eachindex(mutations)
        n_row += length(couplings[m])
    end

    J_tab = (σᵢ_wt=Vector{Char}(undef, n_row),
        i=Vector{Int}(undef, n_row),
        σᵢ=Vector{Char}(undef, n_row),
        j=Vector{Int}(undef, n_row),
        σⱼ=Vector{Char}(undef, n_row),
        J=Vector{Float64}(undef, n_row))

    counter = 0
    for m in eachindex(mutations)
        for (key, value) in couplings[m]
            counter += 1
            J_tab.σᵢ_wt[counter] = mutations[m][1]
            J_tab.i[counter] = parse(Int, mutations[m][2:end-1])
            J_tab.σᵢ[counter] = mutations[m][end]
            J_tab.j[counter] = key[1]
            J_tab.σⱼ[counter] = only(key[2])
            J_tab.J[counter] = value
        end
    end

    return DataFrame(J_tab, copycols=false)

end

function frob_norm(Jtab::DataFrame)

    res_pair = extract_res(Jtab)

    return frob_norm(Jtab, res_pair)

end

function frob_norm(Jtab::DataFrame, res_pair::Vector{Tuple{Int,Int}})

    group_jtab = groupby(Jtab, [:i, :j])
    J_res = Dict{Tuple{Int,Int},Float64}()

    for r in res_pair
        J_res[(r[1], r[2])] = frob_norm(group_jtab, r)
    end

    return sort(J_res, byvalue=true, rev=true)

end

function frob_norm(group_jtab::GroupedDataFrame, res_pair::Tuple{Int,Int})

    g = get(group_jtab, (i=res_pair[1], j=res_pair[2]), nothing)
    g === nothing && return 0.0
    jfrob = sqrt(sum(g.J .^ 2))

    return jfrob

end

function extract_res(Jtab::DataFrame)

    res_pair = unique([(r.i, r.j) for r in eachrow(Jtab)])

    return res_pair

end

function Jmat(Jfrob)

    Jmat = [Union{Vector{Int},Float64}[] for k in 1:length(Jfrob)]
    counter = 0
    for (k, v) in pairs(Jfrob)
        counter += 1
        push!(Jmat[counter], [k...])
        push!(Jmat[counter], v)
    end

    return Jmat

end

# Read PDB files from the data folder
function read_pdbs()

    # Read PDB files
    pdbs = PdbTool.Pdb[]
    for file in readdir("data/PDB/")
        if endswith(file, ".pdb")
            push!(pdbs, PdbTool.parsePdb("data/PDB/" * file))
        end
    end

    return pdbs[1:end-1], pdbs[end] # return PDBs and AF2 structure

end

# Compute energy from structs OptX and QForm
function energy(optx::Vector{OptX}, qform::Vector{QForm})

    ener = 0.0
    for m in eachindex(optx)
        ener += SC2Epistasis.energy(optx[m], qform[m])
    end

    return ener

end

function energy(optx::OptX, qform::QForm)

    return dot(optx.J, qform.A, optx.J) + 2 * dot(qform.f, optx.J) + qform.c

end

# Function computing 3D distances between residues of each putative coupling for a set of PDB's
function threedist(optx::Vector{OptX}, data::Vector{Data}, pdbs::Vector{PdbTool.Pdb}, af_pdb::PdbTool.Pdb)

    M = length(data)

    mut_site = [parse(Int64, data[m].dfit_mut.aa_mut[1][2:end-1]) for m in 1:M]
    j_dim = sum([optx[m].num_j for m in 1:M])

    dist = zeros(Float64, j_dim)

    counter = 0
    for m in 1:M
        i = mut_site[m]
        for j in sort(collect(keys(optx[m].site_aa)))
            for aa in optx[m].site_aa[j]
                dmin = Inf
                counter += 1
                for pdb in pdbs
                    for ch1 in keys(pdb.chain)
                        for ch2 in keys(pdb.chain)
                            if string(i) in keys(pdb.chain[ch1].residue) && string(j) in keys(pdb.chain[ch2].residue)
                                d = PdbTool.residueDist(pdb.chain[ch1].residue[string(i)], pdb.chain[ch2].residue[string(j)])
                            else
                                d = PdbTool.residueDist(af_pdb.chain[ch1].residue[string(i)], af_pdb.chain[ch2].residue[string(j)])
                            end
                            if d < dmin
                                dmin = d
                            end
                        end
                    end
                end
                dist[counter] = dmin
            end
        end
    end

    return dist

end

# Given a set of mutations, interacting mismatching residues and inferred couplings, returns for each mutation a struct with:
# - List of clades 
# - List of fitness effects
# - List of standard deviations of fitness effects
# - List of background amino acids at interacting site
# - List of coupling parameters J for each amino acid
# - Wild-type amino acid at mutated site
# - Arrival amino acid at mutated site
# - Mutated site
# - Interacting site
function coup_dfit(Jtab::DataFrame, dfit::DataFrame, cdiff::DataFrame, muts::Vector{S}, res::Vector{Int}) where {S<:AbstractString}

    j_dfit = []
    for m in eachindex(muts)
        site = parse(Int, muts[m][2:end-1])
        # Extract couplings involving the mutation site
        Jmut = Jtab[(Jtab.σᵢ_wt.==string(muts[m][1])).&(Jtab.i.==site).&(Jtab.σᵢ.==string(muts[m][end])).&(Jtab.j.==res[m]), :]
        # Extract fitness discrepancies for the mutation
        dfit_mut = dfit[dfit.aa_mut.==muts[m], :]
        push!(j_dfit, coup_dfit(Jmut, dfit_mut, cdiff, res[m]))
    end

    return j_dfit

end

function coup_dfit(Jmut::DataFrame, dfit_mut::DataFrame, cdiff::DataFrame, res::Int)

    cp_mut = SC2Epistasis.clade_pair_mut(dfit_mut, cdiff)
    cp_mut = cp_mut[cp_mut.site.==res, :]

    clades = unique(vcat(dfit_mut.clade1, dfit_mut.clade2))
    fit = fill(0.0, length(clades))
    std_fit = fill(0.0, length(clades))
    sj = fill("", length(clades))
    J = fill(0.0, length(clades))

    for (i, c) in enumerate(clades)
        if c in dfit_mut.clade1
            idx_c = findall(dfit_mut.clade1 .== c)
            fit[i] = dfit_mut.fit1[idx_c][1]
            std_fit[i] = dfit_mut.std_fit1[idx_c][1]
        else
            idx_c = findall(dfit_mut.clade2 .== c)
            fit[i] = dfit_mut.fit2[idx_c][1]
            std_fit[i] = dfit_mut.std_fit2[idx_c][1]
        end
        if sum(cp_mut.clade1 .== c) != 0
            aa = cp_mut[cp_mut.clade1.==c, :aa_c1][1]
        else
            aa = cp_mut[cp_mut.clade2.==c, :aa_c2][1]
        end
        sj[i] = aa
        J[i] = Jmut[string.(Jmut.σⱼ).==aa, :J][1]
    end

    return (clades=clades, fit=fit, s_fit=std_fit, sj=sj, J=J, si_wt=Jmut.σᵢ_wt[1], si=Jmut.σᵢ[1], i=Jmut.i[1], j=res)

end