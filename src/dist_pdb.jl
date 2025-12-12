"""
File for handling PDB files and calculating distances between residues
"""

using PdbTool
using BioAlignments

# Function to map all PDB files to Spike reference sequence
function map_pdbs!(pdbs::Vector{PdbTool.Pdb}, ref_seq::Vector{String}; mappedTo::String="ref_seq", dom_shift::Int64=1)

    idx_ref_dic = [Dict{String,Vector{Int64}}() for k in eachindex(pdbs)]

    for (k, p) in enumerate(pdbs)
        idx_ref_dic[k] = SC2Epistasis.map_pdb!(p, ref_seq; mappedTo=mappedTo, dom_shift=dom_shift)
    end

    return idx_ref_dic

end

# Mapping each chain to its reference sequence
function map_pdb!(pdb::PdbTool.Pdb, ref::Vector{String}; mappedTo::String="ref_seq", dom_shift::Int64=1)

    idx_ref = [Int64[] for k in 1:length(pdb.chain)]
    idx_ref_dic = Dict{String,Vector{Int64}}()

    l = 0
    for k in keys(pdb.chain)
        l += 1
        pdb_seq = PdbTool.chainSeq(pdb.chain[k])
        idx_ref[l] = map_pdb!(pdb.chain[k], pdb_seq, ref[l]; mappedTo=mappedTo, dom_shift=dom_shift)
        push!(idx_ref_dic, k => idx_ref[l])
    end

    return idx_ref_dic

end

# Mapping selected chains to reference sequences
function map_pdb!(pdb::PdbTool.Pdb, keys::Vector{String}, ref::Vector{String}; mappedTo::String="ref_seq", dom_shift::Int64=1)

    idx_ref = [Int64[] for k in 1:length(pdb.chain)]
    idx_ref_dic = Dict{String,Vector{Int64}}()

    @assert length(keys) == length(ref)

    l = 0
    for k in keys
        l += 1
        pdb_seq = PdbTool.chainSeq(pdb.chain[k])
        idx_ref[l] = map_pdb!(pdb.chain[k], pdb_seq, ref[l]; mappedTo=mappedTo, dom_shift=dom_shift)
        push!(idx_ref_dic, k => idx_ref[l])
    end

    return idx_ref_dic

end

# Function to map a PDB chain onto a reference sequence
function map_pdb!(chain::PdbTool.Chain, pdb_seq::String, ref::String; mappedTo::String="ref_seq", dom_shift::Int64=1)

    al_model = BioAlignments.OverlapAlignment()
    scoremodel = AffineGapScoreModel(match=5, mismatch=-4, gap_open=-4, gap_extend=-1)

    res = BioAlignments.pairalign(al_model, pdb_seq, ref, scoremodel)
    aln = BioAlignments.alignment(res)

    if length(pdb_seq) >= length(aln)
        idx_ovrlp = fill(1, count_matches(aln) + count_mismatches(aln))
        counter = 0
        for k in 1:length(aln)
            if aln2seq(aln, k)[2] == OP_SEQ_MATCH || aln2seq(aln, k)[2] == OP_SEQ_MISMATCH
                counter += 1
                idx_ovrlp[counter] = k
            end
        end
        idx_ref = collect(1:1:length(ref))
        chain_idx = sort(parse.(Int, collect(keys(chain.residue))))[idx_ovrlp]
        if chain_idx[1] == idx_ovrlp[1]
            idx_match = (chain_idx .== idx_ovrlp)
        elseif chain_idx[1] - dom_shift + 1 == idx_ovrlp[1]
            idx_match = (chain_idx .== idx_ovrlp .+ dom_shift .- 1)
        else
            println("Something's weird")
        end
    else
        idx_ref = fill(1, count_matches(aln) + count_mismatches(aln))
        counter = 0
        for k in 1:length(aln)
            if aln2seq(aln, k)[2] == OP_SEQ_MATCH || aln2seq(aln, k)[2] == OP_SEQ_MISMATCH
                counter += 1
                idx_ref[counter] = k
            end
        end
        idx_ovrlp = collect(1:1:count_matches(aln)+count_mismatches(aln))
        chain_idx = sort(parse.(Int, collect(keys(chain.residue))))[idx_ovrlp]
        if chain_idx[1] == idx_ref[1]
            idx_match = (chain_idx .== idx_ref)
        elseif chain_idx[1] - dom_shift + 1 == idx_ref[1]
            idx_match = (chain_idx .== idx_ref .+ dom_shift .- 1)
        else
            println("Something's weird")
        end
    end

    chain_idx = string.(chain_idx[idx_match])
    idx_ref = idx_ref[idx_match]

    off_set = parse(Int, chain_idx[1])

    for k in chain_idx
        idx_pdb = parse(Int, k)
        idx = idx_pdb - off_set + 1
        chain.residue[k].alignmentPos = idx
        chain.align[idx] = chain.residue[k]
        chain.align[idx].identifier = string(idx)
    end

    chain.mappedTo = mappedTo

    return idx_ref

end

# Given a PDB compute the distances between all possible residue pairs intra/inter chain
function dist_pdb(pdb::PdbTool.Pdb, idx_ref::Dict{String,Vector{Int64}})

    dist = DataFrame(i=Int[], j=Int[], d=Float64[], ch1=String[], ch2=String[])

    for ch1 in keys(idx_ref)
        for ch2 in keys(idx_ref)
            if ch1 == ch2
                L = length(idx_ref[ch1])
                for i in 1:L-1
                    for j in i+1:L
                        res1 = idx_ref[ch1][i]
                        res2 = idx_ref[ch1][j]
                        d = PdbTool.residueDist(pdb.chain[ch1].residue[string(res1)], pdb.chain[ch2].residue[string(res2)])
                        push!(dist, [res1, res2, d, ch1, ch2])
                    end
                end
            else
                for i in idx_ref[ch1]
                    for j in idx_ref[ch2]
                        d = PdbTool.residueDist(pdb.chain[ch1].residue[string(i)], pdb.chain[ch2].residue[string(j)])
                        push!(dist, [i, j, d, ch1, ch2])
                    end
                end
            end
        end
    end

    return dist

end