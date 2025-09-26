""" Utility functions for performing inference and testing """

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