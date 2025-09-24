""" Utility functions for performing inference and testing """

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