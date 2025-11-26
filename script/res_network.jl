""" Script to generate the networks associated to Spike domains: RBD, NTD and S2 unit """

# Import packages
using DataFrames, CSV

# Functions

# Compute the number of connections for all available residues (nodes of the network)
function network_connectivity(Jmat::Vector{Vector{Union{Vector{Int},Float64}}})

    i = unique(map(x -> x[1][1], Jmat))
    j = unique(map(x -> x[1][2], Jmat))

    res = sort(unique(vcat(i, j)))

    g = DataFrame(res=Int[], connect=Int[])

    for r in res
        c = 0
        for n in eachindex(Jmat)
            if r in Jmat[n][1][1] || r in Jmat[n][1][2]
                c += 1
            end
        end
        push!(g, [r, c])
    end

    return g

end

# Associate the coupling magnitude to each residue pair (network links)
function network_edges(Jmat::Vector{Vector{Union{Vector{Int},Float64}}})

    edges = DataFrame(i=Int[], j=Int[], J=Float64[])

    for n in eachindex(Jmat)
        push!(edges, [Jmat[n][1][1], Jmat[n][1][2], Jmat[n][2]])
    end

    return edges

end

# Load inferred couplings
Jtab = CSV.read("results/jcoup_l1_1em3_3d.csv", DataFrame)

# Compute the Frobenius norm of the couplings
Jmat = SC2Epistasis.Jmat(Jtab)

# Interaction network for the RBD
rbd_edge = [319, 528]  # RBD boundaries
J_rbd = Jmat[findall(x -> (rbd_edge[1] <= x[1][1] <= rbd_edge[2]) && (rbd_edge[1] <= x[1][2] <= rbd_edge[2]) && x[2] >= 1.25, Jmat)]
g_rbd = network_connectivity(J_rbd)
edges_rbd = network_edges(J_rbd)
CSV.write("results/network/rbd_connect.csv", g_rbd)
CSV.write("results/network/rbd_edges.csv", edges_rbd)

# Compute nodes connectivity

# Interaction network for the NTD
ntd_edge = [14, 305]  # NTD boundaries
J_ntd = Jmat[findall(x -> (ntd_edge[1] <= x[1][1] <= ntd_edge[2]) && (ntd_edge[1] <= x[1][2] <= ntd_edge[2]) && x[2] >= 1.25, Jmat)]
g_ntd = network_connectivity(J_ntd)
edges_ntd = network_edges(J_ntd)
CSV.write("results/network/ntd_connect.csv", g_ntd)
CSV.write("results/network/ntd_edges.csv", edges_ntd)

# Interaction network for the S2 unit + S1/S2 junction
s2_edge = [529, 1273]  # S2 boundaries
J_s2 = Jmat[findall(x -> (s2_edge[1] <= x[1][1] <= s2_edge[2]) && (s2_edge[1] <= x[1][2] <= s2_edge[2]) && x[2] >= 1.25, Jmat)]
g_s2 = network_connectivity(J_s2)
edges_s2 = network_edges(J_s2)
CSV.write("results/network/s2_connect.csv", g_s2)
CSV.write("results/network/s2_edges.csv", edges_s2)