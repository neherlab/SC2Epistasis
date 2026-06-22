""" 
    Functions to quantify the raw epistatic signal
    contained in the fitness discrepancies data
"""

# Compute the raw epistatic signal for all possible clade pairs
# Input: clades list, dataframe with fitness discrepancies and dataframe with clade mismatches
# Output: dictionaries with z-scores and associated sites
function raw_signal(clades::Vector{S}, dfit::DataFrame, cdiff::DataFrame) where {S<:AbstractString}

    # Initialize dictionaries
    z_dict = Dict{Tuple{S,S},Vector{Float64}}()
    s_dict = Dict{Tuple{S,S},Vector{Int64}}()

    # Loop over all possible clade pairs
    for (clade1, clade2) in collect(combinations(clades, 2))

        # Compute z-scores and sites list
        z, sites = z_fit(dfit, (clade1, clade2))
        # Skip if empty
        if isempty(z) || isempty(sites)
            continue
        end
        # Populate dictionaries
        push!(s_dict, (clade1, clade2) => sites)
        push!(z_dict, (clade1, clade2) => z)

    end

    return z_dict, s_dict

end

# Function for computing the overall epistatic signal between a pair of clades
function z_fit(dfit::DataFrame, cpair::Tuple{S,S}) where {S<:AbstractString}

    # Isolate fitness entries for a specific clusters pair
    dfit_c1c2 = dfit[(dfit.clade1 .== cpair[1]) .& (dfit.clade2 .== cpair[2]), :]
    if isempty(dfit_c1c2)
        dfit_c1c2 = dfit[(dfit.clade1 .== cpair[2]) .& (dfit.clade2 .== cpair[1]), :]
    end

    # List of sites with available fitness measurements
    sites = unique(parse.(Int, map(x -> x[2:(end-1)], dfit_c1c2.aa_mut)))
    sort!(sites)

    # Number of available mutations for each site
    num_aamut = zeros(Int, length(sites))
    for n in eachindex(num_aamut)
        num_aamut[n] = size(dfit_c1c2[findall(x -> parse(Int, x[2:(end-1)]) == sites[n], dfit_c1c2.aa_mut), :], 1)
    end

    z = zeros(length(sites)) #z-score

    for i in eachindex(sites)
        dfit_c1c2_site = dfit_c1c2[findall(x -> parse(Int, x[2:(end-1)]) == sites[i], dfit_c1c2.aa_mut), :]
        @assert size(dfit_c1c2_site, 1) == num_aamut[i]
        std1 = dfit_c1c2_site.std_fit1
        std2 = dfit_c1c2_site.std_fit2
        st_dev = sqrt.(std1 .^ 2 .+ std2 .^ 2)
        f1 = dfit_c1c2_site.fit1
        f2 = dfit_c1c2_site.fit2
        z[i] = sum(abs.(f1 .- f2) ./ (st_dev)) / num_aamut[i]
    end

    return z, sites

end

# Compute the normalized fitness discrepancies over a dataframe
function z_dfit(dfit::DataFrame)

    z = zeros(Float64, size(dfit, 1))

    z .= (dfit.fit1 .- dfit.fit2) ./ sqrt.(dfit.std_fit1 .^ 2 .+ dfit.std_fit2 .^ 2)

    return z

end

# Computes the fraction of sites with epistatic signal above a threshold
# which have a mismatch in a sphere of radius r for a grid of radii and z-thresholds
# This version uses a pre-computed distance matrix to speed up the computation.
function frac_in_sphere(cpair::Tuple{S1,S1}, fit_epi::Dict{Tuple{S2,S2},Vector{Float64}}, sites::Dict{Tuple{S2,S2},Vector{Int64}},
    clade_diff::DataFrame, dist_mat::Matrix{Float64}, radii::Vector{Float64}, z_thr::Vector{Float64}) where {S1<:AbstractString,S2<:AbstractString}

    if !haskey(fit_epi, cpair)
        throw(KeyError)
    end

    frac = zeros(Float64, length(radii), length(z_thr))

    for r in eachindex(radii)
        for z in eachindex(z_thr)
            frac[r, z], _ = nmis_sphere(cpair, fit_epi, sites, clade_diff, dist_mat; z_thr=z_thr[z], d_thr=radii[r])
        end
    end

    return frac

end

# Look for mutations with epistatic signal above threshold and assess
# the presence of background mismatch sites into a sphere of radius r
# This version uses a pre-computed distance matrix to speed up the computation.
function nmis_sphere(cpair::Tuple{S1,S1}, fit_epi::Dict{Tuple{S2,S2},Vector{Float64}}, sites::Dict{Tuple{S2,S2},Vector{Int64}},
    clade_diff::DataFrame, dist_mat::Matrix{Float64};
    z_thr::Float64=3.0,
    d_thr::Float64=10.0) where {S1<:AbstractString,S2<:AbstractString}

    if !haskey(fit_epi, cpair)
        throw(KeyError)
    end

    idx_epi = findall(x -> x >= z_thr, fit_epi[cpair])
    s_mut = sites[cpair][idx_epi]
    s_mis = clade_diff[(clade_diff.clade1 .== cpair[1]) .& (clade_diff.clade2 .== cpair[2]), :site]

    nmut = length(s_mut)

    cnt_neighbour = zeros(Int64, length(idx_epi))
    frac = 0.0

    for i in eachindex(s_mut)
        for j in eachindex(s_mis)
            res1 = s_mut[i]
            res2 = s_mis[j]
            d = dist_mat[res1, res2]
            if d <= d_thr
                cnt_neighbour[i] += 1
            end
        end
    end

    frac = sum(cnt_neighbour .> 0) / nmut

    return frac, cnt_neighbour

end

# Computes the distance between two residues in a set of PDBs.
# If the residues are not found in the PDBs, it uses the AF PDB as a fallback.
function dist_res(res1::Int64, res2::Int64, pdbs::Vector{PdbTool.Pdb}, af_pdb::PdbTool.Pdb)

    dmin = Inf

    for pdb in pdbs
        for ch1 in keys(pdb.chain)
            for ch2 in keys(pdb.chain)
                if haskey(pdb.chain[ch1].residue, string(res1)) && haskey(pdb.chain[ch2].residue, string(res2))
                    d = PdbTool.residueDist(pdb.chain[ch1].residue[string(res1)], pdb.chain[ch2].residue[string(res2)])
                else
                    d = PdbTool.residueDist(af_pdb.chain[ch1].residue[string(res1)], af_pdb.chain[ch2].residue[string(res2)])
                end
                if d < dmin
                    dmin = d
                end
            end
        end
    end

    return dmin

end

# Function to compute 3D distance according to a single AlphaFold structure
function dist_res(i::Int, j::Int, pdb::PdbTool.Pdb; chain::String="A")

    res1 = pdb.chain[chain].residue[string(i)]
    res2 = pdb.chain[chain].residue[string(j)]

    dist = PdbTool.residueDist(res1, res2)

    return dist

end

# Compute the distance matrix for a set of PDBs and an AlphaFold structure
function compute_dist_matrix(pdbs::Vector{PdbTool.Pdb}, af_pdb::PdbTool.Pdb)

    res = collect(1:1:1273) #set of all residues in the Spike protein
    dist_mat = SC2Epistasis.dist_res.(reshape(res, :, 1), reshape(res, 1, :), Ref(pdbs), Ref(af_pdb))

    return dist_mat

end

# Random benchmark: for a given pair of clades generates random mismtaches on the protein chain
# drawing from the pool defined by `var_sites`.
# It then computes the fraction of sites with z-score above threshold for which a mismatch is found
# in spheres of size specified by `radii`.
function rand_frac_in_sphere(cpair::Tuple{S1,S1}, fit_epi::Dict{Tuple{S2,S2},Vector{Float64}}, sites::Dict{Tuple{S2,S2},Vector{Int64}},
    clade_diff::DataFrame, dist_mat::Matrix{Float64}, radii::Vector{Float64}, z_thr::Vector{Float64};
    var_sites::Vector{Int64}=collect(1:size(dist_mat, 1)),
    nsamp::Int64=100) where {S1<:AbstractString,S2<:AbstractString}

    if !haskey(fit_epi, cpair)
        throw(KeyError)
    end

    frac = zeros(Float64, length(radii), length(z_thr))
    f = zeros(Float64, length(radii), length(z_thr), nsamp)
    std_frac = zeros(Float64, length(radii), length(z_thr))

    Threads.@threads for r in eachindex(radii)
        Threads.@threads for z in eachindex(z_thr)
            f[r, z, :], _ = rand_nmis_sphere(cpair, fit_epi, sites, clade_diff, dist_mat; var_sites=var_sites, z_thr=z_thr[z], d_thr=radii[r], nsamp=nsamp)
            frac[r, z] = mean(f[r, z, :])
            std_frac[r, z] = std(f[r, z, :]) / sqrt(nsamp)
        end
    end

    return frac, std_frac

end

# Random background mismatches in a sphere of radius r
function rand_nmis_sphere(cpair::Tuple{S1,S1}, fit_epi::Dict{Tuple{S2,S2},Vector{Float64}}, sites::Dict{Tuple{S2,S2},Vector{Int64}},
    clade_diff::DataFrame, dist_mat::Matrix{Float64};
    var_sites::Vector{Int64}=collect(1:size(dist_mat, 1)),
    z_thr::Float64=1.5,
    d_thr::Float64=10.0,
    nsamp::Int64=100) where {S1<:AbstractString,S2<:AbstractString}

    if !haskey(fit_epi, cpair)
        throw(KeyError)
    end

    idx_epi = findall(x -> x >= z_thr, fit_epi[cpair])
    s_mut = sites[cpair][idx_epi]
    nmut = length(s_mut)
    nmis = length(clade_diff[(clade_diff.clade1 .== cpair[1]) .& (clade_diff.clade2 .== cpair[2]), :site])
    smis = clade_diff[(clade_diff.clade1 .== cpair[1]) .& (clade_diff.clade2 .== cpair[2]), :site]

    frac = zeros(Float64, nsamp) # fraction of sites within sphere for each replicate
    cnt_neighbor = zeros(Float64, nmut, nsamp) # count of neighbors within threshold for each replicate and site
    j_rand = rand(setdiff(var_sites, vcat(s_mut, smis)), nmis, nsamp) # random mismatches

    cnt_neighbor .= dropdims(sum(reshape(dist_mat[s_mut, vec(j_rand)], nmut, nmis, nsamp) .<= d_thr, dims=2), dims=2)
    frac .= vec(mean(cnt_neighbor .>= 1, dims=1))

    return frac, cnt_neighbor

end