"""
    Structs used for optimization and associated functions
"""

# Optimization structs
export Data, OptX, QForm

""" Struct containing data-derived information about fitness effects and mismatches between clade pairs for a specific mutation """
mutable struct Data{S<:AbstractString,T<:Real}

    "`dfit_mut::DataFrame`: dataframe with fitness discrepancies for a specific mutation"
    dfit_mut::DataFrame
    "`clade_pairs::Vector{Vector{S}}`: list of unique clade pairs associated to the mutation"
    clade_pairs::Vector{Vector{S}}
    "`cpair_mut::DataFrame`: dataframe with mismatches between clade pairs"
    cpair_mut::DataFrame
    "`weights::Vector{T}`: weights associated to fitness discrepancies for each clade pair"
    weights::Vector{T}

    @doc "Inner constructor"
    function Data(dfit_mut::DataFrame, clade_pairs::Vector{Vector{S}}, cpair_mut::DataFrame,
        weights::Vector{T}) where {S<:AbstractString,T<:Real}
        new{S,T}(dfit_mut, clade_pairs, cpair_mut, weights)
    end

end

# Function to create a Data struct for a specific mutation
# Inputs:
# - mutation: mutation of interest
# - delta_fit: dataframe with fitness discrepancies for all mutations and clade pairs
# - clade_diff: dataframe with mismatches between clade pairs
function Data(mutation::S1, delta_fit::DataFrame, clade_diff::DataFrame;
    weighted::Bool=true,
    norm::Bool=true) where {S1<:AbstractString}

    dfit_mut = delta_fit[findall(x -> x == mutation, delta_fit.aa_mut), :]
    clade_pairs = unique_clade_pairs(dfit_mut)
    cpair_mut = clade_pair_mut(clade_pairs, clade_diff)
    weights = weight(dfit_mut; weighted=weighted, norm=norm)

    Data(dfit_mut, clade_pairs, cpair_mut, weights)

end

""" Struct containing coupling variables and associated features for a specific mutation """
mutable struct OptX{Q<:Int,T<:Real,S<:AbstractString}

    "`epi::Dict{Vector{Union{Q,S}},T}`: dictionary with site and amino acid as keys and couplings as values"
    epi::Dict{Vector{Union{Q,S}},T}
    "`site_aa::Dict{Q,Vector{S}}`: dictionary with sites as keys and unique amino acids across clade pairs as values"
    site_aa::Dict{Q,Vector{S}}
    "`num_j::Q`: dimension of the couplings vector for the mutation"
    num_j::Q
    "`J::Vector{T}`: vector with couplings"
    J::Vector{T}

    @doc "Inner constructor"
    function OptX(epi::Dict{Vector{Union{Q,S}},T}, site_aa::Dict{Q,Vector{S}}, num_j::Q, J::Vector{T}) where {Q<:Int,T<:Real,S<:AbstractString}
        new{Q,T,S}(epi, site_aa, num_j, J)
    end

end

# Function to create an OptX struct for a specific mutation
# Inputs:
# - data: Data struct for the mutation
function OptX(data::Data{S,T}) where {S<:AbstractString,T<:Real}

    epi = Dict{Vector{Union{Int,S}},T}()
    mut_site = parse(Int64, data.dfit_mut.aa_mut[1][2:end-1])
    site_aa, num_j = coupling_dim(data.cpair_mut, mut_site)
    for s in sort(collect(keys(site_aa)))
        for σ in site_aa[s]
            push!(epi, [s, σ] => 0.0)
        end
    end
    J = fill(T(0), num_j)

    OptX(epi, site_aa, num_j, J)

end

# Function to create an OptX struct for a specific mutation
# initializing couplings with a provided dictionary
# Inputs:
# - data: Data struct for the mutation
# - Jstart: dictionary with initial couplings
function OptX(data::Data, Jstart::Dict{Vector{Union{Q,S}},T}) where {S<:AbstractString,Q<:Int,T<:Real}


    epi = Dict{Vector{Union{Q,S}},T}()
    mut_site = parse(Int64, data.dfit_mut.aa_mut[1][2:end-1])
    site_aa, num_j = coupling_dim(data.cpair_mut, mut_site)
    counter = 0
    for s in sort(collect(keys(site_aa)))
        for σ in site_aa[s]
            counter += 1
            push!(epi, [s, σ] => Jstart[[s, σ]])
            J[counter] = Jstart[[s, σ]]
        end
    end

    OptX(epi, site_aa, num_j, J)

end

""" Structure containing elements of the quadratic form for a specific mutation """
struct QForm{T<:Real}

    "`A::Matrix{T}`: matrix of the quadratic form"
    A::Matrix{T}
    "`f::Vector{T}`: vector of the quadratic form"
    f::Vector{T}
    "`c::T`: constant term of the quadratic form"
    c::T

    @doc "Inner constructor"
    function QForm(A::Matrix{T}, f::Vector{T}, c::T) where {T<:Real}
        new{T}(A, f, c)
    end

end

# Function to create a QForm struct for a specific mutation
# Inputs:
# - data: Data struct for the mutation
# - opt_x: OptX struct for the mutation
function QForm(data::Data{S,T}, opt_x::OptX{Q,T,S}) where {Q<:Int,T<:Real,S<:String}

    w_vec = qf_coeff(data, opt_x)

    A = w_vec * w_vec'
    f = w_vec * ((data.dfit_mut[:, :fit2] .- data.dfit_mut[:, :fit1]) .* sqrt.(data.weights))
    c = offset(data.dfit_mut, data.weights)

    QForm(A, f, c)

end

# Function to create a QForm struct for a specific mutation
# Inputs:
# - w_vec: matrix of coefficients
# - dfit_mut: dataframe with fitness discrepancies for the mutation
# - weights: vector of weights associated to fitness discrepancies
function QForm(w_vec::Matrix{T}, dfit_mut::DataFrame, weights::Vector{T}) where {T<:Real}

    A = w_vec * w_vec'
    f = w_vec * ((dfit_mut[:, :fit2] .- dfit_mut[:, :fit1]) .* sqrt.(weights))
    c = offset(dfit_mut, weights)

    QForm(A, f, c)

end

##############################--------------------------------

# OWL-QN structs

# struct with key algorithm parameters
mutable struct OWLQN
    s                # param_t+1 - param_t [max size of m]
    y                # grad_t+1 - grad_t [max size of m]
    rho              # 1/s]i]'y[i]
    m::Int           # L-BFGS history length 
    t::Int           # iteration
    lambda::Vector{Float64}  # L1 penalty
end

# constructor
function OWLQN(x::Vector{Float64}, λ_vec::Vector{Float64})
    s = []
    y = []
    rho = []
    m = 6
    t = 0
    OWLQN(s, y, rho, m, t, λ_vec)
end

##############################--------------------------------

# Initialization of all structs

function init_all(muts::Vector{S}, delta_fit::DataFrame, clade_diff::DataFrame;
    weighted::Bool=true,
    norm::Bool=true) where {S<:AbstractString}

    # Initialize inference structs for all mutations
    data = Vector{Data}(undef, length(muts))
    optx = Vector{OptX}(undef, length(muts))
    qform = Vector{QForm}(undef, length(muts))

    # Populating data and optx
    Threads.@threads for m in eachindex(muts)
        data[m] = Data(muts[m], delta_fit, clade_diff; weighted=weighted, norm=false) # We need overall normalization later
        optx[m] = OptX(data[m])
    end

    w_norm = sum(vcat([data[m].weights for m in eachindex(muts)]...))

    # Populating qform and normalizing weights
    Threads.@threads for m in eachindex(muts)
        if norm
            data[m].weights ./= w_norm
        end
        qform[m] = QForm(data[m], optx[m])
    end

    return data, optx, qform

end

##############################--------------------------------

# Functions for struct initialization

# Function to define which and how many unknown couplings there are.
# Also defines a dictionary with sites as keys and unique amino acids
# across clade pairs as values.
function coupling_dim(cpair_mut::DataFrame, mut_site::Int64)

    site_aa = Dict{Int,Vector{String}}()

    sites = copy(cpair_mut.site)

    # Excluding self interactions
    idx_nov = findall(x -> x != mut_site, sites)

    for n in idx_nov

        if !(sites[n] in keys(site_aa))
            push!(site_aa, sites[n] => [cpair_mut.aa_c1[n], cpair_mut.aa_c2[n]])
        else
            if !(cpair_mut.aa_c1[n] in site_aa[sites[n]])
                push!(site_aa[sites[n]], cpair_mut.aa_c1[n])
            end
            if !(cpair_mut.aa_c2[n] in site_aa[sites[n]])
                push!(site_aa[sites[n]], cpair_mut.aa_c2[n])
            end
        end

    end

    num_j = 0
    for k in keys(site_aa)
        num_j += length(site_aa[k])
    end

    return site_aa, num_j

end

# Given a mutation it provides the list of unique clade pairs
function unique_clade_pairs(dfit_mut::DataFrame)

    clade1 = dfit_mut[:, :clade1]
    clade2 = dfit_mut[:, :clade2]

    clade_pairs = [Vector{String}(undef, 2) for k in eachindex(clade1)]

    for k in eachindex(clade1)
        clade_pairs[k] = [clade1[k], clade2[k]]
    end

    return clade_pairs

end

# Given a list of clade pairs computes a table with corresponding mismatches
function clade_pair_mut(clade_pairs::Vector{Vector{String}}, clade_diff::DataFrame)

    df_clade_pair_mut = DataFrame(site=Int[], clade1=String[], aa_c1=String[], clade2=String[], aa_c2=String[])
    idx_c1c2 = zeros(Bool, size(clade_diff, 1))

    for c in clade_pairs
        idx_c1c2 .= (clade_diff.clade1 .== c[1]) .& (clade_diff.clade2 .== c[2])
        cdiff_c1c2 = clade_diff[idx_c1c2, 2:end]
        df_clade_pair_mut = vcat(df_clade_pair_mut, cdiff_c1c2)
    end

    return df_clade_pair_mut

end

# Computing the set of weights for a mutation and the associated clade pairs
function weight(dfit_mut::DataFrame; weighted::Bool=true, norm::Bool=true)

    Nc = size(dfit_mut, 1)
    w = zeros(Float64, Nc)

    if weighted
        for n in 1:Nc
            s1 = dfit_mut[n, :std_fit1] * dfit_mut[n, :std_fit1]
            s2 = dfit_mut[n, :std_fit2] * dfit_mut[n, :std_fit2]
            w[n] = 1 / (s1 + s2)
        end
    else
        w .= fill(1.0, Nc)
    end

    if norm
        w ./= sum(w)
    end

    return w

end

# Computing weights for each clade pair associated to a mutation
# according to posterior uncertainty w ∝ 1 / σ²
function weight(s1::Vector{Float64}, s2::Vector{Float64}; weighted::Bool=true, norm::Bool=true)

    @assert length(s1) == length(s2)
    Nc = length(s1)
    w = zeros(Float64, Nc)

    if weighted
        for n in 1:Nc
            s = s1[n] * s1[n] + s2[n] * s2[n]
            w[n] = 1 / s
        end
    else
        w .= fill(1.0, Nc)
    end

    if norm
        w ./= sum(w)
    end

    return w

end

# Function to compute the coefficients of the quadratic form in a vectorized way
function qf_coeff(data::Data{S,T}, opt_x::OptX{Q,T,S}) where {Q<:Int,T<:Real,S<:AbstractString}

    sites = sort(collect(keys(opt_x.site_aa))) # sorted list of sites of background mismatches

    # Variables for vectorized implementation
    aa = vcat([opt_x.site_aa[s] for s in sites]...) # amino acids list
    ss = vcat([fill(s, length(opt_x.site_aa[s])) for s in sites]...) # sites list
    df_cp_mut = [data.cpair_mut[(data.cpair_mut.clade1.==c[1]).&(data.cpair_mut.clade2.==c[2]), [:site, :aa_c1, :aa_c2]] for c in data.clade_pairs] # mismatches for each clade pair

    w_vec = zeros(T, opt_x.num_j, size(df_cp_mut, 1)) # matrix of coefficients

    w_vec .= w_mat.(aa, ss, permutedims(df_cp_mut)) .* reshape(sqrt.(data.weights), 1, :) # weighted coefficients

    return w_vec

end

# Function to compute the constant term of the quadratic form
function offset(dfit_mut::DataFrame, weights::Vector{Float64})

    Δf = dfit_mut[:, :fit2] .- dfit_mut[:, :fit1]
    return sum(weights .* (Δf) .^ 2)

end

# Elements of the matrix wᵃᵇₖ(σₖ)=δ(σₖ,σᵇₖ)-δ(σₖ,σᵃₖ)
function w_mat(a::String, site::Int, cpair_mut_c1c2::DataFrame)

    if site in cpair_mut_c1c2.site
        if a == cpair_mut_c1c2[cpair_mut_c1c2.site.==site, :aa_c1][1]
            return -1
        elseif a == cpair_mut_c1c2[cpair_mut_c1c2.site.==site, :aa_c2][1]
            return +1
        else
            return 0
        end
    else
        return 0
    end

end