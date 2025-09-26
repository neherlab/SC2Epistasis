""" In this file we define code and functions to
determine the coupling parameters  according to
numerical optimization methods """

export learn

"""
    learn(optx, qform, λ_vec; kwargs...)

Infer the coupling parameters through OWL-QN optimization from structs `OptX` and `QForm`.

# Arguments
- `optx::Vector{OptX}`: vector of `OptX` structs,
- `qform::Vector{QForm}`: vector of `QForm` structs,
- `λ_vec::Vector{Float64}`: vector of L1 regularization strengths
- `λ2::Float64=1.0e-6`: L2 regularization strength,
- `epsconv::Float64=1.0e-8`: Convergence threshold,
- `verbose::Bool=false`: Verbosity flag,
- `maxiter::Int=5000`: Maximum number of iterations.

"""
function learn(optx::Vector{OptX}, qform::Vector{QForm}, λ_vec::Vector{Float64};
    λ2::Float64=1.0e-6,
    epsconv::Float64=1.0e-8,
    verbose::Bool=false,
    maxiter::Int=5000)

    @assert length(optx) == length(qform)
    M = length(optx)

    start = 0
    idx_j = fill(1:1, M)

    # Indexes for coupling parameters in the big vector
    for m = 1:M
        idx_j[m] = start+1:(start+optx[m].num_j)
        start += optx[m].num_j
    end

    @assert start == length(λ_vec)

    # Unwrap coupling parameters into a big vector
    J = vcat([optx[m].J for m = 1:M]...)

    # Initialize optimizer struct, gradient and loss
    owl = OWLQN(J, λ_vec)
    g = zeros(Float64, start)
    l = 0.0
    loss = 0.0

    # Pre-computation for first update step
    grad!(g, J, qform, idx_j, λ2)
    ener = energy(J, qform, idx_j, λ2)

    # Main optimization loop
    for i in 1:maxiter

        # Update step for state vector and gradient
        step!(owl, ener, g, J, qform, idx_j, λ2)


        ener = energy(J, qform, idx_j, λ2) # energy contribution
        nrm = sum(abs.(J) .* λ_vec) # L1 norm contribution
        loss = ener + nrm # total loss

        # Convergence check
        if maximum(owl.s[end]) <= epsconv
            verbose && println("status: XTOL_REACHED,\tloss: $loss\n")
            break
        elseif abs(loss - l) <= epsconv
            verbose && println("status: FTOL_REACHED,\tloss: $loss\n")
            break
        end
        l = loss
    end

    # Write coupling into dictionary
    Threads.@threads for m = 1:M
        optx[m].J .= J[idx_j[m]]
        counter = 0
        for s in sort(collect(keys(optx[m].site_aa)))
            for σ in optx[m].site_aa[s]
                counter += 1
                push!(optx[m].epi, [s, σ] => optx[m].J[counter])
            end
        end
    end

    return loss

end

# single update of parameters x
function step!(M::OWLQN, ener::Float64, g::Vector{Float64}, x::Vector{Float64}, θ::Vector{QForm}, idx_j::Vector{UnitRange{Int64}}, λ::Float64)
    s, y, rho, lambda = M.s, M.y, M.rho, M.lambda

    if all(g .== 0.0)
        println("(Local) minimum found: ∇f(x) == 0")
        return
    end

    m = min(M.m, size(s)[1])
    M.t += 1

    x_copy = deepcopy(x)
    g_copy = deepcopy(g)

    pg = pseudo_gradient(g, x, lambda)
    Q = deepcopy(pg)

    if m > 0

        # L-BFGS computation of Hessian-scaled gradient z = H_inv * g
        alpha = []
        for i in m:-1:1
            push!(alpha, rho[i] * (s[i] ⋅ Q))
            Q -= alpha[end] * y[i]
        end
        reverse!(alpha)
        z = Q .* (s[end] ⋅ y[end]) / (y[end] ⋅ y[end])
        for i in 1:m
            z += s[i] * (alpha[i] - rho[i] * (y[i] ⋅ z))
        end

        # zeroing out all elements in z if sign(z[i]) != sign(g[i])
        # that is, if scaling changes gradient sign
        project!(z, pg)

        # fancy way to do x .-= z
        projected_backtracking_line_search_update(ener, pg, x, z, lambda, θ, idx_j, λ)
    else
        # fancy way to do  x .-= g
        projected_backtracking_line_search_update(ener, pg, x, pg, lambda, θ, idx_j, λ; alpha=1 / sqrt(pg ⋅ pg))
    end

    push!(s, x - x_copy)
    push!(y, grad!(g, x, θ, idx_j, λ) - g_copy)
    push!(rho, 1 / (y[end] ⋅ s[end]))

    while length(s) > M.m
        popfirst!(s)
        popfirst!(y)
        popfirst!(rho)
    end
end

# projected gradient based on raw gradient, parameter values, and L1 reg. strength
function pseudo_gradient(g::Vector{Float64}, x::Vector{Float64}, lambda::Vector{Float64})
    pg = zeros(size(g))
    for i in 1:size(g)[1]
        if x[i] > 0
            pg[i] = g[i] + lambda[i]
        elseif x[i] < 0
            pg[i] = g[i] - lambda[i]
        else
            if g[i] + lambda[i] < 0
                pg[i] = g[i] + lambda[i]
            elseif g[i] - lambda[i] > 0
                pg[i] = g[i] - lambda[i]
            end
        end
    end
    return pg
end

# pi alignment operator - projection of a on orthat defined by b
function project!(a::Vector{Float64}, b::Vector{Float64})
    for i in 1:size(a)[1]
        if sign(a[i]) != sign(b[i])
            a[i] = 0.0
        end
    end
end

# projected backtracking line search
function projected_backtracking_line_search_update(ener::Float64, pg::Vector{Float64}, x::Vector{Float64}, z::Vector{Float64}, lambda::Vector{Float64}, θ::Vector{QForm}, idx_j::Vector{UnitRange{Int64}}, λ::Float64;
    alpha::Float64=1.0,
    beta::Float64=0.5,
    gamma::Float64=1e-4)

    y = ener + sum(abs.(x) .* lambda)

    # choose orthant for the new point
    xi = sign.(x)
    for i in 1:size(xi)[1]
        if xi[i] == 0
            xi[i] = sign(-pg[i])
        end
    end

    while true
        # update current point
        xt = x - alpha * z

        # project point onto orthant
        project!(xt, xi)

        # sufficient decrease condition
        if energy(xt, θ, idx_j, λ) + sum(abs.(xt) .* lambda) <= y + gamma * (pg ⋅ (xt - x))
            x .= xt
            break
        end

        # update step size
        alpha *= beta
    end
end

# Compute the energy contribution to the loss function
function energy(x::Vector{Float64}, θ::Vector{QForm}, idx_j::Vector{UnitRange{Int64}}, λ::Float64)

    ener = 0.0
    for m in eachindex(idx_j)
        ener += dot(x[idx_j[m]], θ[m].A, x[idx_j[m]]) + 2 * dot(θ[m].f, x[idx_j[m]]) + θ[m].c + λ * sum(x[idx_j[m]] .^ 2)
    end

    return ener

end

# Compute and update in-place the loss gradient
function grad!(g::Vector{Float64}, x::Vector{Float64}, θ::Vector{QForm}, idx_j::Vector{UnitRange{Int64}}, λ::Float64)

    for m in eachindex(idx_j)
        g[idx_j[m]] .= 2 * (θ[m].A * x[idx_j[m]] .+ θ[m].f .+ λ * x[idx_j[m]])
    end

    return g

end

##############################--------------------------------

# Definition of L1 regularization

# Regularization proportional to the 3D distance between residues
function threeD_l1(λ::Float64, optx::Vector{OptX}, data::Vector{Data}, dist::Vector{Float64})

    reg = threeD_l1(λ, dist)

    # Mutation-wise coupling indices
    M = length(optx)
    start = 0
    idx_j = fill(1:1, M)
    for m = 1:M
        idx_j[m] = start+1:(start+optx[m].num_j)
        start += optx[m].num_j
    end

    idx_jstop = vcat(collect.(idx_j[findall(x -> x[end] == '*', [data[m].dfit_mut.aa_mut[1] for m in 1:M])])...)

    reg[idx_jstop] .= λ # no distance-dependent regularization for couplings to stop codon

    return reg

end

# Regularization proportional to the 3D distance between residues
function threeD_l1(λ::Float64, dist::Vector{Float64})

    dmax = maximum(dist)
    reg = copy(dist) ./ dmax
    for i in eachindex(dist)
        reg[i] = reg[i] * λ
    end

    return reg

end

# Regularization defined by non-linear sigmoid function of the 3D distance between residues
# λ1 and λ2 are the regularization values at the limits of the sigmoid function
# α is the steepness of the sigmoid function
# d0 is the distance at which the regularization is halfway between λ1 and λ2
function sigmoid_l1(λ_low::Float64, λ_high::Float64, optx::Vector{OptX}, data::Vector{Data}, dist::Vector{Float64};
    α::Float64=1.0,
    d0::Float64=15.0)

    reg = sigmoid_l1(λ_low, λ_high, dist; α=α, d0=d0)

    # Mutation-wise coupling indices
    M = length(optx)
    start = 0
    idx_j = fill(1:1, M)
    for m = 1:M
        idx_j[m] = start+1:(start+optx[m].num_j)
        start += optx[m].num_j
    end

    idx_jstop = vcat(collect.(idx_j[findall(x -> x[end] == '*', [data[m].dfit_mut.aa_mut[1] for m in 1:M])])...)

    reg[idx_jstop] .= λ_high # no distance-dependent regularization for couplings to stop codon

    return reg

end

function sigmoid_l1(λ_low::Float64, λ_high::Float64, dist::Vector{Float64};
    α::Float64=1.0,
    d0::Float64=15.0)

    reg = copy(dist)
    for i in eachindex(dist)
        reg[i] = (λ_high - λ_low) / (1 + exp(-α * (dist[i] - d0))) + λ_low
    end

    return reg

end