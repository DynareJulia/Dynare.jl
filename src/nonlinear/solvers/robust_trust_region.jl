struct NewtonTrustRegion end

struct NewtonTrustRegionCache{Tx} <: AbstractSolverCache
    x::Tx
    xold::Tx
    r::Tx
    r_predict::Tx
    p::Tx
    p_c::Tx
    pi::Tx
    d::Tx
end
function NewtonTrustRegionCache(df)
    x = copy(df.x_f) # Current point
    xold = copy(x) # Old point
    r = copy(df.F)       # Current residual
    r_predict = copy(x)  # predicted residual
    p = copy(x)          # Step
    p_c = copy(x)          # Cauchy point
    pi = copy(x)
    d = copy(x)          # Scaling vector
    NewtonTrustRegionCache(x, xold, r, r_predict, p, p_c, pi, d)
end
macro trustregiontrace(stepnorm)
    esc(
        quote
            if tracing
                dt = Dict()
                if extended_trace
                    dt["x"] = copy(cache.x)
                    dt["f(x)"] = copy(value(df))
                    dt["g(x)"] = copy(jacobian(df))
                    dt["delta"] = delta
                    dt["rho"] = rho
                end
                update!(
                    tr,
                    it,
                    maximum(abs, value(df)),
                    $stepnorm,
                    dt,
                    store_trace,
                    show_trace,
                )
            end
        end,
    )
end

function dogleg!(p, p_c, p_i, r, d, J, linsolve, delta::Real)
    T = eltype(d)
    #Jbal = copy(J)
    try
        #=
        @debug "$(now()): start Jbal"
        nn = length(r)
        for i = 1: nn
            col = view(Jbal,:,i)
            col .*= 1/d[i]
        end
        @debug "$(now()): end Jbal"
        =#
        #        copyto!(p_i, J \ vec(r)) # Gauss-Newton step
        #copyto!(p_i, (J \ vec(r)) ./ d)
        # F = lu(J)
        # ldiv!(p_i, F, vec(r))
        linsolve(p_i, J, vec(r))
    catch e
        if isa(e, LAPACKException) || isa(e, SingularException)
            # If Jacobian is singular, compute a least-squares solution to J*x+r=0
            Jfull = convert(Matrix{T}, Jbal)
            @show cond(Matrix(J))
            @show cond(Jfull)
            U, S, V = svd(Jfull) # Convert to full matrix because sparse SVD not implemented as of Julia 0.3
            k = sum(S .> eps())
            mrinv =
                V * Matrix(Diagonal([1 ./ S[1:k]; zeros(eltype(S), length(S) - k)])) * U' # Moore-Penrose generalized inverse of J
            vecpi = vec(p_i)
            mul!(vecpi, mrinv, vec(r))
            vecpi ./= d
        else
            @show "other $e"
            throw(e)
        end
    end
    rmul!(p_i, -one(T))

    # Test if Gauss-Newton step is within the region
    if wnorm(d, p_i) <= delta
        copyto!(p, p_i)   # accepts equation 4.13 from N&W for this step
    else
        # For intermediate we will use the output array `p` as a buffer to hold
        # the gradient. To make it easy to remember which variable that array
        # is representing we make g an alias to p and use g when we want the
        # gradient

        # compute g = J'r ./ (d .^ 2)
        g = p
        mul!(vec(g), transpose(J), vec(r))
        g .= g ./ d .^ 2

        # compute Cauchy point
        p_c .= -wnorm(d, g)^2 / sum(abs2, J * vec(g)) .* g

        if wnorm(d, p_c) >= delta
            # Cauchy point is out of the region, take the largest step along
            # gradient direction
            rmul!(g, -delta / wnorm(d, g))

            # now we want to set p = g, but that is already true, so we're done

        else
            # from this point on we will only need p_i in the term p_i-p_c.
            # so we reuse the vector p_i by computing p_i = p_i - p_c and then
            # just so we aren't confused we name that p_diff
            p_i .-= p_c
            p_diff = p_i

            # Compute the optimal point on dogleg path
            b = 2 * wdot(d, p_c, d, p_diff)
            a = wnorm(d, p_diff)^2
            tau = (-b + sqrt(b^2 - 4a * (wnorm(d, p_c)^2 - delta^2))) / (2a)
            p_c .+= tau .* p_diff
            copyto!(p, p_c)
        end
    end
end

function robust_trust_region_(
    df::OnceDifferentiable,
    initial_x::AbstractArray{T},
    xtol::Real,
    ftol::Real,
    iterations::Integer,
    store_trace::Bool,
    show_trace::Bool,
    extended_trace::Bool,
    factor::Real,
    autoscale::Bool,
    linsolve,
    cache = NewtonTrustRegionCache(df),
) where {T}
    copyto!(cache.x, initial_x)
    value_jacobian!!(df, cache.x)
    cache.r .= value(df)
    check_isfinite(cache.r)

    it = 0
    x_converged, f_converged =
        assess_convergence(initial_x, cache.xold, value(df), NaN, ftol)
    stopped = any(isnan, cache.x) || any(isnan, value(df)) ? true : false

    converged = x_converged || f_converged
    delta = convert(real(T), NaN)
    rho = convert(real(T), NaN)
    if converged
        tr = SolverTrace()
        name = "Trust-region with dogleg"
        if autoscale
            name *= " and autoscaling"
        end

        return SolverResults(
            name,
            #initial_x, reshape(cache.x, size(initial_x)...), norm(cache.r, Inf),
            initial_x,
            copy(cache.x),
            norm(cache.r, Inf),
            it,
            x_converged,
            xtol,
            f_converged,
            ftol,
            tr,
            first(df.f_calls),
            first(df.df_calls),
        )
    end

    tr = SolverTrace()
    tracing = store_trace || show_trace || extended_trace
    @trustregiontrace convert(real(T), NaN)
    nn = length(cache.x)
    if autoscale
        for j = 1:nn
            cache.d[j] = norm(view(jacobian(df), :, j))
            if cache.d[j] == zero(cache.d[j])
                cache.d[j] = one(cache.d[j])
            end
        end
    else
        fill!(cache.d, one(real(T)))
    end

    delta = factor * wnorm(cache.d, cache.x)
    if delta == zero(delta)
        delta = factor
    end

    eta = convert(real(T), 1e-4)

    while !stopped && !converged && it < iterations
        it += 1

        try
            # Compute proposed iteration step
            dogleg!(cache.p, cache.p_c, cache.pi, cache.r, cache.d, jacobian(df), linsolve, delta)
            copyto!(cache.xold, cache.x)
            cache.x .+= cache.p
            value!(df, cache.x)
            # Ratio of actual to predicted reduction (equation 11.47 in N&W)
            mul!(vec(cache.r_predict), jacobian(df), vec(cache.p))
            cache.r_predict .+= cache.r
            rho =
                (sum(abs2, cache.r) - sum(abs2, value(df))) /
                (sum(abs2, cache.r) - sum(abs2, cache.r_predict))

            if rho > eta
                # Successful iteration
                cache.r .= value(df)
                jacobian!(df, cache.x)

                # Update scaling vector
                if autoscale
                    for j = 1:nn
                        cache.d[j] = max(
                            convert(real(T), 0.1) * real(cache.d[j]),
                            norm(view(jacobian(df), :, j)),
                        )
                    end
                end

                x_converged, f_converged =
                    assess_convergence(cache.x, cache.xold, cache.r, xtol, ftol)
                converged = x_converged || f_converged
            else
                cache.x .-= cache.p
                x_converged, converged = false, false

            end

            @trustregiontrace euclidean(cache.x, cache.xold)
            if any(isnan, cache.x) || any(isnan, value(df))
                throw(DomainError(NaN, "model returns NaN"))
            end
        catch e
            cache.x .-= cache.p
            delta = delta / 2
            continue
        end
        # Update size of trust region
        if rho < 0.1
            delta = delta / 2
        elseif rho >= 0.9
            delta = 2 * wnorm(cache.d, cache.p)
        elseif rho >= 0.5
            delta = max(delta, 2 * wnorm(cache.d, cache.p))
        end
    end

    name = "Trust-region with dogleg"
    if autoscale
        name *= " and autoscaling"
    end
    return SolverResults(
        name,
        initial_x,
        copy(cache.x),
        maximum(abs, cache.r),
        it,
        x_converged,
        xtol,
        f_converged,
        ftol,
        tr,
        first(df.f_calls),
        first(df.df_calls),
    )
end

function robust_trust_region(
    df::OnceDifferentiable,
    initial_x::AbstractArray{T},
    xtol::Real,
    ftol::Real,
    iterations::Integer,
    store_trace::Bool,
    show_trace::Bool,
    extended_trace::Bool,
    factor::Real,
    autoscale::Bool;
    linsolve = (x, A, b) -> copyto!(x, A\b),
    cache = NewtonTrustRegionCache(df),
) where {T}
    robust_trust_region_(
        df,
        initial_x,
        convert(real(T), xtol),
        convert(real(T), ftol),
        iterations,
        store_trace,
        show_trace,
        extended_trace,
        convert(real(T), factor),
        autoscale,
        linsolve,
        cache,
    )
end
