function robust_trust_region_(df::OnceDifferentiable,
                              initial_x::AbstractArray{T},
                              xtol::Real,
                              ftol::Real,
                              iterations::Integer,
                              store_trace::Bool,
                              show_trace::Bool,
                              extended_trace::Bool,
                              factor::Real,
                              autoscale::Bool,
                              cache = NewtonTrustRegionCache(df)) where T
    copyto!(cache.x, initial_x)
    value_jacobian!!(df, cache.x)
    cache.r .= value(df)
    check_isfinite(cache.r)

    it = 0
    x_converged, f_converged = assess_convergence(initial_x, cache.xold, value(df), NaN, ftol)
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

        return SolverResults(name,
        #initial_x, reshape(cache.x, size(initial_x)...), norm(cache.r, Inf),
        initial_x, copy(cache.x), norm(cache.r, Inf),
        it, x_converged, xtol, f_converged, ftol, tr,
        first(df.f_calls), first(df.df_calls))
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
            dogleg!(cache.p, cache.p_c, cache.pi, cache.r, cache.d, jacobian(df), delta)
            copyto!(cache.xold, cache.x)
            cache.x .+= cache.p
            value!(df, cache.x)
            # Ratio of actual to predicted reduction (equation 11.47 in N&W)
            mul!(vec(cache.r_predict), jacobian(df), vec(cache.p))
            cache.r_predict .+= cache.r
            rho = (sum(abs2, cache.r) - sum(abs2, value(df))) / (sum(abs2, cache.r) - sum(abs2, cache.r_predict))

            if rho > eta
                # Successful iteration
                cache.r .= value(df)
                jacobian!(df, cache.x)

                # Update scaling vector
                if autoscale
                    for j = 1:nn
                        cache.d[j] = max(convert(real(T), 0.1) * real(cache.d[j]), norm(view(jacobian(df), :, j)))
                    end
                end

                x_converged, f_converged = assess_convergence(cache.x, cache.xold, cache.r, xtol, ftol)
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
            delta = delta/2
            continue
        end
        # Update size of trust region
        if rho < 0.1
            delta = delta/2
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
    return SolverResults(name,
                         initial_x, copy(cache.x), maximum(abs, cache.r),
                         it, x_converged, xtol, f_converged, ftol, tr,
                         first(df.f_calls), first(df.df_calls))
end

function robust_trust_region(df::OnceDifferentiable,
                      initial_x::AbstractArray{T},
                      xtol::Real,
                      ftol::Real,
                      iterations::Integer,
                      store_trace::Bool,
                      show_trace::Bool,
                      extended_trace::Bool,
                      factor::Real,
                      autoscale::Bool,
                      cache = NewtonTrustRegionCache(df)) where T
    robust_trust_region_(df, initial_x, convert(real(T), xtol), convert(real(T), ftol), iterations, store_trace, show_trace, extended_trace, convert(real(T), factor), autoscale, cache)
end
