module DynareSolvers

##
 # Copyright (C) 2021 Dynare Team
 #
 # This file is part of Dynare.
 #
 # Dynare is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # Dynare is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
##

export trustregion, TrustRegionWA

using LinearAlgebra

const p1 = .1
const p5 = .5
const p001 = .001
const p0001 = .0001
const macheps = eps(Float64)

mutable struct TrustRegionWA
    x::Vector{Float64}                   # Vector of unknowns
    xx::Vector{Float64}                  # Vector of unknowns
    fval::Vector{Float64}                # residuals of the non linear equations
    fval0::Vector{Float64}
    fval1::Vector{Float64}               # residuals of the non linear equations (next)
    fjac::Matrix{Float64}                # jacobian matrix of the non linear equations
    fjaccnorm::Vector{Float64}           # norms of the columns of the Jacobian matrix
    fjaccnorm__::Vector{Float64}         # copy of fjaccnorm
    w::Vector{Float64}                   # working array
    p::Vector{Float64}                   # direction
    s::Vector{Float64}
end

function TrustRegionWA(n::Int)
    TrustRegionWA(Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n),
                  Vector{Float64}(undef,n), Matrix{Float64}(undef,n,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n),
                  Vector{Float64}(undef,n), Vector{Float64}(undef,n), Vector{Float64}(undef,n))
end


"""
    dogleg!(x::Vector{Float64}, r::Matrix{Float64}, b::Vector{Float64}, d::Vector{Float64}, δ::Float64)

Given an `n` by `n` matrix `r`, an `n` by 1 vector `d` with non zero entries, an `n` by `1`
vector `b`, and a positive number δ, the problem is to determine the convex combination `x`
of the Gauss-Newton and scaled gradient directions that minimizes (r*x - b) in the least
squares sense, subject to the restriction that the Euclidean norm of d*x be at most delta.
"""
function dogleg!(x::Vector{Float64}, r::Matrix{Float64}, b::Vector{Float64}, d::Vector{Float64}, δ::Float64, s::Vector{Float64})
    n = length(x)
    # Compute the Gauss-Newton direction.
    x .= r\b
    # Compute norm of scaled x.
    qnorm = zero(Float64)
    @inbounds for i=1:n
        qnorm += (d[i]*x[i])^2
    end
    qnorm = sqrt(qnorm)
    if qnorm<=δ
        # Gauss-Newton direction is acceptable. There is nothing to do here.
    else
        # Gauss-Newton direction is not acceptable…
        # Compute the scale gradient direction and its norm
        gnorm = zero(Float64)
        @inbounds for i=1:n
            s[i] = zero(Float64)
            @inbounds for j=1:n
                s[i] += r[j,i]*b[j]
            end
            s[i] /= d[i]
            gnorm += s[i]*s[i]
        end
        gnorm = sqrt(gnorm)
        # Compute the norm of the scaled gradient direction.
        # gnorm = norm(s)
        if gnorm>0
            # Normalize and rescale → gradient direction.
            @inbounds for i = 1:n
                s[i] /= gnorm*d[i]
            end
            temp0 = zero(Float64)
            @inbounds for i=1:n
                temp1 = zero(Float64)
                @inbounds for j=1:n
                    temp1 += r[i,j]*s[j]
                end
                temp0 += temp1*temp1
            end
            sgnorm = gnorm/temp0
            temp0 = sqrt(temp0)
            if sgnorm<=δ
                # The scaled gradient direction is not acceptable…
                # Compute the point along the dogleg at which the
                # quadratic is minimized.
                bnorm = norm(b)
                temp1 = δ/qnorm
                temp2 = sgnorm/δ
                temp0 = bnorm*bnorm*temp2/(gnorm*qnorm)
                temp0 -= temp1*temp2*temp2-sqrt((temp0-temp1)^2+(one(Float64)-temp1*temp1)*(one(Float64)-temp2*temp2))
                α = temp1*(one(Float64)-temp2*temp2)/temp0
            else
                # The scaled gradient direction is acceptable.
                α = zero(Float64)
            end
        else
            # If the norm of the scaled gradient direction is zero.
            α = δ/qnorm
            sgnorm = .0
        end
        # Form the appropriate  convex combination of the Gauss-Newton direction and the
        # scaled gradient direction.
        temp1 = (one(Float64)-α)*min(sgnorm, δ)
        @inbounds for i = 1:n
            x[i] = α*x[i] + temp1*s[i]
        end
    end
end

"""
    trustregion(f!::Function, j!::Function, y0::Vector{Float64}, tolf::Float64, tolx::Float64, maxiter::Int)

Solves a system of nonlinear equations of several variables using a trust region method.

# Arguments
- `f!::Function`: Function evaluating the residuals of the nonlinear equations to be solved. This function admits two arguments, a vector holding the values of the residuals and the vector of unknowns (to be solved for).
- `j!::Function`: Function computing the jacobian of the nonlinear equations (derivatives with respect to the unknowns). This function admits two arguments: the jacobian matrix and the vector of unknowns.
- `x0::Vector{Float64}`: Initial guess for the unknowns.
- `tolf::Float64`: Positive real scalar used in determining the terminination criterion (test on the values of the residuals). Default is 1e-6.
- `tolx::Float64`: Positive real scalar used in determining the terminination criterion (test on the change of the unknowns). Default is 1e-6.
- `maxiter::Int`: Positive integer scalar, maximum number of iterations. Default is 50.
"""
function trustregion(f!::Function, j!::Function, x0::Vector{Float64}, tolx::Float64=1e-6, tolf::Float64=1e-6, maxiter::Int=50)
    wa = TrustRegionWA(length(x0))
    info = trustregion(f!, j!, x0, 1.0, tolx, tolf, maxiter, wa)
    return wa.x, info
end

"""
    trustregion(f!::Function, j!::Function, y0::Vector{Float64}, factor::Float64, tolf::Float64, tolx::Float64, maxiter::Int, wa::TrustRegionWA)

Solves a system of nonlinear equations of several variables using a trust region method. This version requires a last argument of type `TrustRegionWA`
which holds various working arrays used in the function (avoiding array instantiations). Results of the solver are available in `wa.x` if `info=1`.

# Arguments
- `f!::Function`: Function evaluating the residuals of the nonlinear equations to be solved. This function admits two arguments, a vector holding the values of the residuals and the vector of unknowns (to be solved for).
- `j!::Function`: Function computing the jacobian of the nonlinear equations (derivatives with respect to the unknowns). This function admits two arguments: the jacobian matrix and the vector of unknowns.
- `x0::Vector{Float64}`: Initial guess for the unknowns.
- `factor::Float64`: Positive real scalar used in determining the initial step bound. In most cases factor should lie in the interval (.1,100.).
- `tolf::Float64`: Positive real scalar used in determining the terminination criterion (test on the values of the residuals).
- `tolx::Float64`: Positive real scalar used in determining the terminination criterion (test on the change of the unknowns).
- `maxiter::Int`: Positive integer scalar, maximum number of iterations.
- `wa::TrustRegionWA`: Working arrays.
"""
function trustregion(f!::Function, j!::Function, x0::Vector{Float64}, factor::Float64, tolx::Float64, tolf::Float64, maxiter::Int, wa::TrustRegionWA)
    n, iter, info = length(x0), 1, 0
    xnorm, xnorm0 = one(Float64), one(Float64)
    fnorm, fnorm1, fnorm0 = one(Float64), one(Float64), one(Float64)
    wa.x .= x0
    # Initial evaluation of the residuals (and compute the norm of the residuals)
    try
        f!(wa.fval, wa.x)
        fnorm = norm(wa.fval)
    catch
        error("The system of nonlinear equations cannot be evaluated on the initial guess!")
    end
    # Initial evaluation of the Jacobian
    try
        j!(wa.fjac, wa.x)
    catch
        error("The Jacobian of the system of nonlinear equations cannot be evaluated on the initial guess!")
    end
    # Initialize counters.
    ncsucc, nslow1= zero(Int), zero(Int)
    # Initialize scale parameter.
    scale, scale0 = one(Float64), one(Float64)
    # Newton iterations
    δ = 0.0
    while iter<=maxiter && info==0
        # Compute columns norm for the Jacobian matrix.
        @inbounds for i=1:n
            wa.fjaccnorm[i] = zero(Float64)
            @inbounds for j = 1:n
                wa.fjaccnorm[i] += wa.fjac[j,i]*wa.fjac[j,i]
            end
            wa.fjaccnorm[i] = sqrt(wa.fjaccnorm[i])
        end
        if iter==1
            # On the first iteration, calculate the norm of the scaled vector of unknwonws x
            # and initialize the step bound δ. Scaling is done according to the norms of the
            # columns of the initial jacobian.
            @inbounds for i = 1:n
                wa.fjaccnorm__[i] = abs(wa.fjaccnorm[i] < macheps) ? 1.0 : wa.fjaccnorm[i]
            end
            xnorm = zero(Float64)
            @inbounds for i=1:n
                xnorm += (wa.fjaccnorm__[i]*wa.x[i])^2
            end
            δ = sqrt(xnorm)
            if δ<macheps
                δ = one(Float64)
            end
            δ *= factor
        else
            wa.fjaccnorm__ .= max.(0.1.*wa.fjaccnorm__, wa.fjaccnorm)
        end
        # Determine the direction p (with trust region model defined in dogleg routine).
        dogleg!(wa.p, wa.fjac, wa.fval, wa.fjaccnorm__, δ, wa.s)
        # Compute the norm of p.
        pnorm = zero(Float64)
        @inbounds for i=1:n
            pnorm += (wa.fjaccnorm__[i]*wa.p[i])^2
        end
        pnorm = sqrt(pnorm)
        # On first iteration adjust the initial step bound.
        if iter==1
            δ = min(δ, pnorm)
        end
        fwrong, jwrong, siter = true, true, 0
        while (fwrong || jwrong) && scale>.0005
            # Move along the direction p. Set a candidate value for x and predicted improvement for f.
            @inbounds for i=1:n
                wa.xx[i] = wa.x[i]-scale*wa.p[i]
                wa.w[i] = wa.fval[i]
                @inbounds for j = 1:n
                    wa.w[i] -= scale*wa.fjac[i,j]*wa.p[j]
                end
            end
            # Evaluate the function at xx = x+p and calculate its norm.
            try
                f!(wa.fval1, wa.xx)
                fwrong = false
            catch
                # If evaluation of the residuals returns an error, then keep the same
                # direction but reduce the step length.
                scale *= .5
                siter += 1
                continue
            end
            fnorm1 = norm(wa.fval1)
            # Compute the scaled actual reduction.
            if fnorm1<fnorm
                actualreduction = one(Float64)-(fnorm1/fnorm)^2
            else
                actualreduction = -one(Float64)
            end
            # Compute the scaled predicted reduction and the ratio of the actual to the
            # predicted reduction.
            t = norm(wa.w)
            ratio = zero(Float64)
            if t<fnorm
                predictedreduction = one(Float64) - (t/fnorm)^2
                ratio = actualreduction/predictedreduction
            end
            # Update the radius of the trust region if need be.
            δ0 = δ
            ncsucc0 = ncsucc
            if ratio<p1
                # Reduction is much smaller than predicted… Reduce the radius of the trust region.
                ncsucc = 0
                δ *= p5
            else
                ncsucc += 1
                if ratio>=p5 || ncsucc>1
                    δ = max(δ,pnorm/p5)
                end
                if abs(ratio-one(Float64))<p1
                    δ = pnorm/p5
                end
            end
            xnorm0 = xnorm
            fnorm0 = fnorm
            @inbounds for i=1:n
                wa.s[i] = wa.x[i]
                wa.fval0[i] = wa.fval[i]
            end
            if ratio>=1.0e-4
                # Succesfull iteration. Update x, xnorm, fval, fnorm and fjac.
                xnorm = zero(Float64)
                @inbounds for i=1:n
                    # Update of x
                    wa.x[i] = wa.xx[i]
                    xnorm += (wa.fjaccnorm__[i]*wa.x[i])^2
                    # Update fval
                    wa.fval[i] = wa.fval1[i]
                end
                xnorm = sqrt(xnorm)
                fnorm = fnorm1
            end
            # Determine the progress of the iteration.
            nslow0 = nslow1
            nslow1 += 1
            if actualreduction>=p001
                nslow1 = 0
            end
            # Test for convergence.
            if δ<tolx*xnorm || fnorm<tolf
                info = 1
                @goto mainloop
            end
            # Tests for termination and stringent tolerances.
            if p1*max(p1*δ, pnorm)<=macheps*xnorm
                # xtol is too small. no further improvement in
                # the approximate solution x is possible.
                info = 3
                @goto mainloop
            end
            if nslow1==15
                # iteration is not making good progress, as
                # measured by the improvement from the last
                # fifteen iterations.
                info = 4
                @goto mainloop
            end
            # Compute the jacobian for next the iteration.
            try
                j!(wa.fjac, wa.x)
                jwrong = false
                iter += 1
            catch
                # If evaluation of the jacobian returns an error, then keep the same
                # direction but reduce the step length.
                xnorm = xnorm0
                fnorm = fnorm0
                δ = δ0
                ncsucc = ncsucc0
                nslow1 = nslow0
                @inbounds for i=1:n
                    wa.x[i] = wa.s[i]
                    wa.fval[i] = wa.fval0[i]
                end
                scale *= .5
                siter += 1
                jwrong = true
            end
            if fwrong || jwrong
                info = 5
                fill!(wa.x, Inf)
                return info
            end
        end
        # Update the value of the scale parameter.
        if siter>0
            # Something went wrong when evaluating the nonlinear equations or the
            # jacobian matrix, and the scale parameter had to be reduced. The scale
            # parameter is updated with its average across newton iterations (first
            # while-loop). This avoids to use the default value of the scale
            # parameter (1.0) in the following iteration and reduces the number of
            # iterations. The average value of the scale parameter is recursively
            # computed.
            scale = ((iter-1)*scale0+scale)/iter
        else
            # Increase the value of the scale parameter by 5 percent if the previous
            # step provided by the dogleg routine did not cause any trouble...
            scale = min(scale0*1.05, 1.0)
        end
        scale0 = scale
        @label mainloop
    end
    if info==0 && iter>maxiter
        info = 2
        fill!(wa.x, Inf)
        return info
    end
    return info
end

end
