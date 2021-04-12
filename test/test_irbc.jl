using Test
include("../src/perfectforesight/gmres_solver.jl")
include("../src/perfectforesight/perfectforesight_solvers.jl")
#context = @dynare "../test/models/example1/example1.mod"
context = @dynare "test/models/irbc/irbc1.mod" "savemacro" "-DN=2"

function makeA(jacobian::AbstractMatrix{Float64},
               g::AbstractMatrix{Float64},
               n::Int64)
    i, j, v = findnz(jacobian)
    nvar = size(jacobian, 1)
    m = length(i)
    nm = n*m - count(j .<= nvar) - count(j .> 2*nvar)
    i1 = zeros(Int64, nm)
    j1 = zeros(Int64, nm)
    v1 = zeros(nm)
    r = 1
    for el = 1:m
        if j[el] > nvar
            i1[r] = i[el]
            j1[r] = j[el] - nvar
            v1[r] = v[el]
            r += 1
        end
    end        
    offset = nvar
    for k = 2:(n - 1)
        for el = 1:m
            i1[r] = i[el] + offset
            j1[r] = j[el] + offset - nvar
            v1[r] = v[el]
            r += 1
        end
        offset += nvar
    end
    for el = 1:m
        if j[el] <= 2*nvar
            i1[r] = i[el] + offset
            j1[r] = j[el] + offset - nvar
            v1[r] = v[el]
            r += 1
        end
    end
    A = sparse(i1, j1, v1, n*nvar, n*nvar)
    # terminal condition
    k = (n-1)*nvar + 1:n*nvar
    A[k,k] .+= jacobian[:,2*nvar+1:end]*g 
    return A
end

md = context.models[1]
endo_steadystate = context.results.model_results[1].trends.endogenous_steady_state
exo_steadystate = context.results.model_results[1].trends.exogenous_steady_state

periods = 2;
preconditioner_window = 2
res = zeros(periods*md.endogenous_nbr)
res[[4, 3]] .= 0.001
rout = zeros(periods*md.endogenous_nbr)
ws = GmresWs(periods, preconditioner_window, context, "CR")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= repeat(exo_steadystate', periods + 2)

work = context.work
get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
jacobian = Matrix{Float64}(undef, size(work.jacobian))
copy!(jacobian, work.jacobian)
get_abc!(ws.a, ws.b, ws.c, jacobian, md)
P = LREprecond(periods, preconditioner_window, ws.a, ws.b, ws.c, ws.g)

n = md.endogenous_nbr
time_vec = zeros(md.endogenous_nbr)
get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
LREMap = LinearMap(periods*n) do C, B
    @inbounds copyto!(ws.residuals, n + 1, B, 1, periods*n) 
    jacobian_time_vec!(C, ws.dynamic_variables, ws.residuals, ws.endogenous, ws.exogenous,
                       ws.steadystate, ws.presiduals, ws.g, time_vec, context.work, md, periods)
end
A = makeA(hcat(ws.a, ws.b, ws.c), ws.g, 2)
x = randn(n*periods)
y = P\x
@test LREMap*y ≈ x
@test A*x ≈ LREMap*x
target = [P.hh*x; ws.g*P.hh*x + P.hh[:, 1:n]*x[14:26]]
@test P\x ≈ target
x = zeros(n*periods)
x[[4,6]] .= 0.1
y = P\x
@test LREMap*y ≈ x

periods = 10
preconditioner_windows = 2
res = zeros(periods*md.endogenous_nbr)
res[[4, 6]] .= 0.1
rout = zeros(periods*md.endogenous_nbr)
ws = GmresWs(periods, preconditioner_window, context, "CR")

get_abc!(ws.a, ws.b, ws.c, jacobian, md)
P1 = LREprecond(periods, preconditioner_window, ws.a, ws.b, ws.c, ws.g)
get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
LREMap = LinearMap(periods*n) do C, B
    @inbounds copyto!(ws.residuals, n + 1, B, 1, periods*n) 
    jacobian_time_vec!(C, ws.dynamic_variables, ws.residuals, ws.endogenous, ws.exogenous,
                       ws.steadystate, ws.presiduals, ws.g, time_vec, context.work, md, periods)
end
y1 = P1\res

ws = GmresWs(periods, preconditioner_window, context, "CR")
ws.endogenous .= repeat(endo_steadystate, periods + 2)
ws.exogenous .= zeros(periods + 2, md.exogenous_nbr)
gmres!(rout, ws.LREMap, res, log=false,
       verbose=false, Pr=ws.P)
y = P\rout
@test ws.LREMap*y ≈ res

get_jacobian!(work, ws.endogenous, ws.exogenous, endo_steadystate, md, 2)
jacobian = Matrix{Float64}(undef, size(work.jacobian))
copy!(jacobian, work.jacobian)
get_abc!(ws.a, ws.b, ws.c, jacobian, md)
A = makeA(hcat(ws.a, ws.b, ws.c), ws.g, periods)
display([LREMap*rout A*rout res])
display([y1 rout])




