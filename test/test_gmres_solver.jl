include("../src/perfectforesight/gmres_solver.jl")
include("../src/perfectforesight/perfectforesight_solvers.jl")
context = @dynare "test/models/example1/example1.mod"
using LinearRationalExpectations
using Test

m = context.models[1]
lli = m.lead_lag_incidence
dynamic_variables = zeros(nnz(sparse(lli)))
steadystate = context.results.model_results[1].trends.endogenous_steady_state
endogenous = repeat(steadystate,3)
exogenous = zeros(3, m.exogenous_nbr)
get_dynamic_endogenous_variables!(dynamic_variables, endogenous,
                                  lli, m, 2)
target = vcat(steadystate[lli[1,:] .> 0],
              steadystate[lli[2,:] .> 0],
              steadystate[lli[3,:] .> 0])
@test  dynamic_variables ≈ target

y = zeros(5,10)
x = randn(5, 3)
k = [1, 5, 10]

fan_columns!(y, x, k, 0)

@test y[:, k] == x


y = zeros(5,10)
x = randn(5, 5)
k = [1, 5, 10]

fan_columns!(y, x, k, 2)
@time fan_columns!(y, x, k, 2)

@test y[:, k] == x[:, 3:5]

struct M;
    i_bkwrd_b
    i_current
    i_fwrd_b
    n_bkwrd
    n_both
    n_current
end

#m = M([1, 3, 6],
#      collect(2:6),
#      [1 , 3, 5],
#      1, 2, 5)
lead_lag_incidence = vcat([1,  0, 2, 0,  0, 3]',
                          [0,  4, 5, 6,  7, 8]',
                          [9, 10, 0, 0, 11, 0]')
m_orig = context.models[1]
m = Model("test/models/example1/example1",
          m_orig.endogenous_nbr,
          lead_lag_incidence,
          m_orig.exogenous_nbr,
          m_orig.lagged_exogenous_nbr,
          m_orig.exogenous_deterministic_nbr,
          m_orig.parameter_nbr,
          m_orig.maximum_endo_lag,
          m_orig.maximum_endo_lead,
          m_orig.maximum_exo_lag,
          m_orig.maximum_exo_lead,
          m_orig.maximum_exo_det_lag,
          m_orig.maximum_exo_det_lead,
          m_orig.maximum_lag,
          m_orig.maximum_lead,
          m_orig.orig_maximum_endo_lag,
          m_orig.orig_maximum_endo_lead,
          m_orig.orig_maximum_exo_lag,
          m_orig.orig_maximum_exo_lead,
          m_orig.orig_maximum_exo_det_lag,
          m_orig.orig_maximum_exo_det_lead,
          m_orig.orig_maximum_lag,
          m_orig.orig_maximum_lead)

n = 6
jacobian = randn(n, 12)
a = zeros(n, n)
b = zeros(n, n)
c = zeros(n, n)

get_abc!(a, b, c, jacobian, m)
@time get_abc!(a, b, c, jacobian, m)
@test a[:, m.i_bkwrd_b] == jacobian[:, 1:3]
@test b[:, m.i_current] == jacobian[:, 4:8]
@test c[:, m.i_fwrd_b]  == jacobian[:, 9:11]

h0 = zeros(n,n)

g = randn(n, n)

ws_linsolve = LinSolveWs(n)
b1 = copy(b)
work = zeros(n, n)
h0!(h0, b, c, g, work, ws_linsolve)
copy!(b1, b)
@time h0!(h0, b1, c, g, work, ws_linsolve)
@test h0 ≈ inv(b + c*g)

m = 10
hh = zeros(n, n*m)
work1 = zeros(n, n)
work2 = similar(work1)
hf = similar(work1)
hh!(hh, h0, a, hf, m, work1, work2)
hf = -h0*a
@test hh[:, 1:n] ≈ h0
for i = 1:n
    @test hh[:, (i-1)*n+1 : i*n] ≈ hf^(i-1)*h0
end

m = 10
hh = zeros(n, n*m)
work1 = zeros(n, n)
work2 = similar(work1)
hf = similar(work1)
hh!(hh, h0, a, hf, m, work1, work2)
hf = -h0*a
@test hh[:, 1:n] ≈ h0
for i = 1:m
    @test hh[:, (i-1)*n+1 : i*n] ≈ hf^(i-1)*h0
end

m = 3
hh = zeros(n, n*m)
work1 = zeros(n, n)
work2 = similar(work1)
hf = similar(work1)
hh!(hh, h0, a, hf, m, work1, work2)
hf = -h0*a
@test hh[:, 1:n] ≈ h0
for i = 1:m
    @test hh[:, (i-1)*n+1 : i*n] ≈ hf^(i-1)*h0
end

nn = 6
rin = rand(nn*n)
rout = similar(rin)
preconditioner!(rout, rin, g, hh, m, nn)
target = similar(rout)
target[1:6] = hh*rin[1:18]
target[7:12] = g*target[1:6]
target[7:12] += hh*rin[7:24]
target[13:18] = g*target[7:12]
target[13:18] += hh*rin[13:30]
target[19:24] = g*target[13:18]
target[19:24] += hh*rin[19:36]
target[25:30] = g*target[19:24]
target[25:30] += hh[:,1:12]*rin[25:36]
target[31:36] = g*target[25:30]
target[31:36] += hh[:,1:6]*rin[31:36]
@test rout ≈ target


m = context.models[1]
n = m.endogenous_nbr
endogenous = repeat(steadystate,3)
exogenous = zeros(3,m.exogenous_nbr)
work = context.work
get_jacobian!(work, endogenous, exogenous, steadystate, m, 2)
jacobian = work.jacobian

g1_1 = context.results.model_results[1].linearrationalexpectations.g1_1
g = zeros(n, n)
g[:,m.i_bkwrd_b] = g1_1
a = zeros(n, n)
b = zeros(n, n)
c = zeros(n, n)

get_abc!(a, b, c, jacobian, m)

h0 = zeros(n,n)
ws_linsolve = LinSolveWs(n)
bb = copy(b)
work = zeros(n, n)
h0!(h0, bb, c, g, work, ws_linsolve)
k = 6
hh = zeros(n, k*n)
work1 = zeros(n, n)
work2 = similar(work1)
hf = similar(work1)
hh!(hh, h0, c, hf, k, work1, work2)
rin = rand(k*n)
rout = similar(rin)
preconditioner!(rout, rin, g, hh, k, k)
 
jacobian1 = hcat(a, b, c)
A = makeA(jacobian1, 6)
AA = A[1:12, 1:12]
r = vcat(rin[1:6], zeros(6))
x = AA\r

@test x[1:6] ≈ h0*r[1:6] - hh[:,7:12]*r[7:12] - h0*c*hh[:,7:12]*a*x[7:12]

x = A\rin

@test A*x ≈ rin
@test -g*x[25:30] + x[31:36] ≈ h0*rin[31:36] + hh[:,7:12]*a*x[31:36]
@test x[1:6] ≈ hh*rin - h0*c*hh[:,31:36]*a*x[31:36]

periods = 300
A = makeA(jacobian1, periods)
rin = zeros(periods*n)
rin[4] = 0.01
rin[6] = 0.01
x = A\rin
hh = zeros(nn, periods*n)
hh!(hh, h0, c, hf, periods, work1, work2)
@test x[1:6] ≈ hh*rin

preconditioner_window = 6
n = 6
rout = similar(rin)
hh = zeros(n, preconditioner_window*n)
hh!(hh, h0, c, hf, preconditioner_window, work1, work2)
preconditioner!(rout, rin, g, hh, preconditioner_window, periods)
@test rout[1:6] ≈ hh*rin[1:preconditioner_window*n]
@test rout[7:12] ≈ hh*rin[(n+1):(preconditioner_window+1)*n] + g*rout[1:n]
#@test rout[295:300] ≈ h0*rin[295:300] + g*rout[289:294]

preconditioner_window = 60
P = LREprecond(periods, preconditioner_window, a, b, c, g)

ldiv!(rout, P, rin)
ldiv!(P, rin)
rout = P\rin

steadystate = context.results.model_results[1].trends.endogenous_steady_state
work = context.work
md = context.models[1]

k=1
endogenous = repeat(steadystate, k + 2)
exogenous = zeros(k + 2, md.exogenous_nbr)
exogenous[1, :] .= 0.0
residuals = zeros((k+2)*n)
y2 = zeros(k*n)
dynamic_variables = work.dynamic_variables
presiduals = zeros(md.n_bkwrd + md.n_current + md.n_fwrd + 2*md.n_both)
temp_vec = zeros(md.endogenous_nbr)
jacobian_time_vec!(y2, dynamic_variables, residuals, endogenous, exogenous,
               steadystate, presiduals, g, temp_vec, work, md, k)

LRE1 = LinearMap(k*n) do C, B
    copyto!(residuals, n + 1, B, 1, k*n) 
    jacobian_time_vec!(C, dynamic_variables, residuals, endogenous, exogenous,
                   steadystate, presiduals, g, temp_vec, work, md, k)
end

target = (jacobian[:,4:9] + jacobian[:, 10:12]*g[[1, 2, 6], :])*ones(6)
@test LRE1*ones(6) ≈ target


k=300
endogenous = repeat(steadystate, k + 2)
exogenous = zeros(k + 2, md.exogenous_nbr)
exogenous[2, :] .= 0.1
residuals = zeros(1812)
y2 = zeros(1800)
presiduals = zeros(md.n_bkwrd + md.n_current + md.n_fwrd + 2*md.n_both)
dynamic_variables = context.work.dynamic_variables
jacobian_time_vec!(y2, dynamic_variables, residuals, endogenous, exogenous,
               steadystate, presiduals, g, temp_vec, work, md, k)
P = LREprecond(k, 10, a, b, c, g)
LRE = LinearMap(k*n) do C, B
    copyto!(residuals, n + 1, B, 1, k*n) 
    jacobian_time_vec!(C, dynamic_variables, residuals, endogenous, exogenous,
                   steadystate, presiduals, g, temp_vec, work, md, k)
end

res = zeros(1800)
res[[4, 6]] .= 0.01
rout = zeros(1800)
x, h = gmres!(rout, LRE, res, log=true, verbose=true, Pr=P)


periods = 600
preconditioner_window = 50
res = zeros(periods*md.endogenous_nbr)
res[[4, 6]] .= 0.01
res = randn(periods*md.endogenous_nbr)
rout = zeros(periods*md.endogenous_nbr)
ws = GmresWs(periods, preconditioner_window, context, "CR")
ws.endogenous .= repeat(steadystate, periods + 2)
ws.exogenous .= zeros(periods + 2, md.exogenous_nbr)
gmres!(rout, ws.LREMap, res, log=false,
       verbose=true, Pr=ws.P)


#gmres_solver!(rout, res, periods, preconditioner_window, md, work, ws, verbose=true)

