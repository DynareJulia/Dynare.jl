using FastLapackInterface
using FastLapackInterface.QrAlgo
using LinearAlgebra

struct BalancedGrowthWs
    A1::Matrix{Float64}
    A2::Matrix{Float64}
    dd::Vector{Float64}
    ws::QrpWs
    ws1::QrpWs
    function BalancedGrowthWs(n)
        A1 = Matrix{Float64}(undef, n, n)
        A2 = Matrix{Float64}(undef, n, n)
        dd = Vector{Float64}(undef, 2*n)
        qrws = QrpWs(n)
        qrws1 = QrpWs(2*n)
    end
end

function balanced_growth_path(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64}, C::AbstractMatrix{Float64}, d::AbstractVector{Float64}, ws::BalancedGrowthWs)
    vJ = view(jacobian, :, ws.backward_nbr .+ ws.current_nbr .+ (1:ws.forward_nbr)) 
    ws.A1 .+= vj
    vJ = view(jacobian, :, ws.backward_nbr .+ ws.current_dynamic_indices) 
    ws.A1 .+= vj
    ws.B .+= vj
    vJ = view(jacobian, :, 1:ws.backward_nbr)
    ws.A1 .+= vj
    balanced_growth_compute(ws.A1, ws.B, d, ws)
end

function balanced_growth_path(A::AbstractMatrix{Float64}, B::AbstractMatrix{Float64}, d::AbstractVector{Float64}, ws::BalancedGrowthWs)
    ws.A1 .= A .+ B
    balanced_growth_compute(ws.A1, B, d, ws)
end

function balanced_growth_compute(A1::AbstractMatrix{Float64}, A2::AbstractMatrix{Float64}, d::AbstractVector{Float64}, ws::BalancedGrowthWs)
    n = length(d)
    println("m $n")
    ws = QrpWs(A1)
    geqp3!(A1, ws)
    
    M1 =  zeros(2*n, 2*n)
    view(M1, 1:n, 1:n) .= A1
    view(M1, n .+ (1:n), 1:n) .= A1
    view(M1, n .+ (1:n), n .+ (1:n)) .= A2
    ws1 = QrpWs(M1)
    geqp3!(M1, ws1)
    println("triu(M1)")
    display(triu(M1))
    println(ws1.jpvt)
    d1 = zeros(2*n)
    view(d1, 1:n) .= d
    view(d1, n .+ (1:n)) .= d
    dd1 = Matrix{Float64}(undef, 2*n, 1)
    dd1 .= d1
    ormqr_core!('L', transpose(M1), dd1, ws1)
    for i = 1:n
        ws.dd[i] = d1[ws1.pivot[i]]
    end
    for i = n+1 : 2*n
        ws.dd[i] = d1[ws1.pivot[i] - n]
    end
    println(dd1)
    
end

function get_abc!(ws::LinearRationalExpectationsWs, jacobian::AbstractMatrix{Float64})
    i_rows = (ws.static_nbr+1):ws.endogenous_nbr
    fill!(ws.a, 0.0)
    fill!(ws.b, 0.0)
    fill!(ws.c, 0.0)
    ws.a[:, ws.forward_indices_d] .= view(jacobian, i_rows, ws.backward_nbr .+ ws.current_nbr .+ (1:ws.forward_nbr))
    ws.b[:, ws.current_dynamic_indices_d] .= view(jacobian, i_rows, ws.backward_nbr .+ ws.current_dynamic_indices)
    ws.c[:, ws.backward_indices_d] .= view(jacobian, i_rows, 1:ws.backward_nbr)
end

function balance_growth_path_general!(m, work)
    get_dynamic_jacobian!(work,
                          zeros(m.endogenous_nbr),
                          zeros(m.exogenous_nbr),
                          m,
                          2)

end

function balance_growth_path_purely_backward!()
end

function balance_growth_path_purely_forward!()
end

A = [ 1.0 -1 -1 0; 0 0 1 0; 0 1 0 0; 0 0 2 -1]
B = [0 0 0 0; 0 0 -0.9 0; 0 -1 0 0; 0 0 0 0]
d = [0, 0, 2.0, -3.0]
#=
A = [ 1.0 -1 -1 ; 0 0 1 ; 0 1 0 ]
B = [0 0 0 ; 0 0 -0.9 ; 0 -1 0]
d = [0, 0, 2.0]
=#
n = length(d)
AB = A+B
C0 = AB
q0 = qr(C0, Val(true))
println("q0.R")
display(q0.R)
Zn = zeros(n, n)
C1 = vcat(hcat(AB, Zn), hcat(AB, A))
display(C1)
d1 = vcat(d, d)
q1 = qr(C1, Val(true))
qd1 = q1.Q'*d1
k = findall(abs.(diag(q1.R)) .> 1e-12)
Q = q1.Q[k, k]
R = q1.R[k, k]
P = q1.P[:, k]
@show P*(R\qd1[k])

balanced_growth_path(A + B, A, d)
