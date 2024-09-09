n_static(i::LinearRationalExpectations.Indices)     = length(i.static)
n_forward(i::LinearRationalExpectations.Indices)    = length(i.forward)
n_backward(i::LinearRationalExpectations.Indices)   = length(i.backward)
n_both(i::LinearRationalExpectations.Indices)       = length(i.both)
n_current(i::LinearRationalExpectations.Indices)    = length(i.current)
n_dynamic(i::LinearRationalExpectations.Indices)    = length(i.dynamic)
n_endogenous(i::LinearRationalExpectations.Indices) = i.n_endogenous
n_exogenous(i::LinearRationalExpectations.Indices) = length(i.exogenous)

function get_system_jacobian(backward_block, forward_block, exogenous, parameters, steadystate)
    tempterms = []
    dyn_endogenous = repeat(steadystate, 3)
    backward_block.update_jacobian!(tempterms, backward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
    forward_block.update_jacobian!(tempterms, forward_block.jacobian.nzval, dyn_endogenous, exogenous, parameters, steadystate)
    return vcat(backward_block.jacobian, forward_block.jacobian)
end

function block_first_order_approximation(context, sgws)
    backward_block = sgws.backward_block
    forward_block  = sgws.forward_block 
    preamble_block  = sgws.preamble_block
    policy_jacobian = sgws.policy_jacobian
    exogenous      = sgws.sgmodel.exogenous     
    parameters     = sgws.sgmodel.parameters   
    steadystate    = sgws.sgmodel.steadystate
    system_variables = sgws.system_variables
    endogenous_nbr = sgws.sgmodel.endogenous_nbr
    exogenous_nbr = sgws.sgmodel.exogenous_nbr
    
    J = get_system_jacobian(backward_block, forward_block, exogenous, parameters, steadystate)
    @views begin
        state_variables = preamble_block.variables
        ns = length(system_variables)
        nx = length(state_variables)
        A = J[:, system_variables .+ 2*endogenous_nbr]
        B = J[:, system_variables .+ endogenous_nbr]
        C = J[:, system_variables]
        D = J[:, state_variables .+ 2*endogenous_nbr]
        E = J[:, state_variables .+ endogenous_nbr]
        F = J[:, state_variables]
        G = Matrix(preamble_block.jacobian[:, preamble_block.variables])
        H = Matrix(preamble_block.jacobian[:, 3*endogenous_nbr .+ (1:nx)])
        M = zeros(ns, ns)
        N = zeros(ns, nx)
        P = zeros(ns, nx)
    end

    backward_nbr = length(state_variables)
    solve_extended_linear_model!(context, A, B, C, D, E, F, G, H, M, N, P, endogenous_nbr, exogenous_nbr, backward_nbr, state_variables, system_variables)
    return M, N, P
end

"""
solves
```
    Ay_{t+1} + By_t + Cy_{t-1} + Dx_{t+1} + Ex_t + Fx_{t-1} = 0
```
with
```
    x_t = G x_{t_1} + Hu_t
```
The solution has the form
```
    y_t = My_{t-1} + Nx_t + Px_{t-1}
```
"""
function solve_extended_linear_model!(context, A, B, C, D, E, F, G, H, M, N, P, endogenous_nbr, exogenous_nbr, backward_nbr, state_variables, system_variables)
    options = LinearRationalExpectationsOptions()
    results = LinearRationalExpectationsResults(length(system_variables), length(state_variables), backward_nbr)
    model = context.models[1]
    sv = system_variables
    lli = model.lead_lag_incidence[:, sv]
    forward = findall(lli[3,:] .> 0)
    current = findall(lli[2,:] .> 0)
    backward = findall(lli[1,:] .> 0)
    static = findall((lli[2, :] .> 0) .& (lli[1, :] .== 0) .& (lli[3, :] .== 0) )
    ids = LRE.Indices(length(state_variables), forward, current, backward, static)
    ws = LRE.LinearRationalExpectationsWs("GS", ids)
    jacobian = hcat(C, B, A, E)
    ncol = model.n_bkwrd + model.n_current + model.n_fwrd + 2 * model.n_both + model.exogenous_nbr
    compressed_jacobian = zeros(size(jacobian, 1), ncol)
    set_compressed_jacobian!(compressed_jacobian, jacobian, length(system_variables), length(state_variables),  lli)
    solve_M!(context, M, results, compressed_jacobian, options, ws, system_variables)
    solve_P!(context, P, A, B, F, M)
    solve_N!(context, N, A, B, D, E, G, M, P)
    return nothing
end

"""
Solves  AMM + BM + C = 0
"""
function solve_M!(context, M, results, jacobian, options, ws, system_variables)
    model = context.models[1]
    ids = ws.ids
    n_back    = n_backward(ids)
    n_cur     = n_current(ids)
    forward_r = 1:n_forward(ids)
    current_r = 1:n_cur
    back      = ids.backward
    
    if n_static(ids) > 0
        LRE.remove_static!(jacobian, ws)
    end
    
    n_back    = n_backward(ids)
    n_cur     = n_current(ids)
    forward_r = 1:n_forward(ids)
    current_r = 1:n_cur
    back      = ids.backward
    
    @views @inbounds begin
        LRE.copy_jacobian!(ws.solver_ws, jacobian)
        ws.jacobian_current[:, current_r] .= jacobian[:, n_back .+ current_r]
        results.g1, results.gs1 = LRE.solve_g1!(results, ws.solver_ws, options)
        if n_static(ids) > 0
            results.g1, jacobian = LRE.add_static!(results, jacobian, ws)
        end
        
        results.gns1 .= results.g1_1[ids.non_backward, :]
    end
    M[:, back] .= results.g1_1
end

"""
Solves AMP + BP + F = 0
"""
function solve_P!(context, P, A, B, F, M)
    B_ = copy(B)
    mul!(B_, A, M, 1, 1)
    F_ = Matrix(F)
    P = -B_\F_
    return P
end

"""
Solves AMN + ANG + AP + BN + DG + E = 0
"""
function solve_N!(context, N, A, B, D, E, G, M, P)
    nr, nc = size(N)
    N .= -reshape((kron(I(nc), A*M+B) + kron(transpose(G), A))\vec(A*P + D*G + E), nr, nc)
    return N
end
                 
function set_compressed_jacobian!(compressed_jacobian, jacobian, n, nx, lli)
    k = 1
    m = 1
    for i in axes(lli, 1)
        for j = 1:n
            if lli[i, j] > 0
                copyto!(compressed_jacobian, k, jacobian, m, n)
                k += n
            end
            m += n
        end
    end
    # jacobian has 3*n + nx columns
    copyto!(compressed_jacobian, k, jacobian, 3*n*n + 1, n*nx)
    return nothing
end


