using Dynare
using PATHSolver
using SparseArrays

# Using PATHSolver
function solve_path!(F, J, lb, ub, initial_values; kwargs...)

    function function_callback(n::Cint, z::Vector{Cdouble}, F_val::Vector{Cdouble})
        @assert n == length(z) == length(F_val)
        F_val .= F(z)
        @show F_val
        return Cint(0)
    end

    function jacobian_callback(        
        n::Cint,
        nnz::Cint,
        z::Vector{Cdouble},
        col_start::Vector{Cint},
        col_len::Vector{Cint},
        row::Vector{Cint},
        data::Vector{Cdouble}
    )
        @assert n == length(z) == length(col_start) == length(col_len)
        @assert nnz == length(row) == length(data)    

        Jac = J(z)  # SparseMatrixCSC

        # Pouring Jac::SparseMatrixCSC to the format used in PATH
        for i in 1:n
            col_start[i] = Jac.colptr[i]
            col_len[i] = Jac.colptr[i + 1] - Jac.colptr[i]
        end

        rv = rowvals(Jac)
        nz = nonzeros(Jac)
        num_nonzeros = SparseArrays.nnz(Jac)
        for i in 1:num_nonzeros
            row[i] = rv[i]
            data[i] = nz[i]
        end

        return Cint(0)
    end

    # var_name, F_name in RawIndex
    var_name = Vector{String}(undef, 0)
    F_name = Vector{String}(undef, 0)

    # overestimating number of nonzeros in Jacobian
    nnz = max( SparseArrays.nnz(J(lb)), SparseArrays.nnz(J(ub)) )
    nnz = max( nnz, SparseArrays.nnz(J(initial_values)) )
    for i in 1:2
        z_rand = max.(lb, min.(ub, rand(Float64, n)))
        nnz_rand = SparseArrays.nnz(J(z_rand))
        nnz = max( nnz, nnz_rand )
    end
    nnz = min( 2 * nnz, n^2 )


    F_val = zeros(3)
    @show function_callback(Int32(3), initial_values, F_val)
    @show F_val
    # Solve the MCP using PATHSolver
    status, z, info = PATHSolver.solve_mcp(
        function_callback, 
        jacobian_callback, 
        lb, 
        ub, 
        initial_values; 
        nnz = nnz, 
        variable_names = var_name,
        constraint_names = F_name,
        kwargs...
    )

    # This function has changed the content of m already.
    return Int(status), z, info
end

n = 3
A = spdiagm(rand(3))
F(x) = A*x
J(x) = A
lb = repeat([-Inf], n)
lb[3] = 3
ub = repeat([Inf], n)
initial_values = rand(n)

status, z, info = solve_path!(F, J, lb, ub, initial_values)

modeljson = Dynare.parseJSON("/home/michel/projects/julia-packages/dynare.jl/test/models/irreversible/irbc1a")

function ramsey_contraints(context, field)
    for c in field["ramsey_model_constraints"]
        push!(context.work.mcps, c["constraint"])
    end
end
    
function get_tags!(mcps, modeljson)
    for (i, eq) in enumerate(modeljson["model"])
        tag = get(eq["tags"], "mcp", "")
        if !isempty(tag)
            push!(mcps, (i, tag))
        end
    end
end

function mcp!(lb, ub, permutations, mcps, context)
    for m in mcps
        (m1, m2) = m 
        (var, op, expr) = split(m2)
        
        if op[1] == '<'
            ub[m1] = Dynare.dynare_parse_eval(String(expr), context)
        elseif op[1] == '>'
            lb[m1] = Dynare.dynare_parse_eval(String(expr), context)
        else
            error("MCP operator must be <, <=, > or >=")
        end
        iva = context.symboltable[var].orderintype
        if iva != m1
            push!(permutations, (iva, m1))
        end
    end
end

mcps = []
get_tags!(mcps, modeljson)
lb = repeat([-Inf], 3)
ub = repeat([Inf], 3)
permutations = []
mcp!(lb, ub, permutations, mcps, context)

function reorder!(x, permutations, offset)
    for p in permutations
        p1 = p[1] + offset
        p2 = p[2] + offset
        x[p2], x[p1] = x[p1], x[p2]
    end
end

x = collect(1:10)
reorder!(x, [(1, 3), (2, 6)], 0)
@show x
reorder!(x, [(1, 3), (2, 6)], 0)
@show x
x = collect(1:20)
reorder!(x, [(1, 3), (2, 6)], 10)
@show x
