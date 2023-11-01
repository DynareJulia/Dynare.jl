using Dynare
using Test

include("../src/global/blocks.jl")

context = @dynare "models/example3/example3" "notmpterms" "output=second"

matching = get_maximum_cardinality_matching(context)
U = get_incidence_bitmatrix_current_forward(context)

rb = get_recursive_blocks(context, matching, U)

f = Dynare.DFunctions.SparseDynamicResid!

forward_expressions = Expr[]
preamble_expressions = Expr[]
system_expressions = Expr[]

Blocks.analyze_SparseDynamicResid!(forward_expressions, preamble_expressions, system_expressions, rb, matching, context, context.models[1].endogenous_nbr)

Blocks.make_block_functions(context)

steadystate = context.results.model_results[1].trends.endogenous_steady_state
y = repeat(steadystate, 3)
x = context.results.model_results[1].trends.exogenous_steady_state
parameters = context.work.params

x .= [0.5, -0.5]
y[[4, 6]] .= [-0.8, 0.5]

Blocks.preamble_block(y, x, parameters, steadystate)

@test y[11] ≈ parameters[2] * y[4] + parameters[7] * y[6] + x[1]
@test y[13] ≈ y[4] * parameters[7] + parameters[2] * y[6] + x[2]

(equation_xref_list, variable_xref_list) = Blocks.xref_lists(context)

@test equation_xref_list[1] == [8, 9, 12]
@test equation_xref_list[2] == [10, 21]
@test equation_xref_list[3] == [ 3, 8, 11, 12]
@test variable_xref_list[1] == Int[]
@test variable_xref_list[3] == [3, 4, 7]
@test variable_xref_list[9] == [1, 4, 7]
@test variable_xref_list[21] == [2]

@show Blocks.is_block_linear(rb[1], context)
@show Blocks.is_block_linear(rb[2], context)
@show Blocks.is_block_linear(rb[3], context)
