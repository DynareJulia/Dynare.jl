using Dynare
using Test

context = @dynare "models/example3ss/example3ss.mod"

field = Dict("statementName" => "initval",
             "vals" => [Dict("name" => "y", "value" => "1"),
                        Dict("name" => "c", "value" => "1"),
                        Dict("name" => "k", "value" => "1"),
                        Dict("name" => "h", "value" => "1"),
                        Dict("name" => "a", "value" => "1"),
                        Dict("name" => "b", "value" => "1")])

Dynare.initval!(context, field)
y = context.work.initval_endogenous
@test y[1:6] == ones(6)
# auxiliary variable
params = context.work.params
@test y[7] ≈ y[2]*exp(y[6])/(exp(y[6])*0.5*(y[2]+y[2]))*(y[3]*(1-params[4])+y[1]*params[3]*exp(y[6]));

