using CSV
using DataFrames
using Dynare
using ExtendedDates
using TimeDataFrames

tdf = TimeDataFrame(DataFrame(rand(10,3), ["x", "y", "e"]), UndatedDate(1))
CSV.write("data_nl.csv", getfield(tdf, :data))

context = @dynare "models/nonlinear_initialization/nl_init.mod";

aux = context.models[1].initval_endogenous

y = tdf[!, :y]

y*y + a1*y^2 - exp(y) + exp(ln(y)) + sin(y) - cos(y) + y/y(-1)
