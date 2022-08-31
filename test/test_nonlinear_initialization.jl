using Dynare
using Test

tdf = Dynare.TimeDataFrame(
    Dynare.DataFrame(rand(10, 3), ["x", "y", "e"]),
    Dynare.UndatedDate(1),
)
Dynare.CSV.write("data_nl.csv", getfield(tdf, :data))

context = @dynare "models/nonlinear_initialization/nl_init.mod";

y = tdf[!, :y]
a1 = context.work.params[1]
z1 =
    y[2] * y[2] + a1 * y[2]^2 - exp(y[2]) + exp(log(y[2])) + sin(y[2]) - cos(y[2]) +
    y[2] / y[1]
z2 =
    y[3] * y[3] + a1 * y[3]^2 - exp(y[3]) + exp(log(y[3])) + sin(y[3]) - cos(y[3]) +
    y[3] / y[2]

aux = context.work.initval_endogenous[:, 3]

@test aux[3] ≈ z2 - z1
