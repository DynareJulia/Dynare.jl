using Dynare

context = @dynare "fs2000";

D = Dynare.CSV.File("data.csv")
h = Dynare.CSV.Tables.columnnames(D)
m = Dynare.CSV.Tables.matrix(D)
data = Dynare.AxisArrayTable(m, Dynare.Undated(1):Dynare.Undated(size(m, 1)),h )
gy_obs = diff(log.(data.GDPQ))
gp_obs = diff(log.(data.GDPD))
data1 = Dynare.AxisArrayTable([gy_obs gp_obs][2:end,:], Dynare.Undated(1):Dynare.Undated(size(m, 1)-1), [:gy_obs, :gp_obs])
Dynare.estimation!(context, Dict{String,Any}("options" => Dict("data" => data1))         )