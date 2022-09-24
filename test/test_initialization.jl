using Dynare
using ExtendedDates
using Test

period_arg(a, b, c) = (isnothing(a)) ? nothing : "$a$b$c"

function test_period(a, b, c, d)
    if isnothing(a)
        @test isnothing(c)
        return
    end
    if isempty(d)
        @test a == b(c)
    else
        @test a == b(c, d)
    end
end

context = @dynare "models/example3/example3.mod"
context = @dynare "models/example3ss/example3ss.mod"

field = Dict(
    "statementName" => "initval",
    "vals" => [
        Dict("name" => "y", "value" => "1"),
        Dict("name" => "c", "value" => "1"),
        Dict("name" => "k", "value" => "1"),
        Dict("name" => "h", "value" => "1"),
        Dict("name" => "a", "value" => "1"),
        Dict("name" => "b", "value" => "1"),
    ],
)

Dynare.initval!(context, field)
y = context.work.initval_endogenous
@test y[1:6] == ones(6)
# auxiliary variable
params = context.work.params
@test y[7] == 1

let first_obs, last_obs, first_simulation_period, last_simulation_period, nobs, period_type
    options = Dict(
        "first_obs" => "1979Y",
        "last_obs" => "2000Y",
        "first_simulation_period" => "1980Y",
        "last_simulation_period" => "1999Y",
        "nobs" => "22",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == YearDate(1979)
    @test last_obs == YearDate(2000)
    @test first_simulation_period == YearDate(1980)
    @test last_simulation_period == YearDate(1999)
    @test period_type == YearDate

    options = Dict(
        "first_obs" => "1979Q3",
        "last_obs" => "2000Q3",
        "first_simulation_period" => "1980Q3",
        "last_simulation_period" => "1999Q3",
        "nobs" => "85",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == QuarterDate(1979, 3)
    @test last_obs == QuarterDate(2000, 3)
    @test first_simulation_period == QuarterDate(1980, 3)
    @test last_simulation_period == QuarterDate(1999, 3)
    @test period_type == QuarterDate

    options = Dict(
        "first_obs" => "1979",
        "last_obs" => "2000",
        "first_simulation_period" => "1980",
        "last_simulation_period" => "1999",
        "nobs" => "22",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == UndatedDate(1979)
    @test last_obs == UndatedDate(2000)
    @test first_simulation_period == UndatedDate(1980)
    @test last_simulation_period == UndatedDate(1999)
    @test period_type == UndatedDate

    options = Dict(
        "first_obs" => "1979M3",
        "last_obs" => "2000M3",
        "first_simulation_period" => "1980M3",
        "last_simulation_period" => "1999M3",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == MonthDate(1979, 3)
    @test last_obs == MonthDate(2000, 3)
    @test first_simulation_period == MonthDate(1980, 3)
    @test last_simulation_period == MonthDate(1999, 3)
    @test period_type == MonthDate
    @test isnothing(nobs)

    options = Dict(
        "first_obs" => "1979M3",
        "last_obs" => "2000M3",
        "first_simulation_period" => "1980M3",
        "nobs" => "253",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == MonthDate(1979, 3)
    @test last_obs == MonthDate(2000, 3)
    @test first_simulation_period == MonthDate(1980, 3)
    @test isnothing(last_simulation_period)
    @test period_type == MonthDate

    options = Dict(
        "first_obs" => "1979M3",
        "last_obs" => "2000M3",
        "last_simulation_period" => "1999M3",
        "nobs" => "253",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == MonthDate(1979, 3)
    @test last_obs == MonthDate(2000, 3)
    @test isnothing(first_simulation_period)
    @test last_simulation_period == MonthDate(1999, 3)
    @test period_type == MonthDate

    options = Dict(
        "first_obs" => "1979M3",
        "first_simulation_period" => "1980M3",
        "last_simulation_period" => "1999M3",
        "nobs" => "253",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == MonthDate(1979, 3)
    @test isnothing(last_obs)
    @test first_simulation_period == MonthDate(1980, 3)
    @test last_simulation_period == MonthDate(1999, 3)
    @test period_type == MonthDate

    options = Dict(
        "last_obs" => "2000M3",
        "first_simulation_period" => "1980M3",
        "last_simulation_period" => "1999M3",
        "nobs" => "253",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test isnothing(first_obs)
    @test last_obs == MonthDate(2000, 3)
    @test first_simulation_period == MonthDate(1980, 3)
    @test last_simulation_period == MonthDate(1999, 3)
    @test period_type == MonthDate

    options =
        Dict("last_obs" => "2000M3", "last_simulation_period" => "1999M3", "nobs" => "253")
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test isnothing(first_obs)
    @test last_obs == MonthDate(2000, 3)
    @test isnothing(first_simulation_period)
    @test last_simulation_period == MonthDate(1999, 3)
    @test period_type == MonthDate

    options = Dict(
        "first_obs" => "1979M3",
        "first_simulation_period" => "1980M3",
        "nobs" => "253",
    )
    (
        first_obs,
        last_obs,
        first_simulation_period,
        last_simulation_period,
        nobs,
        period_type,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
    @test first_obs == MonthDate(1979, 3)
    @test isnothing(last_obs)
    @test first_simulation_period == MonthDate(1980, 3)
    @test isnothing(last_simulation_period)
    @test period_type == MonthDate

    options = Dict(
        "first_obs" => "1979Y",
        "last_obs" => "2000",
        "first_simulation_period" => "1980Y",
    )
    @test_throws AssertionError (
        first_obs1,
        last_obs1,
        first_simulation_period1,
        last_simulation_period1,
        period_type1,
    ) = Dynare.check_periods_options(options, required_lags = 1, required_leads = 1)
end

field = Dict("statementName" => "histval_file", "options" => Dict("datafile" => "data.csv"))

Dynare.histval_file!(context, field)

field = Dict("statementName" => "initval_file", "options" => Dict("datafile" => "data.csv"))
Dynare.initval_file!(context, field)
