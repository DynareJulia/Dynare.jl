import AxisArrays
using AxisArrayTables
using ExtendedDates

function data!(datafile::AbstractString;
    context::Context = context,
    variables::Vector{<:Union{String,Symbol}} = [],
    start::PeriodsSinceEpoch = Undated(typemin(Int)),
    last::PeriodsSinceEpoch = Undated(typemin(Int)),
    nobs::Integer = 0,
)
    aat = MyAxisArrayTable(datafile)
    Ta = typeof(row_labels(aat)[1])
    if Ta <: Dates.UTInstant
        T = Ta.parameters[1]
    elseif Ta == Int
        T = Int
    else
        error("Unrecognized type")
    end 
    Ts = typeof(start)
    Tl = typeof(last)
    ny = length(variables)
    if typeof(start) != Int || start > typemin(Int)
        # start option is used
        @assert Ts == Ta "error in data!(): start must have the same frequency as the datafile"
        startperiod = start
        if typeof(last) != Int || last > typemin(Int)
            # start option and last option are used
            @assert Tl == Ta "error in data!(): last must have the same frequency as the datafile"
            lastperiod = last
        elseif nobs > 0
            lastperiod = startperiod + T(nobs) - T(1)
        else
            lastperiod = row_labels(aat)[end]
        end
    elseif typeof(last) != Int || last > typemin(Int)
        # start option isn't used but last option is used
        @assert Tl == Ta "error in data!(): last must have the same frequency as the datafile"
        lastperiod = last
        if nobs > 0
            startperiod = last - T(nobs) + T(1)
        else
            startperiod = row_labels(aat)[1]
        end
    else
        # neither start option nor last option are used
        startperiod = row_labels(aat)[1]
        if nobs > 0
            lastperiod = start + T(nobs) - T(1)
        else
            lastperiod = row_labels(aat)[end]
        end
    end
    if length(variables) == 0
        context.work.data = copy(aat[startperiod:lastperiod, :])
    else
        context.work.data = copy(aat[startperiod:lastperiod, Symbol.(variables)])
    end 
    return context.work.data
end

function get_data!(context::Context,
                   datafile::String,
                   data::AxisArrayTable,
                   variables::Vector{<:Union{String, Symbol}},
                   first_obs::PeriodsSinceEpoch,
                   last_obs::PeriodsSinceEpoch,
                   nobs::Int
                   )
    @assert isempty(datafile) || isempty(data) "datafile and data can't be used at the same time"
    
    if !isempty(datafile)
        data!(datafile, 
              context = context,
              variables = variables,
              start = first_obs,
              last = last_obs,
              nobs = nobs)
    elseif !isempty(data)
        first_obs == Undated(typemin(Int)) && (first_obs = row_labels(data)[1])
        last_obs == Undated(typemin(Int)) && (last_obs = row_labels(data)[end])
        context.work.data = copy(data[first_obs:last_obs, Symbol.(variables)])
    else
        error("needs datafile or data argument")
    end
    return context.work.data
end 

function get_detrended_data(context::Context,
    datafile::String,
    data::AxisArrayTable,
    variables::Vector{<:Union{String, Symbol}},
    first_obs::PeriodsSinceEpoch,
    last_obs::PeriodsSinceEpoch,
    nobs::Int
    )

    endogenous_names = get_endogenous(context.symboltable)
    trends = context.results.model_results[1].trends
    steady_state = Vector(AxisArrays.AxisArray(trends.endogenous_steady_state, endogenous_names)[variables])
    linear_trend = Vector(AxisArrays.AxisArray(trends.endogenous_linear_trend, endogenous_names)[variables])
    quadratic_trend = Vector(AxisArrays.AxisArray(trends.endogenous_quadratic_trend, endogenous_names)[variables])

    get_data!(context,
    datafile,
    data,
    variables,
    first_obs,
    last_obs,
    nobs
    )
    aat = copy(context.work.data)
    if context.modfileinfo.has_trends
        if !isempty(quadratic_trend)
            remove_quadratic_trend!(aat, 
            adjoint(steady_state), 
            adjoint(linear_trend), 
            adjoint(quadratic_trend)
            )
        else
            remove_linear_trend!(aat, 
            adjoint(steady_state), 
            adjoint(linear_trend), 
            )
        end
    else
        aat .-= adjoint(steady_state)
    end
    return aat
end 

function find_letter_in_period(period::AbstractString)
    for c in period
        if 'A' <= c <= 'z'
            return uppercase(c)
        end
    end
    return nothing
end

function identify_period_type(period::Union{AbstractString, Number, Date})
    if typeof(period) <: Date
        return Date
    elseif typeof(period) <: Number
        if isinteger(period)
            if period == 1
                return Undated
            else
                return YearSE
            end
        else
            throw(ErrorException)
        end
    else
        c = find_letter_in_period(period)
        let period_type
            if isnothing(c)
                if (cdash = count("-", period)) == 2
                    period_type = DaySE
                elseif cdash == 1
                    period_type = MonthSE
                elseif tryparse(Int64, period) !== nothing
                    period_type = Undated
                else
                    throw(ErrorException)
                end
            else
                if c == 'Y'
                    period_type = YearSE
                elseif c == 'A'
                    period_type = YearSE
                elseif c == 'S'
                    period_type = ExtendedDates.SemesterSE
                elseif c == 'H'
                    period_type = ExtendedDates.SemesterSE
                elseif c == 'Q'
                    period_type = QuarterSE
                elseif c == 'M'
                    period_type = MonthSE
                elseif c == 'W'
                    period_type = WeekSE
                elseif c == 'D'  # day of the year: 1980D364
                    period_type = DaySE
                else
                    throw(ErrorException)
                end
            end
            return period_type
        end
    end
end

function periodparse(period::Union{AbstractString, Number, Date})::ExtendedDates.PeriodsSinceEpoch
    period_type = identify_period_type(period)
    if period_type == YearSE
        if typeof(period) <: Number
            return(YearSE(Int(period)))
        else
            return parse(period_type, period)
        end
    elseif period_type == SemesterSE
        return parse(period_type, period)
    elseif period_type == QuarterSE
        return parse(period_type, period)
    elseif period_type == MonthSE
        return parse(period_type, period)
    elseif period_type == WeekSE
        return parse(period_type, period)
    elseif period_type == DaySE
        return parse(period_type, period)
    elseif period_type == Date
        return DaySE(period)
    elseif period_type == Undated
        return Undated(period)
    else
        throw(ErrorException)
    end
end

function MyAxisArrayTable(filename)
    table = CSV.File(filename)
    cols = AxisArrayTables.Tables.columnnames(table)
    data = (AxisArrayTables.Tables.matrix((;ntuple(i -> Symbol(i) => Tables.getcolumn(table, i), length(cols))...)))
    for (icol, name) in enumerate(cols)
        if uppercase(String(name)) in ["DATE", "DATES", "PERIOD", "PERIODS", "TIME", "COLUMN1"]
            rows = []
            foreach(x -> push!(rows, periodparse(x)), data[:, icol])
            k = union(1:icol-1, icol+1:size(data,2))
            aat = AxisArrayTable(Matrix{Union{Float64, Missing}}(data[:, k]), rows, cols[k])
            return aat
        end
    end
    rows = Undated(1):Undated(size(data,1))
    aat = AxisArrayTable(data, rows, cols)
    return aat
end

function is_continuous(periods::Vector{ExtendedDates.DatePeriod})
    i = periods[1]
    for j = 2:length(periods)
        if periods[j] != i + 1
            return false
        else
            i += 1
        end
    end
end
