using ExtendedDates

function find_letter_in_period(period::AbstractString)
    for c in period
        if 'A' <= c <= 'z'
            return uppercase(c)
        end
    end
    return nothing
end

function identify_period_type(period::AbstractString)
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

function periodparse(period::AbstractString)::ExtendedDates.DatePeriod
    period_type = identify_period_type(period)
    if period_type == YearSE
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("yY"))
    elseif period_type == SemesterSE
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("ySs"))
    elseif period_type == QuarterSE
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("yQq"))
    elseif period_type == MonthSE
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("yMm"))
    elseif period_type == WeekSE
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("yWw"))
    elseif period_type == DaySE
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("y-m-d"))
    elseif period_type == Undated
        return parse(period_type, period, ExtendedDates.SimpleDateFormat("x"))
    else
        throw(ErrorException)
    end
end

# To be moved to TimeDataFrames.jl
function MyTimeDataFrame(filename)
    df = DataFrame(CSV.File(filename))
    local first_period::Any
    periods = []
    continuous = true
    for name in names(df)
        if uppercase(name) in ["DATE", "DATES", "PERIOD", "PERIODS"]
            p = df[!, name]
            foreach(x -> push!(periods, periodparse(x)), p)
            if !is_continuous(periods)
                continuous = false
            end
        else
            periods = range(
                UndatedDate(1),
                length = DataFrames.nrow(df),
                step = ExtendedDates.Undated(1),
            )
        end
        break
    end
    return TimeDataFrame(df, periods, continuous)
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

function get_data(
    filename::String,
    variables::Vector{String};
    start::Int64 = 1,
    last::Int64 = 0,
)
    tdf = MyTimeDataFrame(filename)
    ny = length(variables)

    if last == 0
        last = size(tdf, 1) - start + 1
    end
    nobs = last - start + 1
    Y = Matrix{Union{Missing,Float64}}(undef, ny, nobs)
    for (i, v) in enumerate(variables)
        Y[i, :] .= tdf[!, Symbol(v)][start:last]
    end
    return Y
end
