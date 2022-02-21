function identify_period_frequency(period)::ExtendedDates.Frequency
    period = uppercase(period)
    if 'Y' in period
        frequency = Year
    elseif 'A' in period
        frequency = Year
    elseif 'S' in period
        frequency = ExtendedDates.Semester
    elseif 'H' in period
        frequency = ExtendedDates.Semester
    elseif 'Q' in period
        frequency = Quarter
    elseif 'M' in period
        frequency = Month
    elseif 'W' in period
        frequency = Week
    elseif 'D' in period
        frequency = Day
    elseif (cdash = count("-", period)) == 2
        frequency = Day
    elseif cdash == 1
        frequency = Month
    else
        throw(ErrorException)
    end
end

# To be moved to TimeDataFrames.jl
function MyTimeDataFrame(filename)
    df = DataFrame(CSV.File(filename))
    local first_period::Any
    continuous = true
    for name in names(df)
        if uppercase(name) in ["DATE", "DATES", "PERIOD", "PERIODS"]
            p = df[!, name]
            periods = []
            foreach(x -> push!(periods, parse(x)), p)
            if !is_continuous(periods)
                continuous = false
            end
        else
            periods = range(UndatedDate(1), length=nrow(df), step= ExtendedDates.Undated(1)) 
        end
        break
    end
    return TimeDataFrame(df, periods, continuous)
end

function is_continuous(periods::Vector{ExtendedDates.SimpleDate})
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
        last = size(df, 1) - start + 1
    end
    nobs = last - start + 1
    Y = Matrix{Union{Missing,Float64}}(undef, ny, nobs)
    for (i, v) in enumerate(variables)
        Y[i, :] .= tdf[!, Symbol(v)][start:last]
    end
    return Y
end
