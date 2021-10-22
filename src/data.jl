function identify_period_frequency(period)::Periods.Frequency
    period = uppercase(period)
    if 'Y' in period
        frequency = Year
    elseif 'A' in period
        frequency = Year
    elseif 'S' in period
        frequency = Period.Semester
    elseif 'H' in period
        frequency = Semester
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

function get_data(filename::String, variables::Vector{String};
                  start::Int64 = 1, last::Int64 = 0)
    df = DataFrame(CSV.File(filename))
    if uppercase(names(df)[1]) in ["COLUMN1", "DATE", "DATES",
                                   "PERIOD", "PERIODS"]
        frequency = identify_period_frequency(uppercase(df[1, 1]))
    else
        frequency = Undated
    end
    
    ny = length(variables)

    if last == 0
        last = size(df, 1) - start + 1
    end
    nobs = last - start + 1
    Y = Matrix{Union{Missing,Float64}}(undef, ny, nobs)
    for (i, v) in enumerate(variables)
        Y[i, :] .= df[start:last, v]
    end
    return Y
end

