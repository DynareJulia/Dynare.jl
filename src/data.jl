function identify_period_frequency(period)::Periods.Frequency
    period = uppercase(period)
    if 'Y' in period
        frequency = Periods.Year
    elseif 'A' in period
        frequency = Periods.Year
    elseif 'S' in period
        frequency = Period.Semester
    elseif 'H' in period
        frequency = Periods.Semester
    elseif 'Q' in period
        frequency = Periods.Quarter
    elseif 'M' in period
        frequency = Periods.Month
    elseif 'W' in period
        frequency = Periods.Week
    elseif 'D' in period
        frequency = Periods.Day
    elseif (cdash = count("-", period)) == 2
        frequency = Periods.Day
    elseif cdash == 1
        frequency = Periods.Month
    else
        throw(ErrorException)
    end
end

function get_data(filename::String, variables::Vector{String}, options::Dict{String, Any})
    df = DataFrame(CSV.File(filename))
    if uppercase(names(df)[1]) in ["COLUMN1", "DATE", "DATES",
                                   "PERIOD", "PERIODS"]
        frequency = identify_period_frequency(uppercase(df[1, 1]))
    else
        frequency = Periods.Undated
    end
    
    ny = length(variables)

    start = get(options, "first_obs", 1)
    last = get(options, "last_obs", size(df, 1) - start + 1) 
    nobs = last - start + 1
    Y = Matrix{Union{Missing,Float64}}(undef, ny, nobs)
    for (i, v) in enumerate(variables)
        Y[i, :] .= df[start:last, v]
    end
    return Y
end

