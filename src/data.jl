using AxisArrayTables
using ExtendedDates

function find_letter_in_period(period::AbstractString)
    for c in period
        if 'A' <= c <= 'z'
            return uppercase(c)
        end
    end
    return nothing
end

function identify_period_type(period::Union{AbstractString, Number})
    if typeof(period) <: Number
        if isinteger(period)
            return Undated
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

function periodparse(period::Union{AbstractString, Number})::ExtendedDates.PeriodsSinceEpoch
    period_type = identify_period_type(period)
    if period_type == YearSE
        return parse(period_type, period)
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
    elseif period_type == Undated
        return Int(period)
    else
        throw(ErrorException)
    end
end

function MyAxisArrayTable(filename)
    table = CSV.File(filename)
    cols = AxisArrayTables.Tables.columnnames(table)
    data = AxisArrayTables.Tables.matrix((;ntuple(i -> Symbol(i) => Tables.getcolumn(table, i), length(cols))...))
    for (icol, name) in enumerate(cols)
        if uppercase(String(name)) in ["DATE", "DATES", "PERIOD", "PERIODS", "TIME"]
            rows = []
            foreach(x -> push!(rows, periodparse(x)), data[:, icol])
            k = union(1:icol-1, icol+1:size(data,2))
            aat = AxisArrayTable(data[:, k], rows, cols[k])
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

function get_data(
    filename::String,
    variables::Vector{String};
    start::Int64 = 1,
    last::Int64 = 0,
)
    aat = MyAxisArrayTable(filename)
    ny = length(variables)
    if last == 0
        last = size(aat, 1) - start + 1
    end
    nobs = last - start + 1
    Y = Matrix{Union{Missing,Float64}}(undef, ny, nobs)
    for (i, v) in enumerate(variables)
        Y[i, :] .= Matrix(aat[:, Symbol(v)])[start:last]
    end
    return Y
end
