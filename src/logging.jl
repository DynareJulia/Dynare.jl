function filter_logging_messages(log_args)
    if log_args.level == Logging.Debug && ! occursin("Dynare", get(ENV, "JULIA_DEBUG", ""))
        return false
    end
    return true
end

function set_logging(modname)
    FT = TeeLogger(
        ActiveFilteredLogger(
            filter_logging_messages,
            FormatLogger() do io, args
            println(io, args.message)
            end
        ),
        ActiveFilteredLogger(
            filter_logging_messages,
            FormatLogger("$(modname).log") do io, args
            println(io, args.message)
            end
        )
    )
    global_logger(FT)
end

