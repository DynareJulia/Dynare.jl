function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(parser),String,CommandLineOptions})   # time: 4.983775
    isdefined(Dynare, Symbol("#146#149")) &&
        Base.precompile(Tuple{getfield(Dynare, Symbol("#146#149")),Float64,Int64,Int64})   # time: 0.0036571
end
