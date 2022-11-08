using Dynare

function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing
    Base.precompile(Tuple{typeof(parser),String,CommandLineOptions})   # time: 4.983775
    isdefined(Dynare, Symbol("#146#149")) &&
        Base.precompile(Tuple{getfield(Dynare, Symbol("#146#149")),Float64,Int64,Int64})   # time: 0.0036571
end

@dynare "test/models/example1/example1.mod"
@dynare "test/models/example2/example2.mod"
@dynare "test/models/example3/example3.mod"
@dynare "test/models/example3ss/example3ss.mod"
#@dynare "test/models/example3ss/example3ss_analytical.mod"
#@dynare "test/models/example3ss/example3ss_partial.mod"
@dynare "test/models/example3report/example3report.mod"
@dynare "test/models/cgg/cgg_ramsey.mod"
#@dynare "test/models/example1pf/example1pf.mod"