using DynarePreprocessor_jll
using Pkg

function dynare_preprocess(modfilename::String, args::Vector{Any})
    dynare_args = [basename(modfilename), "language=julia", "json=compute"]
    offset = 0
    for a in args
        astring = string(a)
        if !occursin(r"^json=", astring)
            push!(dynare_args, astring)
        end
    end
    println(dynare_args)
    run_dynare(modfilename, dynare_args)
    println("")
end

function run_dynare(modfilename::String, dynare_args::Vector{String})
    @info "Dynare preprocessor version: $(module_version(DynarePreprocessor_jll))"
    directory = dirname(modfilename)
    if length(directory) > 0
        current_directory = pwd()
        cd(directory)
    end

    dynare_preprocessor_path = dynare_preprocessor()
    run(`$dynare_preprocessor_path $dynare_args`)

    if length(directory) > 0
        cd(current_directory)
    end
end
