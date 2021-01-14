using Artifacts

function dynare_preprocess(modfilename, args)
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
end

function run_dynare(modfilename, dynare_args)
    directory = dirname(modfilename)
    if length(directory) > 0
        current_directory = pwd()
        cd(directory)
    end

    dynare_preprocessor_path = joinpath(artifact"dynare-preprocessor", "dynare-preprocessor")
    run(`$dynare_preprocessor_path $dynare_args`)

    if length(directory) > 0
        cd(current_directory)
    end
end
