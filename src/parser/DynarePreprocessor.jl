DYNARE_BINARY = "/data/projects/dynare/git/preprocessor/src/dynare_m"

function dynare_preprocess(modfilename, args)
    dynare_args = [basename(modfilename), "language=julia", "output=third", "json=compute"]
    offset = 0
    for a in args
        if occursin(r"^output=", a)
            deleteat!(dynare_args, 3)
            offset = 1
        elseif occursin(r"^json=", a)
            deleteat!(dynare_args, 4 - offset)
        end
    end
    append!(dynare_args, args)
    println(dynare_args)
    current_directory = pwd()
    directory = dirname(modfilename)
    cd(directory)
    run(`$DYNARE_BINARY $dynare_args`)
    cd(current_directory)
end
