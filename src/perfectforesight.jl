function perfect_foresight_setup!(context, field)
    context.options["perfect_foresight_setup"] = Dict()
    copy!(context.options["perfect_foresight_setup"], field["options"])
    compute_perfect_foresight_setup(context)
end

function perfect_foresight_solver!(context, field)
    context.options["perfect_foresight_solver"] = Dict()
    copy!(context.options["perfect_foresight_solver"], field["options"])
    compute_perfect_foresight_solver(context)
end

function compute_prefect_foresight_setup(context); end
function compute_perfect_foresight_solver(context); end

