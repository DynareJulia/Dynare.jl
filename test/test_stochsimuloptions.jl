module test_stochsimuloptions
struct StochSimulOptions
    dr_algo::String
    order::Int64
    periods::Int64
    function StochSimulOptions(options::Dict{String, Any})
        ks = keys(options)
        if "dr_cyclic_reduction" in ks && options["dr_cyclic_reduction"]::Bool
            dr_algo = "CR"
        else
            dr_algo = "GS"
        end
        if "order" in ks
            order = options["order"]::Int64
        else
            order = 1
        end
        if "periods" in ks
            periods = options["periods"]::Int64
        else
            periods = 0
        end
        new(dr_algo, order, periods)
    end
end
end
