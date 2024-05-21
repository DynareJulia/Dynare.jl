using FastLapackInterface

struct CalibSmootherOptions
    datafile::String
    first_obs::Int64
    last_obs::Int64
    function CalibSmootherOptions(options::Dict{String,Any})
        datafile = ""
        first_obs = Undated(typemin(Int))
        last_obs = Undated(typemin(Int))
        for (k, v) in pairs(options)
            if k == "datafile"
                datafile = v::String
            elseif k == "first_obs"
                first_obs = v::Int64
            elseif k == "last_obs"
                last_obs = v::Int64
            end
        end
        new(datafile, first_obs, last_obs)
    end
end

function calib_smoother!(context::Context, field::Dict{String,Any})
    options = CalibSmootherOptions(field["options"])
    calibsmoother!(context=context, 
                    datafile = options.datafile,
                    first_obs = options.first_obs,
                    last_obs = options.last_obs)
end

"""
    calibsmoother!(; context::Context=context,
                     datafile::String = "",
                         data::AxisArrayTable = AxisArrayTable(Matrix{Float64}(undef, 0, 0), Undated[], Symbol[]),
                     first_obs::PeriodSinceEpoch = Undated(typemin(Int)),
                     last_obs::PeriodSinceEpoch = Undated(typemin(Int)),
                     nobs::Int = 0
                   )
Compute the smoothed values of the variables for an estimated model

#Keyword arguments
- `periods::Integer`: number of forecasted periods [required]
- `datafile::String`: file with the observations for the smoother
- `data::AxisArrayTable`: AxisArrayTable containing observed variables (can't be used at the same time as `datafile`)
- `first_obs::PeriodSinceEpoch`: first period used by smoother
                                 (default: first observation in the dataset)  
- `last_obs::PeriodSinceEpoch`: last period used by smoother
                                (default: last observation in the dataset)
- `nobs::Int`: number of observations (default: entire dataset)
"""
function calibsmoother!(; context=context,
                        datafile = "",
                        data::AxisArrayTable = AxisArrayTable(Matrix{Float64}(undef, 0, 0), Undated[], Symbol[]),
                        first_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
                        last_obs::PeriodsSinceEpoch = Undated(typemin(Int)),
                        nobs::Int = 0
                        )
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    has_trends = context.modfileinfo.has_trends
    endogenous_vars = get_endogenous(symboltable)
    exogenous_vars = get_exogenous(symboltable)
    varobs_ids =
        [symboltable[v].orderintype for v in varobs if is_endogenous(v, symboltable)]
    model = context.models[1]
    results = context.results.model_results[1]
    lre_results = results.linearrationalexpectations
    observations  = get_transposed_data!(context, datafile, data, varobs, first_obs, last_obs, nobs)
    Y = deepcopy(observations)
    detrend_data!(
        Y,
        context,
        context.work.observed_variables,
        dim = 1
    )

    periods = row_labels(context.work.data)
    statevar_ids = model.i_bkwrd_b
    kalman_statevar_ids = collect(1:model.endogenous_nbr)
    ns = length(kalman_statevar_ids)
    np = model.exogenous_nbr
    ny, nobs = size(Y)
    c = zeros(ny)
    k1 = findall(in(varobs_ids), kalman_statevar_ids)
    k2 = findall(in(statevar_ids), kalman_statevar_ids)
   
    Z = zeros(ny, ns)
#    Y = zeros(Union{Missing, Float64}, nobs, ny)
    endo_names = get_endogenous(context.symboltable)
    for i = 1:ny
        #        Y[:, i] .= Matrix(data[Symbol(varobs[i])])
        Z[i, varobs_ids[i]] = 1.0
    end
    H = zeros(ny, ny)
    d = zeros(ns)
    T = zeros(ns, ns)
    vg1 = view(
        lre_results.g1_1,
        kalman_statevar_ids,
        :,
    )
    T[:, k2] .= vg1
    R = zeros(ns, np)
    vg2 = view(
        lre_results.g1_2,
        kalman_statevar_ids,
        :,
    )
    R .= vg2
    Q = model.Sigma_e
    a0 = zeros(ns, nobs + 1)
    alphah = zeros(ns, nobs)
    att = zeros(ns, nobs)
    epsilonh = zeros(ny, nobs)
    etah = zeros(np, nobs)
    P = zeros(ns, ns, nobs + 1)
    Ptt = zeros(ns, ns, nobs + 1)
    vv = view(
        lre_results.endogenous_variance,
        kalman_statevar_ids,
        kalman_statevar_ids,
    )
    P[:, :, 1] .= vv
    Valpha = zeros(ns, ns, nobs)
    Vepsilon = zeros(ny, ny, nobs)
    Veta = zeros(np, np, nobs)
    start = 1
    last = nobs
    presample = 0
    data_pattern = Vector{Vector{Int64}}(undef, 0)
    for i = 1:nobs
        push!(data_pattern, findall(.!ismissing.(Y[:, i])))
    end
    if count(lre_results.stationary_variables) == model.endogenous_nbr
        kws = KalmanSmootherWs{Float64,Int64}(ny, ns, model.exogenous_nbr, nobs)
        kalman_smoother!(
            Y,
            c,
            Z,
            H,
            d,
            T,
            R,
            Q,
            a0,
            att,
            P,
            Ptt,
            alphah,
            epsilonh,
            etah,
            Valpha,
            Vepsilon,
            Veta,
            start,
            last,
            presample,
            kws,
            data_pattern
        )
    else
        schur_ws = SchurWs(T)
        F = Schur(
            LAPACK.gees!(
                schur_ws,
                'V',
                T,
                select = FastLapackInterface.ed,
                criterium = 1 - 1e-6,
            )...,
        )
        td = transpose(F.Z) * d
        tR = transpose(F.Z) * R
        tZ = Z * F.Z

        P = zeros(ns, ns, nobs + 1)
        k = count(abs.(F.values) .> 1 - 1e-6)
        vT = view(T, (k+1):ns, (k+1):ns)
        vP = view(P, (k+1):ns, (k+1):ns, 1)
        vtR = view(tR, (k+1):ns, :)
        k1 = ns - k
        lyapd_ws = LyapdWs(k1)
        LinearRationalExpectations.extended_lyapd_core!(
            vP,
            vT,
            vtR * Q * transpose(vtR),
            lyapd_ws,
        )
        Pinf = zeros(ns, ns, nobs + 1)
        for i = 1:k
            Pinf[i, i, 1] = 1.0
        end
        Pinftt = zeros(ns, ns, nobs + 1)
        kws = DiffuseKalmanSmootherWs{Float64,Int64}(ny, ns, model.exogenous_nbr, nobs)
        diffuse_kalman_smoother!(
            Y,
            c,
            tZ,
            H,
            td,
            T,
            tR,
            Q,
            a0,
            att,
            Pinf,
            Pinftt,
            P,
            Ptt,
            alphah,
            epsilonh,
            etah,
            Valpha,
            Vepsilon,
            Veta,
            start,
            last,
            presample,
            1e-8,
            kws,
            data_pattern
        )
        a0 = F.Z * a0
        alphah = F.Z * alphah
    end 

    endo_symb = [Symbol(v) for v in endogenous_vars]
    exo_symb = [Symbol(v) for v in exogenous_vars]
    smoother = copy(alphah)
    filter = copy(a0)
    steadystate = context.results.model_results[1].trends.endogenous_steady_state
    if has_trends
        add_linear_trend!(
            filter,
            a0,
            steadystate,
            results.trends.endogenous_linear_trend,
        )
        add_linear_trend!(
            smoother,
            alphah,
            steadystate,
            results.trends.endogenous_linear_trend,
        )
    else
        filter .+= steadystate
        smoother .+= steadystate
    end
    lastperiod = periods[end]
    T = typeof(lastperiod)
    if T <: Dates.UTInstant
        periods1 = vcat(periods, T(lastperiod.periods.value + 1))
    else
        periods1 = vcat(periods, lastperiod + 1)
    end 
    results.filter = AxisArrayTable(transpose(filter), 
                                    periods1, 
                                    endo_symb)
    results.smoother = AxisArrayTable(transpose(vcat(smoother, etah)), 
                                      periods, 
                                      vcat(endo_symb, exo_symb))
    return nothing
end

function make_state_space!(
    A::VecOrMat{Float64},
    B::VecOrMat{Float64},
    Σy::Matrix{Float64},
    g1_1::AbstractVecOrMat{Float64},
    g1_2::AbstractVecOrMat{Float64},
    context::Context,
    ws::LinearRationalExpectationsWs,
)
    vA = view(A, :, context.models[1].i_bckwrd_b)
    vA .= g1_1
    B = g1_2
    extended_lyapd!(Σy, A, B, ws.lyapd_ws)
end
