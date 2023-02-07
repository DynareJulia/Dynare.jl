using FastLapackInterface

struct CalibSmootherOptions
    datafile::String
    first_obs::Int64
    last_obs::Int64
    function CalibSmootherOptions(options::Dict{String,Any})
        datafile = ""
        first_obs = 1
        last_obs = 0
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
    calib_smoother_core!(context, options)
end

function calib_smoother_core!(contex::Context, options::CalibSmootherOptions)
    symboltable = context.symboltable
    varobs = context.work.observed_variables
    varobs_ids =
        [symboltable[v].orderintype for v in varobs if is_endogenous(v, symboltable)]
    model = context.models[1]
    results = context.results.model_results[1]
    if (filename = options.datafile) != ""
        varnames = [v for v in varobs if is_endogenous(v, symboltable)]
        Yorig =
            get_data(filename, varnames, start = options.first_obs, last = options.last_obs)
    else
        error("calib_smoother needs a data file or a TimeDataFrame!")
    end
    @show context.work.model_has_trend[1]
    Y = Matrix{Union{Float64,Missing}}(undef, size(Yorig))
    if context.modfileinfo.has_trends
        remove_linear_trend!(
            Y,
            Yorig,
            results.trends.endogenous_steady_state[varobs_ids],
            results.trends.endogenous_linear_trend[varobs_ids],
        )
    else
        Y = Yorig
    end
    @show Y 
    statevar_ids = model.i_bkwrd_b
    kalman_statevar_ids = collect(1:model.endogenous_nbr)
    ns = length(kalman_statevar_ids)
    np = model.exogenous_nbr
    ny, nobs = size(Y)
    c = zeros(ny)
    k1 = findall(in(varobs_ids), kalman_statevar_ids)
    k2 = findall(in(statevar_ids), kalman_statevar_ids)
    Z = zeros(ny, ns)
    for i = 1:ny
        Z[i, varobs_ids[i]] = 1.0
    end
    H = zeros(ny, ny)
    d = zeros(ns)
    T = zeros(ns, ns)
    vg1 = view(
        context.results.model_results[1].linearrationalexpectations.g1_1,
        kalman_statevar_ids,
        :,
    )
    T[:, k2] .= vg1
    R = zeros(ns, np)
    vg2 = view(
        context.results.model_results[1].linearrationalexpectations.g1_2,
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
        context.results.model_results[1].linearrationalexpectations.endogenous_variance,
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
    if count(results.stationary_variables) == model.endogenous_nbr
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
            data_pattern,
        )
    else
        schur_ws = SchurWs(T)
        F = Schur(
            LAPACK.gees!(
                schur_ws,
                'V',
                T,
                select = (wr, wi) -> wr * wr + wi * wi > 1 - 1e-6,
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
            data_pattern,
        )
        alphah = F.Z * alphah
    end

    results.smoother["alphah"] = Matrix{Float64}(undef, ns, nobs)
    add_linear_trend!(
        results.smoother["alphah"],
        alphah,
        results.trends.endogenous_steady_state,
        results.trends.endogenous_linear_trend,
    )
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
