

function make_state_space!(A, B, Σy, g1_1, g1_2, context, ws)
    vA = view(A, :, context.models[1].i_bckwrd_b]
    vA .= g1_1
    B = g1_2
    extended_lyapd!(Σy, A, B, ws.lyapd_ws)
end
