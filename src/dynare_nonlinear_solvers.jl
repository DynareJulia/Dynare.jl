struct DynareNonlinearSolverFailed <: Exception end
Base.showerror(io::IO, e::DynareNonlinearSolverFailed) = print("""
                                                                      Dynare couldn't solve the nonlinear system.
                                                                      Either there is no solution or the guess values
                                                                      are too far from the solution
                                                                       """)

function dynare_nonlinear_solvers(f!, j!, JJ, lb, ub, guess_values, exogenous, params;
                                  solve_algo = "NonlinearSolve",
                                  linear_solve_algo = nothing,
                                  method = TrustRegion(),
                                  show_trace = false,
                                  tolf = 1e-5,
                                  tolx = 1e-5,
                                  maxit = 100)
    solve_algo == "PATHSolver" && linear_solve_algo == "Pardiso" && error("Pardiso can't be used with PATHSolver")
    y0 = vec(guess_values)
    
    let res
        if solve_algo == "PATHSolver"
            (status, results, info) = solve_path!(f!, j!, JJ, lb, ub, y0)
        elseif solve_algo == "NonlinearSolve"
            fj = NonlinearFunction(f!, jac = j!, jac_prototype = JJ)
            prob = NonlinearProblem(fj, y0, params)
            r = zeros(length(prob.u0))
            prob.f(r, y0, params)
            res = NonlinearSolve.solve(prob, method, show_trace = Val(show_trace), abstol = tolf)
            if res.retcode == NonlinearSolve.ReturnCode.Success
                return res.u
            else
                @debug "Solving nonlinear system failed with $(res.retcode)"
                throw(DynareNonlinearSolverFailed())
            end
        elseif solve_algo == "NLsolve"
            if linear_solve_algo == "Pardiso"
                ls!(x, A, b) = linear_solver!(PardisoLS(), x, A, b)
                res = nlsolve(df, y0, method = method, show_trace = show_trace, ftol = tolf, xtol = tolx, iterations = maxit, linsolve = ls!)    
            else
                res = nlsolve(df, y0, method = method, show_trace = show_trace, ftol = tolf, xtol = tolx, iterations = maxit)    
            end
            if result.retcode == ReturnCode.Success
                return result.u
            else
                @debug "Solving nonlinear system failed with\n $result"
                throw(DynareNonlinearSolverFailed())
            end
        end
    end
end
